"""Adaptador OpenGL moderno basado en shaders y buffers Qt."""

from __future__ import annotations

from array import array

from PyQt6 import sip
from PyQt6.QtGui import QMatrix4x4, QOpenGLContext, QVector3D
from PyQt6.QtOpenGL import (
    QOpenGLBuffer,
    QOpenGLFunctions_4_1_Core,
    QOpenGLShader,
    QOpenGLShaderProgram,
    QOpenGLVertexArrayObject,
)

from .camera import CameraProjection, CameraState
from .mesh import build_scene_mesh

GL_COLOR_BUFFER_BIT = 0x00004000
GL_DEPTH_BUFFER_BIT = 0x00000100
GL_DEPTH_TEST = 0x00000B71
GL_LESS = 0x00000201
GL_TRIANGLES = 0x0004
GL_FLOAT = 0x1406

_VERTEX_SHADER = """#version 410 core
layout(location = 0) in vec3 a_position;
layout(location = 1) in vec3 a_normal;
layout(location = 2) in vec3 a_color;
uniform mat4 u_mvp;
out vec3 v_normal;
out vec3 v_color;
void main() {
    v_normal = a_normal;
    v_color = a_color;
    gl_Position = u_mvp * vec4(a_position, 1.0);
}
"""

_FRAGMENT_SHADER = """#version 410 core
in vec3 v_normal;
in vec3 v_color;
out vec4 frag_color;
void main() {
    vec3 normal = normalize(v_normal);
    vec3 light = normalize(vec3(-0.35, 0.75, 0.55));
    float diffuse = 0.28 + 0.72 * max(dot(normal, light), 0.0);
    frag_color = vec4(v_color * diffuse, 1.0);
}
"""


class MoleculeRenderer:
    """Sube y dibuja una malla de escena sin conocer el modelo químico."""

    def __init__(self) -> None:
        self._gl: QOpenGLFunctions_4_1_Core | None = None
        self._program: QOpenGLShaderProgram | None = None
        self._vbo: QOpenGLBuffer | None = None
        self._vao: QOpenGLVertexArrayObject | None = None
        self._vertex_count = 0
        self._initialized = False
        self._context: QOpenGLContext | None = None

    @property
    def initialized(self) -> bool:
        return self._initialized

    def initialize(self) -> tuple[bool, str]:
        context = QOpenGLContext.currentContext()
        if context is None:
            return False, "No hay un contexto OpenGL actual para inicializar el renderer"
        if self._initialized and self._context == context:
            return True, "OpenGL 4.1 core inicializado"
        self._context = context
        try:
            gl = QOpenGLFunctions_4_1_Core()
            if not gl.initializeOpenGLFunctions():
                self.cleanup()
                return False, "Qt no pudo inicializar funciones OpenGL 4.1 core"
            program = QOpenGLShaderProgram()
            self._program = program
            if not program.addShaderFromSourceCode(QOpenGLShader.ShaderTypeBit.Vertex, _VERTEX_SHADER):
                self.cleanup()
                return False, program.log() or "No se pudo compilar el vertex shader"
            if not program.addShaderFromSourceCode(QOpenGLShader.ShaderTypeBit.Fragment, _FRAGMENT_SHADER):
                self.cleanup()
                return False, program.log() or "No se pudo compilar el fragment shader"
            if not program.link():
                self.cleanup()
                return False, program.log() or "No se pudo enlazar el shader program"
            vbo = QOpenGLBuffer(QOpenGLBuffer.Type.VertexBuffer)
            self._vbo = vbo
            vao = QOpenGLVertexArrayObject()
            self._vao = vao
            if not vbo.create() or not vao.create():
                self.cleanup()
                return False, "No se pudieron crear los buffers OpenGL"
            self._gl = gl
            self._initialized = True
            gl.glEnable(GL_DEPTH_TEST)
            gl.glDepthFunc(GL_LESS)
            return True, "OpenGL 4.1 core inicializado"
        except Exception as exc:
            self.cleanup()
            return False, f"OpenGL no disponible: {exc}"

    def upload(self, scene, selected_atom_id=None, selected_bond_id=None) -> bool:
        if not self._initialized or self._vbo is None or self._vao is None or self._program is None:
            return False
        if QOpenGLContext.currentContext() != self._context:
            return False
        mesh = build_scene_mesh(scene, selected_atom_id, selected_bond_id)
        self._vertex_count = mesh.vertex_count
        self._vao.bind()
        self._vbo.bind()
        payload = array("f", mesh.vertices).tobytes()
        self._vbo.allocate(sip.voidptr(payload), len(payload))
        stride = 9 * 4
        self._program.enableAttributeArray(0)
        self._program.enableAttributeArray(1)
        self._program.enableAttributeArray(2)
        self._program.setAttributeBuffer(0, GL_FLOAT, 0, 3, stride)
        self._program.setAttributeBuffer(1, GL_FLOAT, 3 * 4, 3, stride)
        self._program.setAttributeBuffer(2, GL_FLOAT, 6 * 4, 3, stride)
        self._vbo.release()
        self._vao.release()
        return True

    def render(self, camera: CameraState, width: int, height: int) -> None:
        if not self._initialized or self._gl is None or self._program is None or self._vao is None:
            return
        if QOpenGLContext.currentContext() != self._context:
            return
        if width <= 0 or height <= 0:
            return
        gl = self._gl
        projection = QMatrix4x4()
        aspect = float(width) / float(height)
        near, far = camera.clip_planes()
        if camera.projection == CameraProjection.ORTHOGRAPHIC:
            half_height = camera.ortho_scale * 0.5
            projection.ortho(-half_height * aspect, half_height * aspect, -half_height, half_height, near, far)
        else:
            projection.perspective(camera.fov_y_degrees, aspect, near, far)
        view = QMatrix4x4()
        forward, _right, up = camera.basis
        eye = camera.eye
        view.lookAt(
            QVector3D(*eye),
            QVector3D(eye[0] + forward[0], eye[1] + forward[1], eye[2] + forward[2]),
            QVector3D(*up),
        )
        self._program.bind()
        self._program.setUniformValue("u_mvp", projection * view)
        self._vao.bind()
        gl.glDrawArrays(GL_TRIANGLES, 0, self._vertex_count)
        self._vao.release()
        self._program.release()

    def clear(self) -> None:
        self._vertex_count = 0

    def cleanup(self) -> None:
        """Release Qt GPU wrappers while their owning context is current."""
        context_is_current = QOpenGLContext.currentContext() == self._context
        cleanup_error: Exception | None = None
        try:
            if context_is_current:
                for resource, destroy in (
                    (self._vbo, lambda item: item.destroy() if item.isCreated() else None),
                    (self._vao, lambda item: item.destroy() if item.isCreated() else None),
                    (self._program, lambda item: (item.release(), item.removeAllShaders())),
                ):
                    if resource is None:
                        continue
                    try:
                        destroy(resource)
                    except Exception as exc:
                        if cleanup_error is None:
                            cleanup_error = exc
        finally:
            self._gl = None
            self._program = None
            self._vbo = None
            self._vao = None
            self._vertex_count = 0
            self._initialized = False
            self._context = None
        if cleanup_error is not None:
            raise cleanup_error
