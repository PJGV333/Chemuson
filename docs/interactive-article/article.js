(() => {
  'use strict';

  const D = window.CHEMUSON_DATA;
  const MODULES = window.CHEMUSON_MODULES;
  const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const byId = new Map(D.atoms.map(atom => [atom.id, atom]));
  let selectedAtomId = 1;
  const cleanCoords = D.clean2d.after;
  const SVG_NS = 'http://www.w3.org/2000/svg';
  const fmt = (value, digits = 1) => Number.isFinite(Number(value)) ? Number(value).toFixed(digits) : '—';
  const esc = value => String(value ?? '').replace(/[&<>"']/g, ch => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[ch]);

  function coord(atom, t = 0) {
    const target = cleanCoords[String(atom.id)] || [atom.x, atom.y];
    return [atom.x + (target[0] - atom.x) * t - 140, atom.y + (target[1] - atom.y) * t - 365];
  }

  function bondMarkup(bond, t, graphMode, reveal) {
    const a = byId.get(bond.a1), b = byId.get(bond.a2);
    if (!a || !b) return '';
    const [x1, y1] = coord(a, t), [x2, y2] = coord(b, t);
    const dx = x2 - x1, dy = y2 - y1, length = Math.hypot(dx, dy) || 1;
    const ox = -dy / length * 2.3, oy = dx / length * 2.3;
    const order = graphMode ? 1 : Number(bond.order || 1);
    const delay = reveal ? ` style="--delay:${Math.min(bond.a1, bond.a2) * 58}ms"` : '';
    const revealClass = reveal ? ' reveal' : '';
    let lines = `<line class="bond-line${revealClass}"${delay} x1="${x1.toFixed(2)}" y1="${y1.toFixed(2)}" x2="${x2.toFixed(2)}" y2="${y2.toFixed(2)}"/>`;
    if (order >= 2) {
      lines += `<line class="bond-line secondary${revealClass}"${delay} x1="${(x1 + ox).toFixed(2)}" y1="${(y1 + oy).toFixed(2)}" x2="${(x2 + ox).toFixed(2)}" y2="${(y2 + oy).toFixed(2)}"/>`;
    }
    if (order >= 3 && !graphMode) {
      lines += `<line class="bond-line secondary${revealClass}"${delay} x1="${(x1 - ox).toFixed(2)}" y1="${(y1 - oy).toFixed(2)}" x2="${(x2 - ox).toFixed(2)}" y2="${(y2 - oy).toFixed(2)}"/>`;
    }
    return lines;
  }

  function renderMolecule(container, options = {}) {
    const {
      t = 0, graph = false, interactive = false, selected = null,
      reveal = false, layer = 'structure', numbered = false,
      ariaLabel = 'Estructura molecular representada desde datos del fixture'
    } = options;
    const ringAtoms = new Set((D.layers.motifs.find(m => m.kind === 'ring') || { atoms: [] }).atoms);
    const macrocycle = new Set((D.layers.blocks.find(b => b.kind === 'macrocycle') || { atoms: [] }).atoms);
    const predictedAtoms = new Set(D.atoms.filter(atom => atom.h1 || atom.c13).map(atom => atom.id));
    const bonds = D.bonds.map(bond => bondMarkup(bond, t, graph, reveal)).join('');
    let vectors = '';
    if (layer === 'geometry') {
      vectors = D.atoms.filter(atom => [1, 5, 10, 18].includes(atom.id)).map(atom => {
        const [x1, y1] = coord(atom, 0), [x2, y2] = coord(atom, 1);
        return `<line class="bond-line interaction" x1="${x1.toFixed(2)}" y1="${y1.toFixed(2)}" x2="${x2.toFixed(2)}" y2="${y2.toFixed(2)}"/>`;
      }).join('');
    }
    let overlayLabel = '';
    if (layer === 'naming') overlayLabel = '<text class="atom-label" x="488" y="22">IUPAC-lite · N/D</text>';
    if (layer === 'spectra') overlayLabel = '<text class="atom-label" x="456" y="22">heuristic-v1 · estimación</text>';
    const atomMarkup = D.atoms.map(atom => {
      const [x, y] = coord(atom, t);
      const pointR = interactive ? 4.7 : 4.1;
      let fill = 'var(--accent)';
      if (layer === 'validation') fill = atom.error ? 'var(--orange)' : 'var(--green)';
      else if (layer === 'properties') fill = atom.implicit_h > 1 ? 'var(--orange)' : 'var(--accent)';
      else if (layer === 'motifs') fill = ringAtoms.has(atom.id) ? 'var(--accent)' : 'var(--ink-3)';
      else if (layer === 'blocks') fill = macrocycle.has(atom.id) ? 'var(--blue)' : 'var(--ink-3)';
      else if (layer === 'spectra') fill = predictedAtoms.has(atom.id) ? 'var(--orange)' : 'var(--ink-3)';
      const revealClass = reveal ? ' reveal' : '';
      const delay = reveal ? ` style="--delay:${atom.id * 58}ms"` : '';
      const label = `${atom.element} ${atom.id}; valencia de enlace ${fmt(atom.bond_sum)}; H implícitos ${atom.implicit_h}`;
      if (interactive) {
        return `<g class="atom-node${atom.id === selected ? ' is-selected' : ''}${revealClass}"${delay} data-atom="${atom.id}" role="button" tabindex="0" aria-label="${esc(label)}" aria-pressed="${atom.id === selected}"><circle class="atom-hit" cx="${x.toFixed(2)}" cy="${y.toFixed(2)}" r="13"/><circle class="atom-point" cx="${x.toFixed(2)}" cy="${y.toFixed(2)}" r="${pointR}" style="fill:${fill}"/>${numbered || graph ? `<text class="atom-number" x="${(x + 6).toFixed(2)}" y="${(y - 7).toFixed(2)}">${atom.id}</text>` : ''}${atom.element !== 'C' ? `<text class="atom-label" x="${(x + 7).toFixed(2)}" y="${(y - 6).toFixed(2)}">${esc(atom.element)}</text>` : ''}</g>`;
      }
      return `<g class="atom-node${revealClass}"${delay} aria-hidden="true"><circle class="atom-point" cx="${x.toFixed(2)}" cy="${y.toFixed(2)}" r="${pointR}" style="fill:${fill}"/>${numbered || graph ? `<text class="atom-number" x="${(x + 6).toFixed(2)}" y="${(y - 7).toFixed(2)}">${atom.id}</text>` : ''}${atom.element !== 'C' ? `<text class="atom-label" x="${(x + 7).toFixed(2)}" y="${(y - 6).toFixed(2)}">${esc(atom.element)}</text>` : ''}</g>`;
    }).join('');
    const svgRole = interactive ? 'group' : 'img';
    const html = `<svg xmlns="${SVG_NS}" viewBox="0 0 620 180" role="${svgRole}" aria-label="${esc(ariaLabel)}" focusable="false"><g>${vectors}${bonds}${atomMarkup}</g>${overlayLabel}</svg>`;
    container.innerHTML = html;
    container.setAttribute('role', interactive ? 'group' : 'img');
    container.setAttribute('aria-label', ariaLabel);
    if (interactive) {
      container.querySelectorAll('[data-atom]').forEach(node => {
        const activate = () => selectAtom(Number(node.dataset.atom), container);
        node.addEventListener('click', activate);
        node.addEventListener('keydown', event => {
          if (event.key === 'Enter' || event.key === ' ') {
            event.preventDefault();
            activate();
          }
        });
      });
    }
  }

  function selectAtom(id, container = document.getElementById('stage-molecule')) {
    const atom = byId.get(id);
    if (!atom) return;
    selectedAtomId = id;
    if (container) {
      container.querySelectorAll('[data-atom]').forEach(node => {
        const active = Number(node.dataset.atom) === id;
        node.classList.toggle('is-selected', active);
        node.setAttribute('aria-pressed', String(active));
      });
    }
    const idNode = document.getElementById('atom-id');
    const summary = document.getElementById('atom-summary');
    const props = document.getElementById('atom-properties');
    if (!idNode || !summary || !props) return;
    idNode.textContent = `· ${atom.id}`;
    summary.textContent = `${atom.element} en el fixture; ${atom.motifs.length} motivos y ${atom.blocks.length} bloques asociados en los atributos serializados.`;
    const c13 = atom.c13 ? `${fmt(atom.c13[0])} ppm · ${atom.c13[1]} (${fmt(atom.c13[2] * 100, 0)}%)` : 'Sin señal reportada';
    const h1 = atom.h1 ? `${fmt(atom.h1[0])} ppm · ${atom.h1[2]} (${fmt(atom.h1[3] * 100, 0)}%)` : 'Sin señal reportada';
    const items = [
      ['Elemento', atom.element],
      ['Suma de órdenes', fmt(atom.bond_sum)],
      ['H implícitos', atom.implicit_h],
      ['Valencias permitidas', atom.allowed.join(', ')],
      ['Motivos', atom.motifs.join(', ') || '—'],
      ['Bloques', atom.blocks.join(', ') || '—'],
      ['¹H heurístico', h1],
      ['¹³C heurístico', c13],
      ['Error de valencia', atom.error ? 'sí' : 'no']
    ];
    props.innerHTML = items.map(([key, value]) => `<div><dt>${esc(key)}</dt><dd>${esc(value)}</dd></div>`).join('');
  }

  function setupHero() {
    const target = document.getElementById('hero-molecule');
    if (!target) return;
    renderMolecule(target, { t: 1, reveal: !reducedMotion, ariaLabel: 'Fixture de regresión de veinte carbonos, dibujado a partir de datos moleculares reales' });
  }

  function setupStages() {
    const target = document.getElementById('stage-molecule');
    const buttons = [...document.querySelectorAll('[data-stage]')];
    let stage = 'draw';
    const descriptions = {
      draw: ['Representación almacenada', 'La geometría de entrada conserva las coordenadas originales del fixture. Explora la figura por átomo o selecciona una etapa.'],
      graph: ['Modelo de conectividad', 'Cada vértice es un átomo y cada arista un enlace. Esta vista abstrae órdenes de enlace y geometría; no cambia el MolGraph.'],
      clean: ['Candidato seleccionado', 'Se muestran las coordenadas del candidato ganador de esta ejecución: internal_templates. La conectividad permanece igual.']
    };
    function draw() {
      const t = stage === 'clean' ? 1 : 0;
      const graph = stage === 'graph';
      renderMolecule(target, { t, graph, interactive: true, selected: selectedAtomId, numbered: graph, ariaLabel: graph ? 'Grafo de átomos y enlaces del fixture; selecciona un nodo para inspeccionar sus datos' : 'Estructura del fixture con átomos y enlaces interactivos' });
      document.getElementById('stage-kicker').textContent = descriptions[stage][0];
      document.getElementById('stage-description').textContent = descriptions[stage][1];
      document.getElementById('stage-status').textContent = `20 átomos · 21 enlaces · ${stage === 'clean' ? 'candidato aplicado' : stage === 'graph' ? 'vista abstracta' : 'coordenadas iniciales'}`;
      buttons.forEach(button => {
        const active = button.dataset.stage === stage;
        button.classList.toggle('is-active', active);
        button.setAttribute('aria-pressed', String(active));
      });
    }
    buttons.forEach(button => button.addEventListener('click', () => { stage = button.dataset.stage; draw(); }));
    draw();
    selectAtom(1, target);
  }

  function setupGeometrySlider() {
    const slider = document.getElementById('geometry-slider');
    const target = document.getElementById('clean-molecule');
    const output = document.getElementById('geometry-output');
    if (!slider || !target) return;
    const report = D.clean2d.report;
    function update() {
      const t = Number(slider.value) / 100;
      renderMolecule(target, { t, ariaLabel: `Geometría interpolada, ${(t * 100).toFixed(0)} por ciento hacia el candidato Clean2D` });
      output.value = t < 0.03 ? 'Entrada' : t > 0.97 ? 'Clean2D' : `${Math.round(t * 100)}%`;
      output.textContent = output.value;
      document.getElementById('metric-before').textContent = fmt(report.mean_before + (report.mean_after - report.mean_before) * t, 2);
      document.getElementById('metric-after').textContent = fmt(report.mean_after, 2);
    }
    slider.addEventListener('input', update);
    update();
  }

  const layerInfo = {
    structure: {
      status: 'MODELO', kind: 'model', title: 'Conectividad y atributos',
      description: `${D.atoms.length} átomos, ${D.bonds.length} enlaces. El grafo conserva conectividad y atributos moleculares; esta geometría viene del fixture de regresión.`,
      chips: ['20 átomos', '21 enlaces', '1 componente', '7 motivos', '5 bloques'],
      t: 1, visual: 'structure'
    },
    geometry: {
      status: 'CLEAN2D · APLICADO', kind: 'stable', title: 'Coordenadas de depiction',
      description: `Ganó ${D.clean2d.winner_source}. Longitud media ${fmt(D.clean2d.report.mean_before, 1)} → ${fmt(D.clean2d.report.mean_after, 1)}; distancia no enlazada mínima ${fmt(D.clean2d.report.min_nonbonded_after, 2)}.`,
      chips: ['0 cruces nuevos', `desplazamiento medio ${fmt(D.clean2d.report.mean_displacement, 1)}`, 'reporte: pasó'],
      t: 1, visual: 'geometry'
    },
    validation: {
      status: 'VALIDACIÓN', kind: 'stable', title: 'Consistencia local',
      description: 'En esta ejecución del fixture, el validador de valencias informó cero problemas; las marcas verdes representan átomos sin error reportado en los datos capturados.',
      chips: ['0 problemas', 'valencia por átomo', 'no es revisión experimental'],
      t: 0, visual: 'validation'
    },
    properties: {
      status: 'CHEMCALC', kind: 'stable', title: 'Propiedades calculadas',
      description: 'La fórmula y la masa provienen del cálculo auxiliar sobre el grafo del fixture; no identifican por sí solas un compuesto.',
      chips: [`${Object.entries(D.properties.formula).map(([element, n]) => `${element}${n}`).join('')}`, `${fmt(D.properties.mw, 4)} u`, 'masa promedio'],
      t: 1, visual: 'properties'
    },
    naming: {
      status: 'IUPAC-LITE · MVP', kind: 'mvp', title: 'Nomenclatura con degradación',
      description: `Para este fixture el resultado es “${D.properties.iupac}”. N/D expresa que este caso queda fuera de la cobertura del servicio; no implica que la estructura sea inválida.`,
      chips: ['cobertura parcial', 'no es Blue Book completo', 'fallo seguro'],
      t: 1, visual: 'naming'
    },
    spectra: {
      status: `${D.spectra.source.toUpperCase()} · HEURÍSTICO`, kind: 'heuristic', title: 'Estimaciones espectrales',
      description: `Confianza global ${fmt(D.spectra.confidence * 100, 0)}%. El espectro de masas mostrado es predicho, no una medición.`,
      chips: D.spectra.ms.map(peak => `${fmt(peak[0], 4)} m/z · ${fmt(peak[1], 0)}% · ${peak[2]}`),
      t: 1, visual: 'spectra'
    }
  };

  function setupLayers() {
    const target = document.getElementById('layer-molecule');
    const buttons = [...document.querySelectorAll('[data-layer]')];
    function update(key) {
      const info = layerInfo[key] || layerInfo.structure;
      renderMolecule(target, { t: info.t, layer: info.visual, ariaLabel: `Fixture molecular con capa ${info.title}` });
      const status = document.getElementById('layer-status');
      status.textContent = info.status;
      status.dataset.kind = info.kind;
      document.getElementById('layer-title').textContent = info.title;
      document.getElementById('layer-description').textContent = info.description;
      document.getElementById('layer-data').innerHTML = info.chips.map(chip => `<span class="data-chip">${esc(chip)}</span>`).join('');
      buttons.forEach(button => {
        const active = button.dataset.layer === key;
        button.classList.toggle('is-active', active);
        button.setAttribute('aria-pressed', String(active));
      });
    }
    buttons.forEach(button => button.addEventListener('click', () => update(button.dataset.layer)));
    update('structure');
  }

  function setupPersistenceFormats() {
    const diagram = document.getElementById('transport-diagram');
    const buttons = [...document.querySelectorAll('[data-format]')];
    const formats = {
      cmsn: {
        status: 'DOCUMENTO NATIVO', adapter: 'PersistenceManager', file: 'proyecto .cmsn', scope: 'MolGraph + datos del canvas',
        text: 'Serializa el modelo químico y los datos gráficos del canvas para reconstruir la sesión; no es solo una imagen exportada.'
      },
      smiles: {
        status: 'REPRESENTACIÓN ESTRUCTURAL', adapter: 'ChemIO · SMILES', file: 'cadena SMILES', scope: 'conectividad molecular',
        text: 'Intercambia una representación lineal de la estructura. No conserva anotaciones gráficas del canvas ni reemplaza el documento de proyecto.'
      },
      molfile: {
        status: 'INTERCAMBIO 2D', adapter: 'ChemIO · Molfile / SDF', file: '.mol / .sdf', scope: 'átomos · enlaces · coords 2D',
        text: 'Molfile conserva estructura y coordenadas 2D. SDF se soporta de forma básica/parcial; propiedades avanzadas y algunos extras visuales pueden degradarse.'
      },
      cml: {
        status: 'SUBCONJUNTO CML INICIAL', adapter: 'chemio.cml_io', file: 'documento .cml', scope: 'CML básico + extensiones ChemUSON',
        text: 'El importador/exportador maneja un subconjunto semántico de CML con extensiones propias. Flechas, textos y brackets del dibujo aún no se exportan como CML semántico.'
      },
      xyz: {
        status: 'COMP CHEM · 3D', adapter: 'exportador de geometría 3D', file: 'coordenadas .xyz', scope: 'coords 3D · sin conectividad',
        text: 'Exporta las coordenadas 3D activas, pero XYZ no conserva por sí solo la conectividad química ni el estado 2D del canvas.'
      }
    };
    function update(key) {
      const item = formats[key] || formats.cmsn;
      diagram.dataset.format = key;
      document.getElementById('transport-adapter').textContent = item.adapter;
      document.getElementById('transport-file').textContent = item.file;
      document.getElementById('transport-scope').textContent = item.scope;
      document.getElementById('format-status').textContent = item.status;
      document.getElementById('format-description').textContent = item.text;
      diagram.querySelectorAll('[data-transport-node]').forEach(node => {
        node.classList.toggle('is-active', node.dataset.transportNode === 'model' || node.dataset.transportNode === 'adapter' || node.dataset.transportNode === 'file');
      });
      buttons.forEach(button => {
        const active = button.dataset.format === key;
        button.classList.toggle('is-active', active);
        button.setAttribute('aria-pressed', String(active));
      });
    }
    buttons.forEach(button => button.addEventListener('click', () => update(button.dataset.format)));
    update('cmsn');
  }

  function setupArchitecture() {
    const svg = document.getElementById('module-map');
    if (!svg || !Array.isArray(MODULES)) return;
    const groups = [
      { label: 'MODELO', ids: ['M00'], x: 20, rows: 1 },
      { label: 'DOMINIO', ids: MODULES.map(m => m.id).filter(id => /^M0[1-7]$/.test(id)), x: 235, rows: 7 },
      { label: 'GUI', ids: ['M08', 'M09', 'M10', 'M11', 'M12', 'M13', 'M20'], x: 455, rows: 7 },
      { label: 'SERVICIOS', ids: ['M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M21', 'M22'], x: 675, rows: 8 }
    ];
    const positions = new Map();
    groups.forEach(group => group.ids.forEach((id, index) => positions.set(id, { x: group.x, y: group.ids.length === 1 ? 238 : 40 + index * 54 })));
    const items = new Map(MODULES.map(module => [module.id, module]));
    svg.innerHTML = groups.map(group => `<text x="${group.x}" y="19" fill="var(--ink-3)" font-family="var(--ui)" font-size="9" font-weight="700" letter-spacing="1.1">${group.label}</text>`).join('') + '<g id="map-edges"></g><g id="map-nodes"></g>';
    const nodes = document.getElementById('map-nodes');
    MODULES.forEach(module => {
      const pos = positions.get(module.id) || { x: 675, y: 40 };
      const name = module.name;
      const label = `${module.id} · ${name}`;
      const node = document.createElementNS(SVG_NS, 'g');
      node.setAttribute('class', 'map-node');
      node.setAttribute('data-module', module.id);
      node.setAttribute('role', 'button');
      node.setAttribute('tabindex', '0');
      node.setAttribute('aria-label', `${module.id}: ${module.title}. ${module.deps.length} dependencias directas.`);
      node.innerHTML = `<rect x="${pos.x}" y="${pos.y}" width="190" height="36"/><text class="node-id" x="${pos.x + 8}" y="${pos.y + 14}">${esc(module.id)}</text><text class="node-name" x="${pos.x + 48}" y="${pos.y + 15}">${esc(label.slice(module.id.length + 3, module.id.length + 23))}</text>`;
      node.addEventListener('click', () => selectModule(module.id));
      node.addEventListener('keydown', event => {
        if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); selectModule(module.id); }
        if (event.key === 'ArrowRight' || event.key === 'ArrowDown' || event.key === 'ArrowLeft' || event.key === 'ArrowUp') {
          event.preventDefault();
          const list = [...nodes.querySelectorAll('.map-node')];
          const at = list.indexOf(node);
          const step = (event.key === 'ArrowRight' || event.key === 'ArrowDown') ? 1 : -1;
          list[(at + step + list.length) % list.length].focus();
        }
      });
      nodes.appendChild(node);
    });
    let activeId = 'M02';
    function selectModule(id) {
      const module = items.get(id);
      if (!module) return;
      activeId = id;
      const deps = new Set(module.deps);
      document.querySelectorAll('.map-node').forEach(node => {
        const nodeId = node.dataset.module;
        node.classList.toggle('is-selected', nodeId === id);
        node.classList.toggle('is-dependency', deps.has(nodeId));
        node.classList.toggle('is-unrelated', nodeId !== id && !deps.has(nodeId));
        node.setAttribute('aria-pressed', String(nodeId === id));
      });
      const edges = document.getElementById('map-edges');
      edges.innerHTML = module.deps.map(depId => {
        const from = positions.get(id), to = positions.get(depId);
        if (!from || !to) return '';
        const sx = from.x + 95, sy = from.y + 18, ex = to.x + 95, ey = to.y + 18;
        const bend = (sx + ex) / 2 + (sx < ex ? 16 : -16);
        return `<path class="dep-edge" d="M${sx} ${sy} Q${bend} ${(sy + ey) / 2} ${ex} ${ey}"/>`;
      }).join('');
      document.getElementById('module-kicker').textContent = `${module.id} · ${module.name.toUpperCase()}`;
      document.getElementById('module-title').textContent = module.title;
      document.getElementById('module-responsibility').textContent = module.responsibility;
      const status = document.getElementById('module-status');
      status.textContent = module.status;
      status.className = `status-pill ${module.status.toLowerCase().replace(/[^a-z]/g, '')}`;
      const risk = document.getElementById('module-risk');
      risk.textContent = `riesgo ${module.risk}`;
      risk.className = `risk-pill ${module.risk}`;
      document.getElementById('module-deps').innerHTML = module.deps.length ? module.deps.map(dep => `<span>${esc(dep)}</span>`).join('') : '<span>sin dependencias declaradas</span>';
      document.getElementById('module-api').innerHTML = (module.public_api.length ? module.public_api : ['(sin API pública listada)']).map(api => `<li>${esc(api)}</li>`).join('');
    }
    selectModule(activeId);
  }

  function setupWorker() {
    const diagram = document.getElementById('worker-diagram');
    const buttons = [...document.querySelectorAll('[data-worker-step]')];
    const text = [
      'La ventana coordina una operación; no carga por defecto el backend nativo en su propio proceso.',
      'La tarea se serializa como una solicitud acotada; no se envía estado arbitrario de la interfaz.',
      'El worker importa RDKit dentro de otro proceso, ejecuta el trabajo y queda sujeto a timeout.',
      'La UI interpreta una respuesta, un error o una interrupción y puede degradar a una alternativa.'
    ];
    buttons.forEach(button => button.addEventListener('click', () => {
      const step = Number(button.dataset.workerStep);
      diagram.dataset.step = String(step);
      document.getElementById('worker-explainer').textContent = text[step];
      buttons.forEach(other => {
        const active = other === button;
        other.classList.toggle('is-active', active);
        other.setAttribute('aria-pressed', String(active));
      });
    }));
  }

  function setupScreenshots() {
    const image = document.getElementById('ui-screenshot');
    const shots = {
      reference: ['assets/ui/foundation-light.png', 'Captura real de ChemUSON en la etapa Foundation de la modernización visual; ventana clara con la interfaz previa al panel lateral moderno.'],
      current: ['assets/ui/production-light.png', 'Captura real de producción de ChemUSON tras la convergencia visual, con el tema claro y el panel lateral moderno.'],
      dark: ['assets/ui/production-dark.png', 'Captura real de producción de ChemUSON tras la convergencia visual, con el tema oscuro y el panel lateral moderno.']
    };
    document.querySelectorAll('[data-shot]').forEach(button => button.addEventListener('click', () => {
      const [src, alt] = shots[button.dataset.shot];
      image.src = src;
      image.alt = alt;
      document.querySelectorAll('[data-shot]').forEach(other => {
        const active = other === button;
        other.classList.toggle('is-active', active);
        other.setAttribute('aria-pressed', String(active));
      });
    }));
  }

  function setupTheme() {
    const key = 'chemuson-article-theme';
    const saved = localStorage.getItem(key);
    const initial = saved || (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    document.documentElement.dataset.theme = initial;
    const button = document.getElementById('theme-toggle');
    const icon = button.querySelector('span');
    function setTheme(theme) {
      document.documentElement.dataset.theme = theme;
      localStorage.setItem(key, theme);
      button.setAttribute('aria-label', theme === 'dark' ? 'Cambiar a tema claro' : 'Cambiar a tema oscuro');
      button.title = theme === 'dark' ? 'Cambiar a tema claro' : 'Cambiar a tema oscuro';
      icon.textContent = theme === 'dark' ? '☼' : '◐';
    }
    setTheme(initial);
    button.addEventListener('click', () => setTheme(document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark'));
  }

  function setupReadProgress() {
    const tocLinks = [...document.querySelectorAll('.toc a')];
    const sections = tocLinks.map(link => document.querySelector(link.getAttribute('href'))).filter(Boolean);
    if (!('IntersectionObserver' in window)) return;
    const observer = new IntersectionObserver(entries => {
      const visible = entries.filter(entry => entry.isIntersecting).sort((a, b) => b.intersectionRatio - a.intersectionRatio)[0];
      if (!visible) return;
      tocLinks.forEach(link => {
        const current = link.getAttribute('href') === `#${visible.target.id}`;
        link.setAttribute('aria-current', current ? 'location' : 'false');
        link.classList.toggle('is-current', current);
      });
    }, { rootMargin: '-18% 0px -68% 0px', threshold: [0, .2, .5, 1] });
    sections.forEach(section => observer.observe(section));
  }

  setupTheme();
  setupHero();
  setupStages();
  setupGeometrySlider();
  setupLayers();
  setupPersistenceFormats();
  setupArchitecture();
  setupWorker();
  setupScreenshots();
  setupReadProgress();
})();
