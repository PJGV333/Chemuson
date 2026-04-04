"""Comandos de edición con soporte de deshacer/rehacer (QUndoCommand)."""
from __future__ import annotations

from .atom_commands import (
    AddAtomCommand,
    ChangeAtomCommand,
    ChangeChargeCommand,
    ChangeNoImplicitCommand,
    ChangeAtomLabelScaleCommand,
    SetCoordinationCenterCommand,
    ChangeCoordinationSphereStyleCommand,
)

from .bond_commands import (
    AddBondCommand,
    ChangeBondCommand,
    ChangeBondLengthCommand,
    ChangeBondStrokeCommand,
    ChangeBondColorCommand,
    ChangeDoubleBondOrientationCommand,
    AddRingCommand,
    AddChainCommand,
)

from .text_commands import (
    MoveTextItemsCommand,
    ScaleTextItemsCommand,
    FormatTextItemsCommand,
    AddTextItemCommand,
)

from .annotation_commands import (
    ChangeArrowStrokeCommand,
    ChangeBracketStrokeCommand,
    AddArrowCommand,
    AddBracketCommand,
    AddImageItemCommand,
    AddWavyAnchorCommand,
)

from .diagram_commands import (
    AddEnergyDiagramItemCommand,
    TransformEnergyDiagramItemsCommand,
    ConfigureEnergyDiagramItemsCommand,
    AddCompositeDiagramItemCommand,
    EditSemanticDiagramCommand,
    AddOrbitalItemCommand,
    TransformOrbitalItemsCommand,
    StyleOrbitalItemsCommand,
)

from .plate_commands import (
    AddPlateItemCommand,
    MovePlateItemsCommand,
    TransformPlateItemsCommand,
    AddSpotBandCommand,
    RemoveSpotBandCommand,
    MoveSpotBandCommand,
    ChangeSpotBandPropertyCommand,
    ChangePlateLabelsCommand,
)

from .transform_commands import (
    ChangeCanvasOpacityCommand,
    MoveAtomsCommand,
    MoveArrowItemsCommand,
    MoveBracketItemsCommand,
    TransformImageItemsCommand,
    ScaleArrowItemsCommand,
    ScaleBracketItemsCommand,
    DeleteSelectionCommand,
)

__all__ = [
    "AddAtomCommand",
    "ChangeAtomCommand",
    "ChangeChargeCommand",
    "ChangeNoImplicitCommand",
    "ChangeAtomLabelScaleCommand",
    "SetCoordinationCenterCommand",
    "ChangeCoordinationSphereStyleCommand",
    "AddBondCommand",
    "ChangeBondCommand",
    "ChangeBondLengthCommand",
    "ChangeBondStrokeCommand",
    "ChangeBondColorCommand",
    "ChangeDoubleBondOrientationCommand",
    "AddRingCommand",
    "AddChainCommand",
    "MoveTextItemsCommand",
    "ScaleTextItemsCommand",
    "FormatTextItemsCommand",
    "AddTextItemCommand",
    "ChangeArrowStrokeCommand",
    "ChangeBracketStrokeCommand",
    "AddArrowCommand",
    "AddBracketCommand",
    "AddImageItemCommand",
    "AddWavyAnchorCommand",
    "AddEnergyDiagramItemCommand",
    "TransformEnergyDiagramItemsCommand",
    "ConfigureEnergyDiagramItemsCommand",
    "AddCompositeDiagramItemCommand",
    "EditSemanticDiagramCommand",
    "AddOrbitalItemCommand",
    "TransformOrbitalItemsCommand",
    "StyleOrbitalItemsCommand",
    "AddPlateItemCommand",
    "MovePlateItemsCommand",
    "TransformPlateItemsCommand",
    "AddSpotBandCommand",
    "RemoveSpotBandCommand",
    "MoveSpotBandCommand",
    "ChangeSpotBandPropertyCommand",
    "ChangePlateLabelsCommand",
    "ChangeCanvasOpacityCommand",
    "MoveAtomsCommand",
    "MoveArrowItemsCommand",
    "MoveBracketItemsCommand",
    "TransformImageItemsCommand",
    "ScaleArrowItemsCommand",
    "ScaleBracketItemsCommand",
    "DeleteSelectionCommand",
]
