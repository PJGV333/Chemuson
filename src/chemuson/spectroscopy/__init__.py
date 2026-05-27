"""Predictores espectrales desacoplados."""

from .service import (
    CarbonNmrPeak,
    MassPeak,
    ProtonNmrPeak,
    SpectralPrediction,
    SpectrumPredictor,
    predict_spectra,
    register_predictor,
)

__all__ = [
    "CarbonNmrPeak",
    "MassPeak",
    "ProtonNmrPeak",
    "SpectralPrediction",
    "SpectrumPredictor",
    "predict_spectra",
    "register_predictor",
]
