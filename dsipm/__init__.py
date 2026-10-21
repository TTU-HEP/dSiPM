"""Small DREAM/FLASH photon-response and depth-reconstruction core.

Use geometry, sensors, depth and io directly. No plotting, ROOT import or
production submission occurs at package import time.
"""
from .depth import DepthContract, PhotonObservations, fit_calibration, reconstruct

__all__ = ["DepthContract", "PhotonObservations", "fit_calibration", "reconstruct"]
