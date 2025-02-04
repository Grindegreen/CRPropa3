import os
import sys

# Import the compiled CRPropa module (generated from SWIG or pybind11)
try:
    from .crpropa_py import *  # Imports all functions, classes, etc.
except ImportError as e:
    raise ImportError(
        f"Could not import CRPropa. Ensure it is installed correctly. Original error: {e}"
    )

# Define the package metadata
__version__ = "3.2.0"  # Update this dynamically if needed
__author__ = "CRPropa Collaboration"
__license__ = "GPL-3.0"

# Optional: Expose key functionalities in a structured way
__all__ = [
    "Module1",
    "Module2",
    "PropagationEngine",
    # List other main modules/classes/functions here
]

# Print a message to confirm successful import (useful for debugging)
if os.getenv("CRPROPA_DEBUG"):
    print("CRPropa successfully imported!")