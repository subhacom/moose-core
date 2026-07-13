# base.py ---
#
# Filename: base.py
# Description:
# Author: Subhasis Ray
# Created: Mon Jan 19 00:03:40 2026 (+0530)
#

# Code:
from typing import List, Optional, Protocol

import moose


class FormatHandler(Protocol):
    """Protocol for model format handlers."""

    extensions: tuple[str, ...]

    def read(self, filepath: str, loadpath: str, **options) -> moose.ObjId: ...

    def write(self, modelpath: str, filepath: str, **options) -> moose.ObjId: ...


class ModelLoadError(Exception):
    """Failed to load model."""

    def __init__(self, message, filepath=None, loadpath=None):
        super().__init__(message)
        self.filepath = filepath
        self.loadpath = loadpath


#
# base.py ends here
