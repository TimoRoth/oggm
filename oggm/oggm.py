from . import utils

from ._cfg import CFG
from ._workflow import Workflow


class OGGM(CFG, Workflow):
    def __init__(self, **kwargs):
        CFG.__init__(self, **kwargs)
        Workflow.__init__(self, **kwargs)

        self.utils = self.u = utils.Utils(self)


    # Convenience-Wrapper
    # Allows accessing the entirety of utils via the OGGM object
    # without needing to double-check where exactly each function is
    # each and every time.
    def __getattr__(self, name):
        try:
            return getattr(self.utils, name)
        except AttributeError:
            return getattr(utils, name)
