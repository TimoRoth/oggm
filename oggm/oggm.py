from .utils import Utils

from ._cfg import CFG
from ._workflow import Workflow


class OGGM(CFG, Workflow):
    def __init__(self, **kwargs):
        CFG.__init__(self, **kwargs)
        Workflow.__init__(self, **kwargs)

        self.utils = self.u = Utils(self)

