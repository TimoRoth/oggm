# flake8: noqa

from ._funcs import *
from ._downloads import *
from ._workflow import *

from ._funcs import Funcs
from ._downloads import Downloads
from ._workflow import Workflow

class Utils(Funcs, Downloads, Workflow):
    def __init__(self, oggm, **kwargs):
        self.oggm = oggm
        
        Funcs.__init__(self, oggm, **kwargs)
        Downloads.__init__(self, oggm, **kwargs)
        Workflow.__init__(self, oggm, **kwargs)
