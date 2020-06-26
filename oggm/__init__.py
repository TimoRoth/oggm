""" OGGM package.

Copyright: OGGM e.V. and OGGM Contributors

License: BSD-3-Clause
"""
# flake8: noqa
from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
finally:
    del get_distribution, DistributionNotFound

from .oggm import OGGM
