"""  Configuration file and options

A number of globals are defined here to be available everywhere.
"""
import os
import numpy as np
from scipy.signal import gaussian

# Machine epsilon
FLOAT_EPS = np.finfo(float).eps

# Constants
SEC_IN_YEAR = 365*24*3600
SEC_IN_DAY = 24*3600
SEC_IN_HOUR = 3600
SEC_IN_MONTH = 2628000
DAYS_IN_MONTH = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

G = 9.80665  # gravity

GAUSSIAN_KERNEL = dict()
for ks in [5, 7, 9]:
    kernel = gaussian(ks, 1)
    GAUSSIAN_KERNEL[ks] = kernel / kernel.sum()

# Path to the config file
DEFAULT_CONFIG_FILE = os.path.join(os.path.expanduser('~'), '.oggm_config')

# Path to the cache directory
DEFAULT_CACHE_DIR = os.path.join(os.path.expanduser('~'), '.oggm')


class DocumentedDict(dict):
    """Quick "magic" to document the BASENAMES entries."""

    def __init__(self):
        self._doc = dict()

    def _set_key(self, key, value, docstr=''):
        if key in self:
            raise ValueError('Cannot overwrite a key.')
        dict.__setitem__(self, key, value)
        self._doc[key] = docstr

    def __setitem__(self, key, value):
        # Overrides the original dic to separate value and documentation
        try:
            self._set_key(key, value[0], docstr=value[1])
        except BaseException:
            raise ValueError('DocumentedDict accepts only tuple of len 2')

    def info_str(self, key):
        """Info string for the documentation."""
        return '    {}'.format(self[key]) + '\n' + '        ' + self._doc[key]

    def doc_str(self, key):
        """Info string for the documentation."""
        return '        {}'.format(self[key]) + '\n' + '            ' + \
               self._doc[key]


BASENAMES = DocumentedDict()

_doc = ('A geotiff file containing the DEM (reprojected into the local grid).'
        'This DEM is not smoothed or gap filles, and is the closest to the '
        'original DEM source.')
BASENAMES['dem'] = ('dem.tif', _doc)

_doc = ('A glacier mask geotiff file with the same extend and projection as '
        'the `dem.tif`. This geotiff has value 1 at glaciated grid points and '
        ' value 0 at unglaciated points.')
BASENAMES['glacier_mask'] = ('glacier_mask.tif', _doc)

_doc = ('The glacier outlines in the local map projection (Transverse '
        'Mercator).')
BASENAMES['outlines'] = ('outlines.shp', _doc)

_doc = ('The glacier intersects in the local map projection (Transverse '
        'Mercator).')
BASENAMES['intersects'] = ('intersects.shp', _doc)

_doc = ('Each flowline has a catchment area computed from flow routing '
        'algorithms: this shapefile stores the catchment outlines (in the '
        'local map projection (Transverse Mercator).')
BASENAMES['flowline_catchments'] = ('flowline_catchments.shp', _doc)

_doc = ('The intersections between cathments (shapefile) in the local map '
        'projection (Transverse Mercator).')
BASENAMES['catchments_intersects'] = ('catchments_intersects.shp', _doc)

_doc = 'A ``salem.Grid`` handling the georeferencing of the local grid.'
BASENAMES['glacier_grid'] = ('glacier_grid.json', _doc)

_doc = 'A dictionary containing runtime diagnostics useful for debugging.'
BASENAMES['diagnostics'] = ('diagnostics.json', _doc)

_doc = ('A netcdf file containing several gridded data variables such as '
        'topography, the glacier masks, the interpolated 2D glacier bed, '
        'and more.')
BASENAMES['gridded_data'] = ('gridded_data.nc', _doc)

_doc = ('A dictionary containing the shapely.Polygons of a glacier. The '
        '"polygon_hr" entry contains the geometry transformed to the local '
        'grid in (i, j) coordinates, while the "polygon_pix" entry contains '
        'the geometries transformed into the coarse grid (the i, j elements '
        'are integers). The "polygon_area" entry contains the area of the '
        'polygon as computed by Shapely. The "catchment_indices" entry'
        'contains a list of len `n_centerlines`, each element containing '
        'a numpy array of the indices in the glacier grid which represent '
        'the centerlines catchment area.')
BASENAMES['geometries'] = ('geometries.pkl', _doc)

_doc = ('A dictionary containing the downsteam line geometry as well as the '
        'bed shape computed from a parabolic fit.')
BASENAMES['downstream_line'] = ('downstream_line.pkl', _doc)

_doc = 'A text file with the source of the topo file (GIMP, SRTM, ...).'
BASENAMES['dem_source'] = ('dem_source.txt', _doc)

_doc = ('A hypsometry file computed by OGGM and provided in the same format '
        'as the RGI (useful for diagnostics).')
BASENAMES['hypsometry'] = ('hypsometry.csv', _doc)

_doc = 'A list of :py:class:`oggm.Centerline` instances, sorted by flow order.'
BASENAMES['centerlines'] = ('centerlines.pkl', _doc)

_doc = ('A "better" version of the Centerlines, now on a regular spacing '
        'i.e., not on the gridded (i, j) indices. The tails of the '
        'tributaries are cut out to make more realistic junctions. '
        'They are now "1.5D" i.e., with a width.')
BASENAMES['inversion_flowlines'] = ('inversion_flowlines.pkl', _doc)

_doc = 'The historical monthly climate timeseries stored in a netCDF file.'
BASENAMES['climate_historical'] = ('climate_historical.nc', _doc)

_doc = 'Deprecated: old name for `climate_historical`.'
BASENAMES['climate_monthly'] = ('climate_monthly.nc', _doc)

_doc = ('Some information (dictionary) about the mass '
        'balance parameters for this glacier.')
BASENAMES['climate_info'] = ('climate_info.json', _doc)

_doc = 'The monthly GCM climate timeseries stored in a netCDF file.'
BASENAMES['gcm_data'] = ('gcm_data.nc', _doc)

_doc = "A dict containing the glacier's t*, bias, and the flowlines' mu*"
BASENAMES['local_mustar'] = ('local_mustar.json', _doc)

_doc = 'List of dicts containing the data needed for the inversion.'
BASENAMES['inversion_input'] = ('inversion_input.pkl', _doc)

_doc = 'List of dicts containing the output data from the inversion.'
BASENAMES['inversion_output'] = ('inversion_output.pkl', _doc)

_doc = 'Dict of fs and fd as computed by the inversion optimisation.'
BASENAMES['inversion_params'] = ('inversion_params.pkl', _doc)

_doc = 'List of flowlines ready to be run by the model.'
BASENAMES['model_flowlines'] = ('model_flowlines.pkl', _doc)

_doc = ('When using a linear mass-balance for the inversion, this dict stores '
        'the optimal ela_h and grad.')
BASENAMES['linear_mb_params'] = ('linear_mb_params.pkl', _doc)

_doc = ('A netcdf file containing enough information to reconstruct the '
        'entire flowline glacier along the run (can be data expensive).')
BASENAMES['model_run'] = ('model_run.nc', _doc)

_doc = ('A netcdf file containing the model diagnostics (volume, '
        'mass-balance, length...).')
BASENAMES['model_diagnostics'] = ('model_diagnostics.nc', _doc)

_doc = "A dict containing the glacier's t*, bias, mu*. Analogous " \
       "to 'local_mustar.json', but for the volume/area scaling model."
BASENAMES['vascaling_mustar'] = ('vascaling_mustar.json', _doc)


def add_to_basenames(basename, filename, docstr=''):
    """Add an entry to the list of BASENAMES.

    BASENAMES are access keys to files available at the gdir level.

    Parameters
    ----------
    basename : str
        the key (e.g. 'dem', 'model_flowlines')
    filename : str
        the associated filename (e.g. 'dem.tif', 'model_flowlines.pkl')
    docstr : str
        the associated docstring (for documentation)
    """
    global BASENAMES
    if '.' not in filename:
        raise ValueError('The filename needs a proper file suffix!')
    BASENAMES[basename] = (filename, docstr)
