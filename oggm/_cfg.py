import logging
import os
import shutil
import sys
import glob
import json
from collections import OrderedDict
from distutils.util import strtobool

import numpy as np
import pandas as pd
from configobj import ConfigObj, ConfigObjError
try:
    import geopandas as gpd
except ImportError:
    pass
try:
    import salem
except ImportError:
    pass

from oggm.exceptions import InvalidParamsError
from oggm.cfg import DEFAULT_CONFIG_FILE, DEFAULT_CACHE_DIR

# Local logger
log = logging.getLogger(__name__)


class ResettingOrderedDict(OrderedDict):
    """OrderedDict wrapper that resets our multiprocessing on set"""

    def __init__(self, oggm, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.oggm = oggm

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        self.oggm.CONFIG_MODIFIED = True


class PathOrderedDict(ResettingOrderedDict):
    """Quick "magic" to be sure that paths are expanded correctly."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        # Overrides the original dic to expand the path
        try:
            value = os.path.expanduser(value)
        except AttributeError:
            raise InvalidParamsError('The value you are trying to set does '
                                     'not seem to be a valid path: '
                                     '{}'.format(value))

        ResettingOrderedDict.__setitem__(self, key, value)


class ParamsLoggingDict(ResettingOrderedDict):
    """Quick "magic" to log the parameter changes by the user."""

    do_log = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        # Overrides the original dic to log the change
        if self.do_log:
            self._log_param_change(key, value)
        ResettingOrderedDict.__setitem__(self, key, value)

    def _log_param_change(self, key, value):

        prev = self.get(key)
        if prev is None:
            if key in ['baseline_y0', 'baseline_y1']:
                raise InvalidParamsError('The `baseline_y0` and `baseline_y1` '
                                         'parameters have been removed. '
                                         'You now have to set them explicitly '
                                         'in your call to '
                                         '`process_climate_data`.')

            log.warning('WARNING: adding an unknown parameter '
                        '`{}`:`{}` to PARAMS.'.format(key, value))
            return

        if prev == value:
            return

        log.workflow("PARAMS['{}'] changed from `{}` to `{}`.".format(key,
                                                                      prev,
                                                                      value))


class CFG:
    def __init__(self, config_file=None, cache_dir=None, **kwargs):
        self.CONFIG_FILE = config_file or DEFAULT_CONFIG_FILE
        self.CACHE_DIR = cache_dir or DEFAULT_CACHE_DIR

        self.CONFIG_MODIFIED = False

        self.IS_INITIALIZED = False
        self.PARAMS = ParamsLoggingDict(self)
        self.PATHS = PathOrderedDict(self)
        self.LRUHANDLERS = ResettingOrderedDict(self)
        self.DATA = ResettingOrderedDict(self)

        # TODO
        self.DL_VERIFIED = dict()


    def oggm_static_paths(self):
        """Initialise the OGGM paths from the config file."""

        # Create cache dir if it doesn't exist
        if not os.path.exists(self.CACHE_DIR):
            os.makedirs(self.CACHE_DIR)

        # See if the file is there, if not create it
        if not os.path.exists(self.CONFIG_FILE):
            dldir = os.path.join(os.path.expanduser('~'), 'OGGM')
            config = ConfigObj()
            config['dl_cache_dir'] = os.path.join(dldir, 'download_cache')
            config['dl_cache_readonly'] = False
            config['tmp_dir'] = os.path.join(dldir, 'tmp')
            config['rgi_dir'] = os.path.join(dldir, 'rgi')
            config['test_dir'] = os.path.join(dldir, 'tests')
            config['has_internet'] = True
            config.filename = self.CONFIG_FILE
            config.write()

        # OK, read in the file
        try:
            config = ConfigObj(self.CONFIG_FILE, file_error=True)
        except (ConfigObjError, IOError) as e:
            log.critical('Config file could not be parsed (%s): %s',
                         self.CONFIG_FILE, e)
            raise RuntimeError('Config file could not be parsed')

        # Check that all keys are here
        for k in ['dl_cache_dir', 'dl_cache_readonly', 'tmp_dir',
                  'rgi_dir', 'test_dir', 'has_internet']:
            if k not in config:
                raise InvalidParamsError('The oggm config file ({}) should have '
                                         'an entry for {}.'.format(self.CONFIG_FILE, k))

        # Override defaults with env variables if available
        if os.environ.get('OGGM_DOWNLOAD_CACHE_RO') is not None:
            ro = bool(strtobool(os.environ.get('OGGM_DOWNLOAD_CACHE_RO')))
            config['dl_cache_readonly'] = ro
        if os.environ.get('OGGM_DOWNLOAD_CACHE') is not None:
            config['dl_cache_dir'] = os.environ.get('OGGM_DOWNLOAD_CACHE')
        if os.environ.get('OGGM_EXTRACT_DIR') is not None:
            # This is for the directories where OGGM needs to extract things
            # On the cluster it might be useful to do it on a fast disc
            edir = os.path.abspath(os.environ.get('OGGM_EXTRACT_DIR'))
            config['tmp_dir'] = os.path.join(edir, 'tmp')
            config['rgi_dir'] = os.path.join(edir, 'rgi')

        # Fill the PATH dict
        for k, v in config.iteritems():
            if not k.endswith('_dir'):
                continue
            self.PATHS[k] = os.path.abspath(os.path.expanduser(v))

        # Other
        self.PARAMS.do_log = False
        self.PARAMS['has_internet'] = config.as_bool('has_internet')
        self.PARAMS['dl_cache_readonly'] = config.as_bool('dl_cache_readonly')
        self.PARAMS.do_log = True

        # Create cache dir if possible
        if not os.path.exists(self.PATHS['dl_cache_dir']):
            if not self.PARAMS['dl_cache_readonly']:
                os.makedirs(self.PATHS['dl_cache_dir'])


    def set_logging_config(self, logging_level='INFO'):
        """Set the global logger parameters.

        Logging levels:

        DEBUG
            Print detailed information, typically of interest only when diagnosing
            problems.
        INFO
            Print confirmation that things are working as expected, e.g. when
            each task is run correctly (this is the default).
        WARNING
            Indication that something unexpected happened on a glacier,
            but that OGGM is still working on this glacier.
        WORKFLOW
            Print only high level, workflow information (typically, one message
            per task). Errors and warnings will still be printed.
        ERROR
            Print errors only, e.g. when a glacier cannot run properly.
        CRITICAL
            Print nothing but fatal errors.

        Parameters
        ----------
        logging_level : str or None
            the logging level. See description above for a list of options. Setting
            to `None` is equivalent to `'CRITICAL'`, i.e. no log output will be
            generated.
        """

        # Add a custom level - just for us
        logging.addLevelName(25, 'WORKFLOW')

        def workflow(self, message, *args, **kws):
            """Standard log message with a custom level."""
            if self.isEnabledFor(25):
                # Yes, logger takes its '*args' as 'args'.
                self._log(25, message, args, **kws)

        logging.WORKFLOW = 25
        logging.Logger.workflow = workflow

        # Remove all handlers associated with the root logger object.
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        # Spammers
        logging.getLogger("Fiona").setLevel(logging.CRITICAL)
        logging.getLogger("fiona").setLevel(logging.CRITICAL)
        logging.getLogger("shapely").setLevel(logging.CRITICAL)
        logging.getLogger("rasterio").setLevel(logging.CRITICAL)
        logging.getLogger("matplotlib").setLevel(logging.CRITICAL)
        logging.getLogger("numexpr").setLevel(logging.CRITICAL)

        # Basic config
        if logging_level is None:
            logging_level = 'CRITICAL'
        logging_level = logging_level.upper()
        logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=getattr(logging, logging_level))


    def initialize_minimal(self, file=None, logging_level='INFO'):
        """Same as initialise() but without requiring any download of data.

        This is useful for "flowline only" OGGM applications

        Parameters
        ----------
        file : str
            path to the configuration file (default: OGGM params.cfg)
        logging_level : str
            set a logging level. See :func:`set_logging_config` for options.
        """

        self.set_logging_config(logging_level=logging_level)

        if file is None:
            file = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                'params.cfg')

        try:
            cp = ConfigObj(file, file_error=True)
        except (ConfigObjError, IOError) as e:
            log.critical('Config file could not be parsed (%s): %s', file, e)
            raise RuntimeError('Config file could not be parsed')

        log.workflow('Using configuration file: %s', file)

        self.oggm_static_paths()

        # Paths
        self.PATHS['working_dir'] = cp['working_dir']
        self.PATHS['dem_file'] = cp['dem_file']
        self.PATHS['climate_file'] = cp['climate_file']

        # Do not spam
        self.PARAMS.do_log = False

        # Size of LRU cache
        try:
            lru_maxsize = int(os.environ['LRU_MAXSIZE'])
            log.workflow('Size of LRU cache set to {} '.format(lru_maxsize) +
                         'according to the ENV variable LRU_MAXSIZE')
        except KeyError:
            lru_maxsize = cp.as_int('lru_maxsize')
        self.PARAMS['lru_maxsize'] = lru_maxsize

        # Some non-trivial params
        self.PARAMS['continue_on_error'] = cp.as_bool('continue_on_error')
        self.PARAMS['grid_dx_method'] = cp['grid_dx_method']
        self.PARAMS['topo_interp'] = cp['topo_interp']
        self.PARAMS['use_intersects'] = cp.as_bool('use_intersects')
        self.PARAMS['use_compression'] = cp.as_bool('use_compression')
        self.PARAMS['border'] = cp.as_int('border')
        self.PARAMS['use_multiple_flowlines'] = cp.as_bool('use_multiple_flowlines')
        self.PARAMS['filter_min_slope'] = cp.as_bool('filter_min_slope')
        self.PARAMS['auto_skip_task'] = cp.as_bool('auto_skip_task')
        self.PARAMS['correct_for_neg_flux'] = cp.as_bool('correct_for_neg_flux')
        self.PARAMS['filter_for_neg_flux'] = cp.as_bool('filter_for_neg_flux')
        self.PARAMS['run_mb_calibration'] = cp.as_bool('run_mb_calibration')
        self.PARAMS['rgi_version'] = cp['rgi_version']
        self.PARAMS['use_rgi_area'] = cp.as_bool('use_rgi_area')
        self.PARAMS['compress_climate_netcdf'] = cp.as_bool('compress_climate_netcdf')
        self.PARAMS['use_tar_shapefiles'] = cp.as_bool('use_tar_shapefiles')
        self.PARAMS['clip_mu_star'] = cp.as_bool('clip_mu_star')
        self.PARAMS['clip_tidewater_border'] = cp.as_bool('clip_tidewater_border')
        self.PARAMS['dl_verify'] = cp.as_bool('dl_verify')
        self.PARAMS['calving_line_extension'] = cp.as_int('calving_line_extension')
        k = 'use_kcalving_for_inversion'
        self.PARAMS[k] = cp.as_bool(k)
        self.PARAMS['use_kcalving_for_run'] = cp.as_bool('use_kcalving_for_run')
        self.PARAMS['calving_use_limiter'] = cp.as_bool('calving_use_limiter')
        k = 'error_when_glacier_reaches_boundaries'
        self.PARAMS[k] = cp.as_bool(k)

        # Climate
        self.PARAMS['baseline_climate'] = cp['baseline_climate'].strip().upper()
        self.PARAMS['hydro_month_nh'] = cp.as_int('hydro_month_nh')
        self.PARAMS['hydro_month_sh'] = cp.as_int('hydro_month_sh')
        self.PARAMS['climate_qc_months'] = cp.as_int('climate_qc_months')
        self.PARAMS['temp_use_local_gradient'] = cp.as_bool('temp_use_local_gradient')
        self.PARAMS['tstar_search_glacierwide'] = cp.as_bool('tstar_search_glacierwide')

        k = 'temp_local_gradient_bounds'
        self.PARAMS[k] = [float(vk) for vk in cp.as_list(k)]
        k = 'tstar_search_window'
        self.PARAMS[k] = [int(vk) for vk in cp.as_list(k)]
        self.PARAMS['use_bias_for_run'] = cp.as_bool('use_bias_for_run')
        k = 'free_board_marine_terminating'
        self.PARAMS[k] = [float(vk) for vk in cp.as_list(k)]

        # Inversion
        k = 'use_shape_factor_for_inversion'
        self.PARAMS[k] = cp[k]

        # Flowline model
        k = 'use_shape_factor_for_fluxbasedmodel'
        self.PARAMS[k] = cp[k]

        # Delete non-floats
        ltr = ['working_dir', 'dem_file', 'climate_file', 'use_tar_shapefiles',
               'grid_dx_method', 'run_mb_calibration', 'compress_climate_netcdf',
               'climate_qc_months', 'temp_use_local_gradient', 'temp_local_gradient_bounds',
               'topo_interp', 'use_compression', 'bed_shape', 'continue_on_error',
               'use_multiple_flowlines', 'tstar_search_glacierwide',
               'border', 'hydro_month_nh', 'clip_mu_star',
               'tstar_search_window', 'use_bias_for_run', 'hydro_month_sh',
               'use_intersects', 'filter_min_slope', 'clip_tidewater_border',
               'auto_skip_task', 'correct_for_neg_flux', 'filter_for_neg_flux',
               'rgi_version', 'dl_verify', 'use_mp_spawn', 'calving_use_limiter',
               'use_shape_factor_for_inversion', 'use_rgi_area',
               'use_shape_factor_for_fluxbasedmodel', 'baseline_climate',
               'calving_line_extension', 'use_kcalving_for_run', 'lru_maxsize',
               'free_board_marine_terminating', 'use_kcalving_for_inversion',
               'error_when_glacier_reaches_boundaries']
        for k in ltr:
            cp.pop(k, None)

        # Other params are floats
        for k in cp:
            self.PARAMS[k] = cp.as_float(k)
        self.PARAMS.do_log = True

        # Empty defaults
        self.set_intersects_db()
        self.IS_INITIALIZED = True


    def initialize(self, file=None, logging_level='INFO'):
        """Read the configuration file containing the run's parameters.

        This should be the first call, before using any of the other OGGM modules
        for most (all?) OGGM simulations.

        Parameters
        ----------
        file : str
            path to the configuration file (default: OGGM params.cfg)
        logging_level : str
            set a logging level. See :func:`set_logging_config` for options.
        """
        self.initialize_minimal(file=file, logging_level=logging_level)

        # Do not spam
        self.PARAMS.do_log = False

        # Make sure we have a proper cache dir
        self.utils.download_oggm_files()

        # Read-in the reference t* data for all available models types (oggm, vas)
        model_prefixes = ['oggm_', 'vas_']
        for prefix in model_prefixes:
            fns = ['ref_tstars_rgi5_cru4', 'ref_tstars_rgi6_cru4',
                   'ref_tstars_rgi5_histalp', 'ref_tstars_rgi6_histalp']
            for fn in fns:
                fpath = self.utils.get_demo_file(prefix + fn + '.csv')
                self.PARAMS[prefix + fn] = pd.read_csv(fpath)
                fpath = self.utils.get_demo_file(prefix + fn + '_calib_params.json')
                with open(fpath, 'r') as fp:
                    mbpar = json.load(fp)
                self.PARAMS[prefix + fn + '_calib_params'] = mbpar

        # Read in the demo glaciers
        file = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'data', 'demo_glaciers.csv')
        self.DATA['demo_glaciers'] = pd.read_csv(file, index_col=0)

        # Add other things
        if 'dem_grids' not in self.DATA:
            grids = {}
            for grid_json in ['gimpdem_90m_v01.1.json',
                              'arcticdem_mosaic_100m_v3.0.json',
                              'Alaska_albers_V3.json',
                              'AntarcticDEM_wgs84.json',
                              'REMA_100m_dem.json']:
                if grid_json not in grids:
                    fp = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                      'data', grid_json)
                    try:
                        grids[grid_json] = salem.Grid.from_json(fp)
                    except NameError:
                        pass
            self.DATA['dem_grids'] = grids

        # Trigger a one time check of the hash file
        self.utils.get_dl_verify_data('dummy_section')

        # OK
        self.PARAMS.do_log = True


    def get_lru_handler(self, tmpdir=None, maxsize=None, ending='.tif'):
        """LRU handler for a given temporary directory (singleton).

        Parameters
        ----------
        tmpdir : str
            path to the temporary directory to handle. Default is
            ``oggm.PATHS['tmp_dir']``.
        maxsize : int
            the max number of files to keep in the directory
        ending : str
            consider only the files with a certain ending
        """

        # see if we're set up
        if tmpdir is None:
            tmpdir = self.PATHS['tmp_dir']
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)

        # one handler per directory and file ending
        # (in practice not very useful, but a dict is easier to handle)
        k = (tmpdir, ending)
        if k in self.LRUHANDLERS:
            # was already there
            lru = self.LRUHANDLERS[k]
            # possibility to increase or decrease the cachesize if need be
            if maxsize is not None:
                lru.maxsize = maxsize
                lru.purge()
            return lru
        else:
            # we do a new one
            from oggm.utils import LRUFileCache
            # the files already present have to be counted, too
            l0 = list(glob.glob(os.path.join(tmpdir, '*' + ending)))
            l0.sort(key=os.path.getctime)
            lru = LRUFileCache(self, l0, maxsize=maxsize)
            self.LRUHANDLERS[k] = lru
            return lru


    def set_intersects_db(self, path_or_gdf=None):
        """Set the glacier intersection database for OGGM to use.

        It is now set automatically by the
        :func:`oggm.workflow.init_glacier_directories` task, but setting it
        manually can be useful for a slightly faster run initialization.

        See :func:`oggm.utils.get_rgi_intersects_region_file` for how to obtain
        such data.

        Parameters
        ----------
        path_or_gdf : str of geopandas.GeoDataframe
            the intersects file to use
        """

        self.PARAMS.do_log = False

        if self.PARAMS['use_intersects'] and path_or_gdf is not None:
            if isinstance(path_or_gdf, str):
                self.PARAMS['intersects_gdf'] = gpd.read_file(path_or_gdf)
            else:
                self.PARAMS['intersects_gdf'] = path_or_gdf
        else:
            self.PARAMS['intersects_gdf'] = pd.DataFrame()
        self.PARAMS.do_log = True


    def reset_working_dir(self):
        """Deletes the content of the working directory. Careful: cannot be undone!
        """
        if self.PATHS['working_dir']:
            if os.path.exists(self.PATHS['working_dir']):
                shutil.rmtree(self.PATHS['working_dir'])
            os.makedirs(self.PATHS['working_dir'])
