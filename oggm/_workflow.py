# Built ins
import logging
import os
import inspect
from shutil import rmtree
from collections.abc import Sequence
from functools import partial, wraps
# External libs
import multiprocessing
import numpy as np

# Locals
from oggm import utils
from oggm.exceptions import InvalidParamsError

# Module logger
log = logging.getLogger(__name__)


class Workflow:
    def __init__(self, **kwargs):
        pass


    def execute_entity_task(self, task, gdirs):
        """Execute a task on gdirs.

        If ``task`` has more arguments than `gdir` they have to be keyword
        arguments.

        Returns the processed (delayed) gdirs.

        Parameters
        ----------
        task : various
             the entity task to apply
             Can either be an instance of EntityTask or GlobalTask.
             If neither of the two, it will be assumed that calling
             task(oggm) will result in such an instance.

             Can optionally also be a list of the above,
             which will chain the execution of the tasks.
        gdirs : list of :py:class:`oggm.GlacierDirectory` objects
            the glacier directories to process
        """

        if not task:
            raise InvalidParamsError('task is invalid')

        # Should be iterable
        gdirs = utils.tolist(gdirs)
        tasks = utils.tolist(task)

        if len(gdirs) == 0:
            return

        for task in tasks:
            if not isinstance(task, (utils.GlobalTask, utils.EntityTask)):
                task = task(self)

            if not isinstance(task, (utils.GlobalTask, utils.EntityTask)):
                raise InvalidParamsError('Invalid Entity or Global Task provided')

            log.workflow('Execute entity task %s on %d glaciers', task.task_name, len(gdirs))

            if isinstance(task, utils.GlobalTask):
                gdirs = task.run(gdirs)
            else:
                gdirs = utils.EntityTaskResultList(None, [task.run(gdir) for gdir in gdirs])

        return gdirs


    def execute_parallel_tasks(self, gdir, tasks):
        """Execute a list of task on a single gdir (experimental!).

        This is useful when running a non-sequential list of task on a gdir,
        mostly for e.g. different experiments with different output files.

        Parameters
        ----------
        gdir : :py:class:`oggm.GlacierDirectory`
             the directory to process.
        tasks : list
             the the list of entity tasks to apply.
             Optionally, each list element can be a tuple, with the first element
             being the task, and the second element a dict that
             will be passed to the task function as ``**kwargs``.
        """

        res = []

        for task in tasks:
            # execute_entity_task calls tolist on the gdir, no need to turn into list
            out = self.execute_entity_task(task, gdir)
            res.append(out)

        return res


    def gdir_from_prepro(self, entity, from_prepro_level=None,
                         prepro_border=None, prepro_rgi_version=None,
                         check_demo_glacier=False, base_url=None):

        if prepro_border is None:
            prepro_border = int(self.PARAMS['border'])
        if prepro_rgi_version is None:
            prepro_rgi_version = self.PARAMS['rgi_version']
        try:
            rid = entity.RGIId
        except AttributeError:
            rid = entity

        if check_demo_glacier and base_url is None:
            demo_id = utils.demo_glacier_id(rid)
            if demo_id is not None:
                rid = demo_id
                entity = demo_id
                base_url = self.utils.DEMO_GDIR_URL

        tar_base = self.utils.get_prepro_gdir(prepro_rgi_version, rid, prepro_border,
                                              from_prepro_level, base_url=base_url)
        from_tar = os.path.join(tar_base.replace('.tar', ''), rid + '.tar.gz')
        return utils.GlacierDirectory(self, entity, from_tar=from_tar)


    def init_glacier_regions(self, rgidf=None, *, reset=False, force=False,
                             from_prepro_level=None, prepro_border=None,
                             prepro_rgi_version=None, prepro_base_url=None,
                             from_tar=False, delete_tar=False,
                             use_demo_glaciers=None):
        """DEPRECATED: Initializes the list of Glacier Directories for this run.

        This is the very first task to do (always). If the directories are already
        available in the working directory, use them. If not, create new ones.

        Parameters
        ----------
        rgidf : GeoDataFrame or list of ids, optional for pre-computed runs
            the RGI glacier outlines. If unavailable, OGGM will parse the
            information from the glacier directories found in the working
            directory. It is required for new runs.
        reset : bool
            delete the existing glacier directories if found.
        force : bool
            setting `reset=True` will trigger a yes/no question to the user. Set
            `force=True` to avoid this.
        from_prepro_level : int
            get the gdir data from the official pre-processed pool. See the
            documentation for more information
        prepro_border : int
            for `from_prepro_level` only: if you want to override the default
            behavior which is to use `oggm.PARAMS['border']`
        prepro_rgi_version : str
            for `from_prepro_level` only: if you want to override the default
            behavior which is to use `oggm.PARAMS['rgi_version']`
        prepro_base_url : str
            for `from_prepro_level` only: if you want to override the default
            URL from which to download the gdirs. Default currently is
            https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.1/
        use_demo_glaciers : bool
            whether to check the demo glaciers for download (faster than the
            standard prepro downloads). The default is to decide whether or
            not to check based on simple criteria such as glacier list size.
        from_tar : bool, default=False
            extract the gdir data from a tar file. If set to `True`,
            will check for a tar file at the expected location in `base_dir`.
        delete_tar : bool, default=False
            delete the original tar file after extraction.
        delete_tar : bool, default=False
            delete the original tar file after extraction.

        Returns
        -------
        gdirs : list of :py:class:`oggm.GlacierDirectory` objects
            the initialised glacier directories

        Notes
        -----
        This task is deprecated in favor of the more explicit
        init_glacier_directories. Indeed, init_glacier_directories is very
        similar to init_glacier_regions, but it does not process the DEMs:
        a glacier directory is valid also without DEM.
        """

        if reset and not force:
            reset = utils.query_yes_no('Delete all glacier directories?')

        if prepro_border is None:
            prepro_border = int(self.PARAMS['border'])

        if from_prepro_level and prepro_border not in [10, 80, 160, 250]:
            if 'test' not in self.utils.GDIR_URL:
                raise InvalidParamsError("prepro_border or oggm.PARAMS['border'] "
                                         "should be one of: 10, 80, 160, 250.")

        # if reset delete also the log directory
        if reset:
            fpath = os.path.join(self.PATHS['working_dir'], 'log')
            if os.path.exists(fpath):
                rmtree(fpath)

        gdirs = []
        new_gdirs = []
        if rgidf is None:
            if reset:
                raise ValueError('Cannot use reset without setting rgidf')
            log.workflow('init_glacier_regions by parsing available folders '
                         '(can be slow).')
            # The dirs should be there already
            gl_dir = os.path.join(self.PATHS['working_dir'], 'per_glacier')
            for root, _, files in os.walk(gl_dir):
                if files and ('dem.tif' in files):
                    gdirs.append(utils.GlacierDirectory(self, os.path.basename(root)))
        else:

            # Check if dataframe or list of strs
            try:
                entities = []
                for _, entity in rgidf.iterrows():
                    entities.append(entity)
            except AttributeError:
                entities = utils.tolist(rgidf)

            # Check demo
            if use_demo_glaciers is None:
                use_demo_glaciers = len(entities) < 100

            if from_prepro_level is not None:
                log.workflow('init_glacier_regions from prepro level {} on '
                             '{} glaciers.'.format(from_prepro_level,
                                                   len(entities)))
                # Read the hash dictionary before we use multiproc
                if self.PARAMS['dl_verify']:
                    self.utils.get_dl_verify_data('cluster.klima.uni-bremen.de')
                gdirs = self.execute_entity_task(gdir_from_prepro, entities,
                                                 from_prepro_level=from_prepro_level,
                                                 prepro_border=prepro_border,
                                                 prepro_rgi_version=prepro_rgi_version,
                                                 check_demo_glacier=use_demo_glaciers,
                                                 base_url=prepro_base_url)
            else:
                # We can set the intersects file automatically here
                if (self.PARAMS['use_intersects'] and
                        len(self.PARAMS['intersects_gdf']) == 0):
                    rgi_ids = np.unique(np.sort([entity.RGIId for entity in
                                                 entities]))
                    rgi_version = rgi_ids[0].split('-')[0][-2:]
                    fp = self.utils.get_rgi_intersects_entities(rgi_ids,
                                                                version=rgi_version)
                    self.set_intersects_db(fp)

                gdirs = self.execute_entity_task(utils.GlacierDirectory, entities,
                                                 reset=reset,
                                                 from_tar=from_tar,
                                                 delete_tar=delete_tar)

                for gdir in gdirs:
                    if not os.path.exists(gdir.get_filepath('dem')):
                        new_gdirs.append(gdir)

        if len(new_gdirs) > 0:
            # If not initialized, run the task in parallel
            self.execute_entity_task(tasks.define_glacier_region, new_gdirs)

        return gdirs


    def init_glacier_directories(self, rgidf=None, *, reset=False, force=False,
                                 from_prepro_level=None, prepro_border=None,
                                 prepro_rgi_version=None, prepro_base_url=None,
                                 from_tar=False, delete_tar=False,
                                 use_demo_glaciers=None):
        """Initializes the list of Glacier Directories for this run.

        This is the very first task to do (always). If the directories are already
        available in the working directory, use them. If not, create new ones.

        Parameters
        ----------
        rgidf : GeoDataFrame or list of ids, optional for pre-computed runs
            the RGI glacier outlines. If unavailable, OGGM will parse the
            information from the glacier directories found in the working
            directory. It is required for new runs.
        reset : bool
            delete the existing glacier directories if found.
        force : bool
            setting `reset=True` will trigger a yes/no question to the user. Set
            `force=True` to avoid this.
        from_prepro_level : int
            get the gdir data from the official pre-processed pool. See the
            documentation for more information
        prepro_border : int
            for `from_prepro_level` only: if you want to override the default
            behavior which is to use `oggm.PARAMS['border']`
        prepro_rgi_version : str
            for `from_prepro_level` only: if you want to override the default
            behavior which is to use `oggm.PARAMS['rgi_version']`
        prepro_base_url : str
            for `from_prepro_level` only: if you want to override the default
            URL from which to download the gdirs. Default currently is
            https://cluster.klima.uni-bremen.de/~fmaussion/gdirs/oggm_v1.1/
        use_demo_glaciers : bool
            whether to check the demo glaciers for download (faster than the
            standard prepro downloads). The default is to decide whether or
            not to check based on simple criteria such as glacier list size.
        from_tar : bool, default=False
            extract the gdir data from a tar file. If set to `True`,
            will check for a tar file at the expected location in `base_dir`.
        delete_tar : bool, default=False
            delete the original tar file after extraction.

        Returns
        -------
        gdirs : list of :py:class:`oggm.GlacierDirectory` objects
            the initialised glacier directories

        Notes
        -----
        This task is very similar to init_glacier_regions, with one main
        difference: it does not process the DEMs for this glacier.
        Eventually, init_glacier_regions will be deprecated and removed from the
        codebase.
        """

        if reset and not force:
            reset = utils.query_yes_no('Delete all glacier directories?')

        if prepro_border is None:
            prepro_border = int(self.PARAMS['border'])

        if from_prepro_level and prepro_border not in [10, 80, 160, 250]:
            if 'test' not in self.utils.GDIR_URL:
                raise InvalidParamsError("prepro_border or oggm.PARAMS['border'] "
                                         "should be one of: 10, 80, 160, 250.")

        # if reset delete also the log directory
        if reset:
            fpath = os.path.join(self.PATHS['working_dir'], 'log')
            if os.path.exists(fpath):
                rmtree(fpath)

        if rgidf is None:
            # Infer the glacier directories from folders available in working dir
            if reset:
                raise ValueError('Cannot use reset without setting rgidf')
            log.workflow('init_glacier_directories by parsing all available '
                         'folders (this takes time: if possible, provide rgidf '
                         'instead).')
            # The dirs should be there already
            gl_dir = os.path.join(self.PATHS['working_dir'], 'per_glacier')
            gdirs = []
            for root, _, files in os.walk(gl_dir):
                if files and ('outlines.shp' in files or
                              'outlines.tar.gz' in files):
                    gdirs.append(utils.GlacierDirectory(self, os.path.basename(root)))
        else:
            # TODO: prepare for rgidf to be dask.delayed!

            # Create glacier directories from input
            # Check if dataframe or list of str
            try:
                entities = []
                for _, entity in rgidf.iterrows():
                    entities.append(entity)
            except AttributeError:
                entities = utils.tolist(rgidf)

            # Check demo
            if use_demo_glaciers is None:
                use_demo_glaciers = len(entities) < 100

            if from_prepro_level is not None:
                log.workflow('init_glacier_directories from prepro level {} on '
                             '{} glaciers.'.format(from_prepro_level,
                                                   len(entities)))
                gdirs = self.execute_entity_task(gdir_from_prepro, entities,
                                                 from_prepro_level=from_prepro_level,
                                                 prepro_border=prepro_border,
                                                 prepro_rgi_version=prepro_rgi_version,
                                                 check_demo_glacier=use_demo_glaciers,
                                                 base_url=prepro_base_url)
            else:
                # We can set the intersects file automatically here
                if (self.PARAMS['use_intersects'] and
                        len(self.PARAMS['intersects_gdf']) == 0):
                    try:
                        rgi_ids = np.unique(np.sort([entity.RGIId for entity in
                                                     entities]))
                        rgi_version = rgi_ids[0].split('-')[0][-2:]
                        fp = self.utils.get_rgi_intersects_entities(rgi_ids,
                                                                    version=rgi_version)
                        self.set_intersects_db(fp)
                    except AttributeError:
                        # List of str
                        pass

                gdirs = self.execute_entity_task(utils.GlacierDirectory, entities,
                                                 reset=reset,
                                                 from_tar=from_tar,
                                                 delete_tar=delete_tar)

        return gdirs


    def gis_prepro_tasks(self, gdirs):
        """Shortcut function: run all flowline preprocessing tasks.

        Parameters
        ----------
        gdirs : list of :py:class:`oggm.GlacierDirectory` objects
            the glacier directories to process
        """

        from oggm import tasks

        task_list = [
            tasks.define_glacier_region,
            tasks.glacier_masks,
            tasks.compute_centerlines,
            tasks.initialize_flowlines,
            tasks.compute_downstream_line,
            tasks.compute_downstream_bedshape,
            tasks.catchment_area,
            tasks.catchment_intersections,
            tasks.catchment_width_geom,
            tasks.catchment_width_correction
        ]
        for task in task_list:
            self.execute_entity_task(task, gdirs)


    def climate_tasks(self, gdirs):
        """Shortcut function: run all climate related tasks.

        Parameters
        ----------
        gdirs : list of :py:class:`oggm.GlacierDirectory` objects
            the glacier directories to process
        """

        # Process climate data
        self.execute_entity_task(tasks.process_climate_data, gdirs)

        # Then, calibration?
        if self.PARAMS['run_mb_calibration']:
            tasks.compute_ref_t_stars(gdirs)

        # Mustar and the apparent mass-balance
        self.execute_entity_task(tasks.local_t_star, gdirs)
        self.execute_entity_task(tasks.mu_star_calibration, gdirs)


    def inversion_tasks(self, gdirs):
        """Shortcut function: run all ice thickness inversion tasks.

        Parameters
        ----------
        gdirs : list of :py:class:`oggm.GlacierDirectory` objects
            the glacier directories to process
        """

        if self.PARAMS['use_kcalving_for_inversion']:
            # Differentiate between calving and non-calving glaciers
            gdirs_nc = []
            gdirs_c = []
            for gd in gdirs:
                if gd.is_tidewater:
                    gdirs_c.append(gd)
                else:
                    gdirs_nc.append(gd)

            if gdirs_nc:
                self.execute_entity_task(tasks.prepare_for_inversion, gdirs_nc)
                self.execute_entity_task(tasks.mass_conservation_inversion, gdirs_nc)
                self.execute_entity_task(tasks.filter_inversion_output, gdirs_nc)

            if gdirs_c:
                self.execute_entity_task(tasks.find_inversion_calving, gdirs_c)
        else:
            self.execute_entity_task(tasks.prepare_for_inversion, gdirs)
            self.execute_entity_task(tasks.mass_conservation_inversion, gdirs)
            self.execute_entity_task(tasks.filter_inversion_output, gdirs)


    def merge_glacier_tasks(self, gdirs, main_rgi_id=None, return_all=False, buffer=None,
                            **kwargs):
        """Shortcut function: run all tasks to merge tributaries to a main glacier

        Parameters
        ----------
        gdirs : list of :py:class:`oggm.GlacierDirectory`
            all glaciers, main and tributary. Preprocessed and initialised
        main_rgi_id: str
            RGI ID of the main glacier of interest. If None is provided merging
            will start based uppon the largest glacier
        return_all : bool
            if main_rgi_id is given and return_all = False: only the main glaicer
            is returned
            if main_rgi_is given and return_all = True, the main glacier and every
            remaining glacier from the initial gdirs list is returned, possible
            merged as well.
        buffer : float
            buffer around a flowline to first better find an overlap with another
            flowline. And second assure some distance between the lines at a
            junction. Will default to `oggm.PARAMS['kbuffer']`.
        kwargs: keyword argument for the recursive merging

        Returns
        -------
        merged_gdirs: list of all merged :py:class:`oggm.GlacierDirectory`
        """

        from oggm.core import flowline

        if len(gdirs) > 100:
            raise InvalidParamsError('this could take time! I should include an '
                                     'optional parameter to ignore this.')

        # sort all glaciers descending by area
        gdirs.sort(key=lambda x: x.rgi_area_m2, reverse=True)

        # if main glacier is asked, put it in first position
        if main_rgi_id is not None:
            gdir_main = [gd for gd in gdirs if gd.rgi_id == main_rgi_id][0]
            gdirs.remove(gdir_main)
            gdirs = [gdir_main] + gdirs

        merged_gdirs = []
        while len(gdirs) > 1:
            # main glacier is always the first: either given or the largest one
            gdir_main = gdirs.pop(0)
            gdir_merged, gdirs = self._recursive_merging(gdirs, gdir_main, **kwargs)
            merged_gdirs.append(gdir_merged)

        # now we have gdirs which contain all the necessary flowlines,
        # time to clean them up
        for gdir in merged_gdirs:
            flowline.clean_merged_flowlines(gdir, buffer=buffer)

        if main_rgi_id is not None and return_all is False:
            return [gd for gd in merged_gdirs if main_rgi_id in gd.rgi_id][0]

        # add the remaining glacier to the final list
        merged_gdirs = merged_gdirs + gdirs

        return merged_gdirs


    def _recursive_merging(self, gdirs, gdir_main, glcdf=None,
                           filename='climate_historical', input_filesuffix=''):
        """ Recursive function to merge all tributary glaciers.

        This function should start with the largest glacier and then be called
        upon all smaller glaciers.

        Parameters
        ----------
        gdirs : list of :py:class:`oggm.GlacierDirectory`
            all glaciers, main and tributary. Preprocessed and initialised
        gdir_main: :py:class:`oggm.GlacierDirectory`
            the current main glacier where the others are merge to
        glcdf: geopandas.GeoDataFrame
            which contains the main glaciers, will be downloaded if None
        filename: str
            Baseline climate file
        input_filesuffix: str
            Filesuffix to the climate file

        Returns
        -------
        merged_gdir: :py:class:`oggm.GlacierDirectory`
            the mergeed current main glacier
        gdirs : list of :py:class:`oggm.GlacierDirectory`
            updated list of glaciers, removed the already merged ones
        """
        # find glaciers which intersect with the main

        from oggm.core import centerlines, flowline

        tributaries = centerlines.intersect_downstream_lines(gdir_main,
                                                             candidates=gdirs)
        if len(tributaries) == 0:
            # if no tributaries: nothing to do
            return gdir_main, gdirs

        # seperate those glaciers which are not already found to be a tributary
        gdirs = [gd for gd in gdirs if gd not in tributaries]

        gdirs_to_merge = []

        for trib in tributaries:
            # for each tributary: check if we can merge additional glaciers to it
            merged, gdirs = self._recursive_merging(gdirs, trib, glcdf=glcdf,
                                                    filename=filename,
                                                    input_filesuffix=input_filesuffix)
            gdirs_to_merge.append(merged)

        # create merged glacier directory
        gdir_merged = utils.initialize_merged_gdir(
            gdir_main, tribs=gdirs_to_merge, glcdf=glcdf, filename=filename,
            input_filesuffix=input_filesuffix)

        flowline.merge_to_one_glacier(gdir_merged, gdirs_to_merge,
                                      filename=filename,
                                      input_filesuffix=input_filesuffix)

        return gdir_merged, gdirs
