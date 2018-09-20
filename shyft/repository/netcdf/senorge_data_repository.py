# This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
# See file COPYING for more details **/
# questions: sven.decker@geo.uio.no
import re
import os
import numpy as np
from netCDF4 import Dataset
from shyft import api
from shyft import shyftdata_dir
from .. import interfaces
from .time_conversion import convert_netcdf_time
from .utils import _limit_2D, _slice_var_2D, _numpy_to_geo_ts_vec, dummy_var, create_geo_ts_type_map

class SeNorgeDataRepositoryError(Exception):
    pass

class SeNorgeDataRepository(interfaces.GeoTsRepository):
    """
    Repository for senorge TEMP and PREC netcdf data files available per:


    Cristian Lussana, & Ole Einar Tveito. (2017). seNorge2 dataset [Data set]. Zenodo.
    http://doi.org/10.5281/zenodo.845733

    """

    def __init__(self, epsg, directory, filename=None, elevation_file=None, padding=5000., allow_subset=False):
        """
        Construct the netCDF4 dataset reader for data for seNorge v2 datasets,

        Information on datasets are available from: http://doi.org/10.5281/zenodo.845733

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates.
            Currently "32632" and "32633" are supported.
        directory: string
            Path to directory holding the data files.
            os.path.isdir(directory) should be true, or exception is raised.
        filename: string, optional
            Name of netcdf file in directory that contains spatially
            distributed input data. Can be a regex pattern as well
        padding: float, optional
            padding in meters
        allow_subset: bool
            Allow extraction of a subset of the given source fields
            instead of raising exception.
        """
        directory = os.path.expandvars(directory)
        self._directory = directory
        self._filename = filename
        self.allow_subset = allow_subset
        if not os.path.isdir(directory):
            raise SeNorgeDataRepositoryError("No such directory '{}'".format(directory))

        if elevation_file is not None:
            self.elevation_file = os.path.join(self._directory, elevation_file)
            if not os.path.isfile(self.elevation_file):
                raise SeNorgeDataRepositoryError(
                    "Elevation file '{}' not found".format(self.elevation_file))
        else:
            self.elevation_file = os.path.join(shyftdata_dir, "repository/senorge_data_repository/seNorge2_dem_UTM33.nc")

        self.shyft_cs = "+init=EPSG:{}".format(epsg)
        self._padding = padding

        # Field names and mappings
        self.senorge_shyft_map = {
            "mean_temperature": "temperature",
            "precipitation_amount": "precipitation",
            "radiation": "radiation",
            "wind_speed": "wind_speed",
            "relative_humidity": "relative_humidity"
        }
        self.series_type = {
            "relative_humidity": api.POINT_AVERAGE_VALUE,
            "temperature": api.POINT_AVERAGE_VALUE,
            "precipitation": api.POINT_AVERAGE_VALUE,
            "radiation": api.POINT_AVERAGE_VALUE,
            "wind_speed": api.POINT_AVERAGE_VALUE
        }

    def get_timeseries(self, input_source_types, utc_period, geo_location_criteria=None):
        """
        see shyft.repository.interfaces.GeoTsRepository
        """

        # for these variables we use 'dummy' data
        sources = {}
        input_source_types_dummy = \
            [ist for ist in input_source_types if ist in ('radiation',
                                                          'wind_speed',
                                                          'relative_humidity')]  # for these variables we use 'dummy' data
        input_source_types_real = \
            [ist for ist in input_source_types if ist in ('temperature',
                                                          'precipitation')]
        input_source_types_remaining = [ist for ist in input_source_types if
                                   ist not in input_source_types_dummy and ist not in input_source_types_real]
        if input_source_types_remaining:
            raise SeNorgeDataRepositoryError("Input source types {} not supported".format(input_source_types_remaining))
        if input_source_types_dummy:
            sources_dummy = dummy_var(input_source_types_dummy, utc_period, geo_location_criteria)
            sources.update(sources_dummy)
        if input_source_types_real:
            filename = os.path.join(self._directory, self._filename)
            if not os.path.isfile(filename):
                #if re.compile(self._filename).groups > 0:  # check if it is a filename-pattern
                #    filename = _get_files(self._directory, self._filename, utc_period.start, SeNorgeDataRepositoryError)
                #else:
                raise SeNorgeDataRepositoryError("File '{}' not found".format(filename))
            with Dataset(filename) as dataset:
                sources_real = self._get_data_from_dataset(dataset, input_source_types_real, utc_period, geo_location_criteria)
            sources.update(sources_real)
        return sources

    def _check_and_get_coord_vars(self, dataset, var_types):
        cs = []
        coord_names = []
        for k, v in self.senorge_shyft_map.items():
            if v in var_types and k in dataset.variables:
                cs.append(dataset.variables[k].getncattr('grid_mapping'))
                coord_names.append([d for d in dataset.variables[k].dimensions if d in ['time', 'X', 'Y', 'latitude', 'longitude']])
        if not all(elem == cs[0] for elem in cs):
            SeNorgeDataRepositoryError('Requested vars have different coord_sys. Do index extraction per var.')
        if not all(elem == coord_names[0] for elem in coord_names):
            SeNorgeDataRepositoryError('Requested vars have different coords. Do index extraction per var.')
        time = dataset.variables.get("time", None)
        if not time:
            raise SeNorgeDataRepositoryError("Time variable not found in dataset.")
        time = convert_netcdf_time(time.units, time)

        if 'Y' in coord_names[0]:
            x = dataset.variables.get("X", None)
            y = dataset.variables.get("Y", None)
        elif 'latitude' in coord_names[0]:
            x = dataset.variables.get("longitude", None)
            y = dataset.variables.get("latitude", None)
        else:
            SeNorgeDataRepositoryError('No recognized coordinate dimension names found.')

        if not all([x, y]):
            raise SeNorgeDataRepositoryError("Spatial Coordinate variables not found in dataset.")
        if 'Y' in coord_names[0]:
            if not all([var.units in ['km', 'm', 'meters'] for var in [x, y]]) and x.units == y.units:
                raise SeNorgeDataRepositoryError("The unit for x and y coordinates should be either m or km.")
        else:
            if not (y.units == 'degrees_north' and x.units == 'degrees_east'):
                raise SeNorgeDataRepositoryError("The unit for latitude and longitude coordinates should be "
                                               "'degrees_north' and 'degrees_east' repectively.")
        coord_conv = 1.
        if y.units == 'km':
            coord_conv = 1000.

        data_cs = dataset.variables.get(cs[0], None)
        if data_cs is None:
            raise SeNorgeDataRepositoryError("No coordinate system information in dataset.")
        return time, x, y, data_cs, coord_conv


    def _get_data_from_dataset(self, dataset, input_source_types, utc_period, geo_location_criteria, ensemble_member=None):

        #copied from met_netcdf repo
        # Check for presence and consistency of coordinate variables
        time, x_var, y_var, data_cs, coord_conv = self._check_and_get_coord_vars(dataset, input_source_types)

        x, y, (x_inds, y_inds), (x_slice, y_slice) = _limit_2D(x_var[:] * coord_conv, y_var[:] * coord_conv,
                                                               data_cs.proj4, self.shyft_cs, geo_location_criteria,
                                                               self._padding, SeNorgeDataRepositoryError)

        # "The daily precipitation (comment: and temperature) for day D has been defined as the accumulated
        # precipitation between 06:00 UTC of day D − 1  and 06:00 UTC of day D".
        # https://www.earth-syst-sci-data.net/10/235/2018/essd-10-235-2018.pdf
        # Here: time points are interpreted in datasert as UTC 00:00 of day D.
        #time = time - 18. * 3600.  # Correcting dataset time to true utc time (18.00 hours earlier than
                                   # stated in dataset, e.g. 1957-01-01 00:00 -> 1956-12-31 06:00).
        time = time + 6. * 3600.  #
        # Make temporal slilce
        time_slice, issubset = self._make_time_slice(time, utc_period, SeNorgeDataRepositoryError)


        raw_data = {}
        for k in dataset.variables.keys():

            if k in self.senorge_shyft_map.keys():
                var = self.senorge_shyft_map.get(k, None)
                if var in input_source_types:
                # if k in self._shift_fields and issubset:  # Add one to time slice
                #     data_time_slice = slice(time_slice.start, time_slice.stop + 1)  # to supply the extra value that is needed for accumulated variables
                # else:
                #     data_time_slice = time_slice
                    data = dataset.variables[k]
                    pure_arr = _slice_var_2D(data, x_var.name, y_var.name, x_slice, y_slice, x_inds, y_inds, SeNorgeDataRepositoryError,
                                         # slices={'time': data_time_slice, 'ensemble_member': ensemble_member})
                                         slices={'time': time_slice, 'ensemble_member': ensemble_member})
                    raw_data[self.senorge_shyft_map[k]] = pure_arr, k, data.units

        if self.elevation_file is not None:
            _x, _y, z = self._read_elevation_file(self.elevation_file, x_var.name, y_var.name,
                                                  geo_location_criteria)
            assert np.linalg.norm(x - _x) < 1.0e-10  # x/y coordinates should match
            assert np.linalg.norm(y - _y) < 1.0e-10
        else:
            raise SeNorgeDataRepositoryError("No elevation file given.")

        # Make sure requested fields are valid, and that dataset contains the requested data.

        if not self.allow_subset and not (set(raw_data.keys()).issuperset(input_source_types)):
            raise SeNorgeDataRepositoryError("Could not find all data fields")

        if ensemble_member is None and 'ensemble_member' in data.dimensions:
            dims_flat = [d for d in data.dimensions if d in ['time', 'ensemble_member', x_var.name]]
            ens_dim_idx = dims_flat.index('ensemble_member')
            ens_slice = len(dims_flat) * [slice(None)]
            returned_data = []
            for i in range(dataset.dimensions['ensemble_member'].size):
                ens_slice[ens_dim_idx] = i
                ensemble_raw = {k: (raw_data[k][0][ens_slice], raw_data[k][1], raw_data[k][2]) for k in raw_data.keys()}
                returned_data.append(_numpy_to_geo_ts_vec(self._transform_raw(ensemble_raw, time[time_slice]),#, issubset=issubset)
                                                          x, y, z, SeNorgeDataRepositoryError))
        else:
            returned_data = _numpy_to_geo_ts_vec(self._transform_raw(raw_data, time[time_slice]),#, issubset=issubset),
                                                 x, y, z, SeNorgeDataRepositoryError)
            data = self._transform_raw(raw_data, time[time_slice])
            geo_pts = api.GeoPointVector.create_from_x_y_z(*[api.DoubleVector_FromNdArray(arr) for arr in [x, y, z]])
            shape = list(data.values())[0][0].shape
            ndim = len(shape)
            geo_ts = {key: create_geo_ts_type_map[key](ta, geo_pts, arr[:, :].transpose(), self.series_type[key])
                      for key, (arr, ta) in data.items()}
        return returned_data

    def _make_time_slice(self, time, utc_period, err):
        UTC = api.Calendar()
        del_t = int(time[1]-time[0])
        if utc_period is None:
            # If period is None set period to be from start to end of dataset time dimension
            utc_period = api.UtcPeriod(int(time[0]), int(time[-1])+del_t)
        idx_min = np.argmin(time-del_t <= utc_period.start) - 1  # raise error if result is -1
        idx_max = np.argmax(time >= utc_period.end)  # raise error if result is 0
        if idx_min < 0:
            raise err(
                "The earliest time in repository ({}) is later than the start of the period for which data is "
                "requested ({})".format(UTC.to_string(int(time[0]-del_t)), UTC.to_string(utc_period.start)))
        if idx_max == 0:
            raise err(
                "The latest time in repository ({}) is earlier than the end of the period for which data is "
                "requested ({})".format(UTC.to_string(int(time[-1])), UTC.to_string(utc_period.end)))

        issubset = True if idx_max < len(time) - 1 else False
        time_slice = slice(idx_min, idx_max + 1)
        return time_slice, issubset


    def _read_elevation_file(self, filename, x_var_name, y_var_name, geo_location_criteria):
        with Dataset(filename) as dataset:
            elev = dataset.variables["elevation"]
            if "elevation" not in dataset.variables.keys():
                raise interfaces.InterfaceError(
                    "File '{}' does not contain altitudes".format(filename))
            x, y, (x_inds, y_inds), (x_slice, y_slice) = _limit_2D(dataset.variables['easting'][:], dataset.variables['northing'][:],
                                                                  elev.projection, self.shyft_cs, geo_location_criteria,
                                                                   self._padding, SeNorgeDataRepositoryError)
            x_var_name = 'easting'
            y_var_name = 'northing'
            z = _slice_var_2D(elev, x_var_name, y_var_name, x_slice, y_slice, x_inds, y_inds, SeNorgeDataRepositoryError)
            return x, y, z


    def _transform_raw(self, data, time, issubset=False):
        """
        We need full time if deaccumulating
        """

        def dacc_time(t):
            del_t = t[1] - t[0]
            t0 = t[0]-del_t
            return api.TimeAxis(t0, del_t, len(t))

        def noop_space(x):
            return x

        def air_temp_conv(T):
            return T

        def prec_conv(p):
            return p/24. # mm/day -> mm/h

        convert_map = {"mean_temperature": lambda x, t, u: (air_temp_conv(x), dacc_time(t)),
                       "precipitation_amount": lambda x, t, u: (prec_conv(x), dacc_time(t))}

        res = {}
        for k, (v, ak, unit) in data.items():
            res[k] = convert_map[ak](v, time, unit)
        return res


