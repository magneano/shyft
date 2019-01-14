# This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
# See file COPYING for more details **/
import os
import time
from shyft import api
from netCDF4 import Dataset
from .time_conversion import convert_netcdf_time
from shyft.repository.interfaces import GeoTsRepository, ForecastSelectionCriteria
from shyft.repository.netcdf.concat_data_repository import ConcatDataRepository
from shyft.repository.netcdf.met_netcdf_data_repository import MetNetcdfDataRepository
# import numpy as np
import datetime
from shyft.repository.netcdf.utils import _clip_ensemble_of_geo_timeseries
from shyft.repository.netcdf.utils import parallelize_geo_timeseries
from shyft.repository.netcdf.utils import merge_ensemble_geo_timeseries

UTC = api.Calendar()

class WXRepositoryError(Exception):
    pass

class WXRepository(GeoTsRepository):

    def __init__(self, epsg, filename, padding=15000., flattened=True, allow_year_shift=True, cache_data=True):
        """
        Repository for reading ensemble weather scenarios from netcdf file.

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates. Currently "32632" and "32633" are supported.
        filename: string
            Path to netcdf file containing weather data
        flattened: bool
            Flags whether grid_points are flattened
        allow_year_shift: bool
            Flags whether shift of years is allowed
        cache_data: bool
            Use cache data if True
        """
        self.allow_year_shift = allow_year_shift
        self.cache_data = cache_data
        self.cache = None
        if flattened:
            self.wx_repo = ConcatDataRepository(epsg, filename, padding=padding)
        elif not flattened:
            self.wx_repo = MetNetcdfDataRepository(epsg, None, filename, padding=padding)
            filename = os.path.expandvars(filename)
            with Dataset(filename) as dataset:
                time = dataset.variables.get("time", None)
                time = convert_netcdf_time(time.units, time)
                self.wx_repo.time = time

        self.source_type_map = {"relative_humidity": api.RelHumSource,
                                "temperature": api.TemperatureSource,
                                "precipitation": api.PrecipitationSource,
                                "radiation": api.RadiationSource,
                                "wind_speed": api.WindSpeedSource}

        self.source_vector_map = {"relative_humidity": api.RelHumSourceVector,
                                "temperature": api.TemperatureSourceVector,
                                "precipitation": api.PrecipitationSourceVector,
                                "radiation": api.RadiationSourceVector,
                                "wind_speed": api.WindSpeedSourceVector}

    def get_forecast_ensemble(self, input_source_types, utc_period, t_c, geo_location_criteria=None):
        """
        Same as get_timeseries since no time_stamp structure to filename

        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        print("WXRepository.get_forecast_ensembles")
        t_total = time.time()
        res = self.get_timeseries_ensemble(input_source_types, utc_period, geo_location_criteria=geo_location_criteria)
        elapsed_time = time.time() - t_total
        print("Total time for WXRepository.get_forecast_ensembles: {}".format(elapsed_time))
        return res

    def get_timeseries_ensemble(self, input_source_types, utc_period, geo_location_criteria=None):
        """
            Get ensemble of shyft source vectors of time series covering utc_period
            for input_source_types.

            Parameters
            ----------
            see interfaces.GeoTsRepository

            Returns
            -------
            see interfaces.GeoTsRepository
        """
        wx_repo = self.wx_repo
        if self.allow_year_shift and utc_period is not None:
            utc_start_date = datetime.datetime.utcfromtimestamp(utc_period.start)
            repo_date_start = datetime.datetime.utcfromtimestamp(int(wx_repo.time[0]))
            if utc_start_date.timetuple().tm_yday + utc_start_date.timetuple().tm_hour >= \
                    repo_date_start.timetuple().tm_yday + repo_date_start.timetuple().tm_hour:
                utc_start_shifted_year = repo_date_start.timetuple().tm_year
            else:
                utc_start_shifted_year = repo_date_start.timetuple().tm_year + 1
            utc_start_shifted = UTC.time(utc_start_shifted_year, utc_start_date.month, utc_start_date.day,  utc_start_date.hour)
            d_t = utc_period.start - utc_start_shifted
            utc_end_shifted = utc_period.end - d_t
            utc_period_shifted = api.UtcPeriod(utc_start_shifted, utc_end_shifted)
        else:
            d_t = 0
            utc_period_shifted = utc_period
        if self.cache_data:
            if self.cache is None:
                self.cache = wx_repo.get_timeseries_ensemble(input_source_types, None, geo_location_criteria)
            raw_ens = self.cache
        else:
            raw_ens = wx_repo.get_timeseries_ensemble(input_source_types, utc_period_shifted, geo_location_criteria)
        res = [{key: self.source_vector_map[key]([self.source_type_map[key](src.mid_point(), src.ts.time_shift(d_t))
                    for src in geo_ts]) for key, geo_ts in ens.items()} for ens in raw_ens]
        return _clip_ensemble_of_geo_timeseries(res, utc_period, WXRepositoryError)


class WXParallelizationRepositoryError(Exception):
    pass

class WXParallelizationRepository(GeoTsRepository):

    def __init__(self, epsg, filename, truth_file=None, padding=15000., cache_data=True, numb_years=None):
        """
        Reads weather scenario from file. Get_timeseire_ensemble and get_forecast_ensemble returns
        ensemble of weather scenarios using parallelization.

        TODO: describe requirements for netcfdf-file

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates. Currently "32632" and "32633" are supported.
        filename: string
            Path to netcdf file containing weather data
        truth_file: file_path
            File path to netcdf file containing truth scenario. If not None, truth scenario
            will be concatenated onto start of ensemble. Note that netcdf file containing
            truth scenarios must have same format as synthetic scenarios
        cache_data: bool
            Use cache data if True
        numb_years: int
            Limits number of years to return. If None return as many as possible.
        """
        self.cache_data = cache_data
        self.cache = None
        self.truth_file = truth_file
        self.numb_years = numb_years
        self.syn_repo = WXRepository(epsg, filename, padding=padding, flattened=True, allow_year_shift=True, cache_data=False)
        if truth_file is not None:
            self.truth_repo = WXRepository(epsg, truth_file, padding=padding, flattened=True, allow_year_shift=False, cache_data=False)

    def get_timeseries_ensemble(self, input_source_types, utc_period, geo_location_criteria=None):
        """
             Get ensemble of shyft source vectors of time series covering utc_period
             for input_source_types.

             Parameters
             ----------
             see interfaces.GeoTsRepository

             Returns
             -------
             see interfaces.GeoTsRepository
        """

        if self.cache is None:
            if self.truth_file is not None:
                truth_end_time = self.truth_repo.wx_repo.time[0] + self.truth_repo.wx_repo.lead_times_in_sec[-1]
                # Shift time period for synthetic scenario so they overlap truth scenario
                syn_start_time = self.syn_repo.wx_repo.time[0]
                syn_end_time = self.syn_repo.wx_repo.time[0] + self.syn_repo.wx_repo.lead_times_in_sec[-1]
                truth_end_date = datetime.datetime.utcfromtimestamp(int(truth_end_time))
                syn_start_date = datetime.datetime.utcfromtimestamp(int(syn_start_time))
                if truth_end_date.timetuple().tm_yday + truth_end_date.timetuple().tm_hour >= \
                        syn_start_date.timetuple().tm_yday + syn_start_date.timetuple().tm_hour:
                    syn_start_shifted_year = truth_end_date.timetuple().tm_year
                else:
                    syn_start_shifted_year = truth_end_date.timetuple().tm_year - 1
                syn_start_time_shifted = UTC.time(syn_start_shifted_year, syn_start_date.month, syn_start_date.day,
                                                  syn_start_date.hour)
                d_t = syn_start_time - syn_start_time_shifted
                syn_end_time_shifted = syn_end_time - d_t
                syn_period_shifted = api.UtcPeriod(syn_start_time_shifted, syn_end_time_shifted)
                # get (shifted) time series for synthetic series
                scen = self.syn_repo.get_timeseries_ensemble(input_source_types, syn_period_shifted, geo_location_criteria=geo_location_criteria)
                # get truth and merge
                truth_scen = self.truth_repo.get_timeseries_ensemble(input_source_types, None, geo_location_criteria=geo_location_criteria)
                scen = merge_ensemble_geo_timeseries(truth_scen, scen)
            else:
                scen = self.syn_repo.get_timeseries_ensemble(input_source_types, None, geo_location_criteria=geo_location_criteria)
            if self.cache_data:
                self.cache = scen
        else:
            scen = self.cache
        self.parallelized_years, wx_scen_parallelized = parallelize_geo_timeseries(scen[0], utc_period, self.numb_years)
        return wx_scen_parallelized

    def get_forecast_ensemble(self, input_source_types, utc_period, t_c, geo_location_criteria=None):
        """
        Same as get_timeseries since no time_stamp structure to filename

        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        print("WXParallelizationRepository.get_forecast_ensembles")
        t_total = time.time()
        res = self.get_timeseries_ensemble(input_source_types, utc_period, geo_location_criteria=geo_location_criteria)
        elapsed_time = time.time() - t_total
        print("Total time for WXParallelizationRepository.get_forecast_ensembles: {}".format(elapsed_time))
        return res



