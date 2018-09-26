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
import datetime as dt
from shyft.repository.netcdf.utils  import _clip_ensemble_of_geo_timeseries

UTC = api.Calendar()

class WXRepositoryError(Exception):
    pass

class WXRepository(GeoTsRepository):

    def __init__(self, epsg, filename, padding=15000., flattened=False, allow_year_shift=True, cache_data=True):
        """
        Construct the netCDF4 dataset reader for concatenated gridded forecasts and initialize data retrieval.

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates. Currently "32632" and "32633" are supported.
        filename: string
            Path to netcdf file containing concatenated forecasts
        flattened: bool
            Flags whether grid_points are flattened
        allow_year_shift: bool
            Flags whether shift of years is allowed
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
            utc_start_date = dt.datetime.utcfromtimestamp(utc_period.start)
            repo_date_start = dt.datetime.utcfromtimestamp(int(wx_repo.time[0]))
            utc_start_tt = utc_start_date.timetuple()
            if utc_start_tt.tm_yday + utc_start_tt.tm_sec >= repo_date_start.timetuple().tm_yday + repo_date_start.timetuple().tm_sec:
                utc_start_shifted_year = repo_date_start.timetuple().tm_year
            else:
                utc_start_shifted_year = repo_date_start.timetuple().tm_year + 1
            utc_start_shifted = UTC.time(utc_start_shifted_year, utc_start_tt.tm_mon, utc_start_tt.tm_mday,  utc_start_tt.tm_hour)
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