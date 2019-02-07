# This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
# See file COPYING for more details **/
import os
import time
from typing import Optional, List, Any, Dict
from shyft import api
from netCDF4 import Dataset
from .time_conversion import convert_netcdf_time
from shyft.repository.interfaces import GeoTsRepository, ForecastSelectionCriteria
from shyft.repository.netcdf.concat_data_repository import ConcatDataRepository
from shyft.repository.netcdf.met_netcdf_data_repository import MetNetcdfDataRepository
from shyft.repository.netcdf.utils import _clip_ensemble_of_geo_timeseries
from shyft.repository.netcdf.utils import parallelize_geo_timeseries
from shyft.repository.netcdf.utils import merge_ensemble_geo_timeseries

UTC = api.Calendar()

class WXRepositoryError(Exception):
    pass

class WXRepository(GeoTsRepository):

    def __init__(self, epsg: int, filename: str, padding: float = 15000., flattened: bool = True,
                 allow_year_shift: bool = True, cache_data: bool = True, ensemble_member: int = 0):
        """
        Repository for reading ensemble weather scenarios from netcdf file.

        Parameters
        ----------
        epsg: int
            Unique coordinate system id for result coordinates. Currently "32632" and "32633" are supported.
        filename: string
            Path to netcdf file containing weather data
        flattened: bool
            Flags whether grid_points are flattened
        allow_year_shift: bool
            Flags whether shift of years is allowed
        cache_data: bool
            Use cache data if True
        ensemble_member: int, optional
            Ensemble member returned by get_timeseries if dataset is of ensemble type (has dimension 'ensemble_member').
            Must be non-negative integer less than dimension size of 'ensemble_member'


        NetCDF4 dataset assumptions depends on boolean value of 'flattened':
            flattened = False:
                see requirements of MetNetcdfDataRepository
            flattened = True:
                see requirements of ConcatDataRepository with the additional restriction that time-dimension size is 1
        """
        self.allow_year_shift = allow_year_shift
        self.cache_data = cache_data
        self.ensemble_member = ensemble_member
        self.cache = None
        if flattened:
            self.wx_repo = ConcatDataRepository(epsg, filename, padding=padding, ensemble_member=ensemble_member)
        elif not flattened:
            self.wx_repo = MetNetcdfDataRepository(epsg, None, filename, padding=padding,
                                                   ensemble_member=ensemble_member)
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

    def get_forecast_ensemble(self, input_source_types: List[str], utc_period: Any, t_c: int,
                              geo_location_criteria: Optional[Any] = None) -> List[Dict[str, Any]]:

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

    def get_timeseries_ensemble(self, input_source_types: List[str], utc_period: Any,
                                geo_location_criteria: Optional[Any] = None) -> List[Dict[str, Any]]:
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
            utc_start_date = UTC.calendar_units(utc_period.start)
            repo_start_date = UTC.calendar_units(int(wx_repo.time[0]))
            if utc_period.start - UTC.time(utc_start_date.year, 1, 1) >= \
                    int(wx_repo.time[0]) - UTC.time(repo_start_date.year, 1, 1):
                utc_start_shifted_year = repo_start_date.year
            else:
                utc_start_shifted_year = repo_start_date.year + 1
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

    def get_timeseries(self, input_source_types: List[str], utc_period: Any,
                       geo_location_criteria: Optional[Any] = None) -> Dict[str, Any]:
        """
            Parameters
            ----------
            see interfaces.GeoTsRepository

            Returns
            -------
            see interfaces.GeoTsRepository
        """
        return self.get_timeseries_ensemble(input_source_types, utc_period, geo_location_criteria=None)[self.ensemble_member]

class WXParallelizationRepositoryError(Exception):
    pass

class WXParallelizationRepository(GeoTsRepository):

    def __init__(self, epsg: int, filename: str, truth_file: Optional[str] = None, padding: float = 15000.,
                 cache_data: bool = True, first_year: Optional[int] = None, numb_years: Optional[int] = None):
        """
        Reads weather scenario from flattened netcdf-file(s) of geo-located weather data.

        Get_timeseries_ensemble and get_forecast_ensemble returns ensemble of weather scenarios using parallelization.
        Get_timeseries returns a single scenario (first ensemble_member) in file

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates. Currently "32632" and "32633" are supported.
        filename: string
            Path to netcdf file containing synthetic weather date
        truth_file: file_path
            File path to netcdf file containing truth scenario. If not None, truth scenario
            will be concatenated onto start of synthetic scenario before parallelization. Note
            that netcdf file containing truth scenarios must have same format as synthetic scenarios
        cache_data: bool
            Use cache data if True
        first_year: int
            First year to start parallelization from. I None start as early as possible
        numb_years: int
            Limits number of years (i.e. scenarios) to return from ensemble methods. If None return as many as possible.


        NetCDF4 dataset assumptions:
            * Dimensions:
                * time (size = 1)
                * lead_time
                * ensemble_member (optional) - only first ensemble member is used
                * grid_point
            * Variables:
                * time:(float) array with periodic forecast creation timestamps in seconds since (1970.01.01 00:00, UTC)
                * lead_time:(float) array with hours since creation time
                * x: (float) array of latitudes of dims (grid_point)
                * y: (float) array of longitudes of dims (grid_point)
                * z: (float) array of altitudes [m] of dims (grid_point)
                * forecast_is_complete: flag array of dims (time)
                * crs: has attribute proj4, a string describing the coordinate system
            * Optional variables:
                * dew_point_temperature_2m: [K],
                * surface_air_pressure: [Pa],
                * relative_humidity_2m: [1],
                * air_temperature_2m: [K],
                * precipitation_amount: [kg/m^2],
                * precipitation_amount_acc: [kg/m^2],
                * x_wind_10m: [m/s],
                * y_wind_10m: [m/s],
                * windspeed_10m: [m/s],
                * integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time: [W s/m^2]}
                * All optional variables are (float) array with dims (time, lead_time, [ensemble_member], grid_point)
        """
        self.cache_data = cache_data
        self.cache = None
        self.truth_file = truth_file
        self.first_year = first_year
        self.numb_years = numb_years
        self.syn_repo = WXRepository(epsg, filename, padding=padding, flattened=True, allow_year_shift=True, cache_data=False)
        if truth_file is not None:
            self.truth_repo = WXRepository(epsg, truth_file, padding=padding, flattened=True, allow_year_shift=False, cache_data=False)

    def get_timeseries_ensemble(self, input_source_types: List[str], utc_period: Any,
                       geo_location_criteria: Optional[Any] = None) -> List[Dict[str, Any]]:
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
        scen = self._read_from_file_and_merge(input_source_types, geo_location_criteria=geo_location_criteria)
        self.parallelized_years, wx_scen_parallelized = parallelize_geo_timeseries(scen[0], utc_period,
                                                                first_year=self.first_year, numb_years=self.numb_years)
        return wx_scen_parallelized

    def get_forecast_ensemble(self, input_source_types: List[str], utc_period: Any, t_c: int,
                              geo_location_criteria: Optional[Any] = None) -> List[Dict[str, Any]]:

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

    def get_timeseries(self, input_source_types: List[str], utc_period: Any,
                       geo_location_criteria: Optional[Any] = None) -> Dict[str, Any]:
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
        scen = self._read_from_file_and_merge(input_source_types, geo_location_criteria=geo_location_criteria)
        return _clip_ensemble_of_geo_timeseries(scen, utc_period, WXParallelizationRepositoryError)[0]

    def _read_from_file_and_merge(self, input_source_types: List, geo_location_criteria: Optional[Any] = None) -> Any:
        if self.cache is None:
            if self.truth_file is not None:
                truth_end_time = self.truth_repo.wx_repo.time[0] + self.truth_repo.wx_repo.lead_times_in_sec[-1]
                # Shift time period for synthetic scenario so they overlap truth scenario
                syn_start_time = self.syn_repo.wx_repo.time[0]
                syn_end_time = self.syn_repo.wx_repo.time[0] + self.syn_repo.wx_repo.lead_times_in_sec[-1]
                truth_end_date = UTC.calendar_units(int(truth_end_time))
                syn_start_date =  UTC.calendar_units(int(syn_start_time))
                if int(truth_end_time) - UTC.time(truth_end_date.year, 1, 1) >= \
                        int(syn_start_time) - UTC.time(syn_start_date.year, 1, 1):
                    syn_start_shifted_year = truth_end_date.year
                else:
                    syn_start_shifted_year = truth_end_date.year - 1
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
        return scen



