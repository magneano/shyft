import unittest
from os import path
from netCDF4 import Dataset
from pyproj import Proj
from pyproj import transform
import numpy as np

from shyft import shyftdata_dir
from shyft import api
from shyft.repository.netcdf.wx_repository import WXRepository, WXRepositoryError
from shyft.repository.netcdf.wx_repository import WXParallelizationRepository, WXParallelizationRepositoryError
from shapely.geometry import box


class WXParallelizationRepositoryTestCase(unittest.TestCase):
    def test_get_timeseries_ensemble(self):
        """
        Simple regression test of WX data repository.
        """
        EPSG, bbox, bpoly = self.wx_epsg_bbox

        # Period
        utc = api.Calendar()
        start = utc.time(2018, 10, 1)
        end = utc.time(2019, 10, 15)
        utc_period = api.UtcPeriod(start, end)

        # number of scenarios
        numb_years = 10

        f1 = path.join(shyftdata_dir, "repository", "wx_repository", f"flattened_nordic_scenarios_truth_sim_Vik.nc")

        wx_repo = WXParallelizationRepository(EPSG, filename=f1, numb_years=numb_years)
        source_types = ('temperature', 'precipitation', 'wind_speed', 'relative_humidity', 'radiation')
        wx_scenarios = wx_repo.get_timeseries_ensemble(source_types, utc_period, geo_location_criteria=bpoly)
        self.assertTrue(len(wx_scenarios) == numb_years)
        [self.assertTrue(set(scen) == set(source_types)) for scen in wx_scenarios]
        [self.assertTrue(wx_scenarios[0][src][0].ts.total_period().start <= start) for src in source_types]
        [self.assertTrue(wx_scenarios[0][src][0].ts.total_period().end >= end) for src in source_types]

        wx_scenarios_get_fcst = wx_repo.get_forecast_ensemble(source_types, utc_period, None, geo_location_criteria=bpoly)
        self.assertTrue(len(wx_scenarios_get_fcst) == numb_years)
        [self.assertTrue(set(scen) == set(source_types)) for scen in wx_scenarios_get_fcst]
        [self.assertTrue(wx_scenarios_get_fcst[0][src][0].ts.total_period().start <= start) for src in source_types]
        [self.assertTrue(wx_scenarios_get_fcst[0][src][0].ts.total_period().end >= end) for src in source_types]

        self.assertEqual(wx_scenarios[0]['temperature'][0].ts, wx_scenarios_get_fcst[0]['temperature'][0].ts)

    @property
    def wx_epsg_bbox(self):
        """A slice of test-data located in shyft-data repository/wrf."""
        EPSG = 32633
        x0 = 27000.0  # lower left
        y0 = 6780000.0  # lower right
        nx = 31
        ny = 21
        dx = 1000.0
        dy = 1000.0
        #return EPSG, ([x0, x0 + nx * dx, x0 + nx * dx, x0], [y0, y0, y0 + ny * dy, y0 + ny * dy])
        return EPSG, ([x0, x0 + nx * dx, x0 + nx * dx, x0], [y0, y0, y0 + ny * dy, y0 + ny * dy]), box(x0, y0, x0 + dx * nx, y0 + dy * ny)




if __name__ == "__main__":
    unittest.main()
