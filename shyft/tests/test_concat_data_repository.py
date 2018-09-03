from __future__ import print_function
import unittest
from os import path

from shyft import shyftdata_dir
from shyft import api
from shyft.repository.netcdf.concat_data_repository import ConcatDataRepository
#from shyft.repository.netcdf.concant_data_repository import ConcatDataRepositoryError
from shapely.geometry import box
#import netCDF4
import numpy as np


class ConcatDataRepositoryTestCase(unittest.TestCase):
    """
    TODO: This test needs to be written, along with utility that generates concat from sources
    """

    @property
    def _epsg_bbox(self):
        """A slice of test-data located in shyft-data repository/arome."""
        EPSG = 32632
        x0 = 436100.0  # lower left
        y0 = 6823000.0  # lower right
        nx = 74
        ny = 24
        dx = 1000.0
        dy = 1000.0
        return EPSG, ([x0, x0 + nx*dx, x0 + nx*dx, x0], [y0, y0, y0 + ny*dy, y0 + ny*dy]), box(x0, y0, x0 + dx*nx, y0 + dy*ny)

    def test_transform_functions_fixed_interval(self):
        """
        test the _transform_raw function.
        """
        return
        # TODO: add concat file in shyft-data, then implement the tests
        EPSG, bbox, bpoly = self._epsg_bbox

        # Period start
        t0 = api.YMDhms(2015, 8, 24, 0)
        date_str = "{}{:02}{:02}_{:02}".format(t0.year, t0.month, t0.day, t0.hour)
        utc = api.Calendar()  # No offset gives Utc

        f1 = path.join(shyftdata_dir, "repository", "arome_data_repository",f"arome_metcoop_red_default2_5km_{date_str}_diff_time_unit.nc")
        ar1 = ConcatDataRepository(epsg=EPSG, filename=f1)
        np_raw_array = np.array(
                    [  # 0  # 1 #  2 #  3
                        [1.0, 2.0, 3.0, 4.0],
                        [1.1, 2.1, 3.1, 4.1],
                        [1.2, 2.2, 3.2, 4.2],
                        [1.4, 2.5, 3.6, 4.7]
                    ], dtype=np.float64
                )
        raw_values = {
            'wind_speed': (np_raw_array, 'wind_speed', 'm/s'),
            'rel_hum': (np_raw_array, 'relative_humidity_2m', '?'),
            'temperature': (273.15 + np_raw_array, 'air_temperature_2m', 'K'),
            'radiation': (3600.0*np_raw_array, 'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time', 'W s/m2'),
            'prepitation_acc': (np_raw_array, 'precipitation_amount_acc', 'Mg/m^2'),
            'prepitation': (np_raw_array, 'precipitation_amount', 'mm')
        }
        raw_time = np.array([0, 3600, 7200, 10800], dtype=np.int64)

        rd = ar1._transform_raw(raw_values, raw_time)
        ta3 = api.TimeAxis(api.time(0), api.time(3600), 3)
        ta4 = api.TimeAxis(api.time(0), api.time(3600), 4)
        e_precip_acc = np.array(
                    [  # 0  # 1 #  2 #  3
                        [100.0, 100.0, 100.0, 100.0],
                        [100.0, 100.0, 100.0, 100.0],
                        [200.0, 300.0, 400.0, 500.0],
                    ], dtype=np.float64
        )
        e_precip = np.array(
                    [  # 0  # 1 #  2 #  3
                        [1.1, 2.1, 3.1, 4.1],
                        [1.2, 2.2, 3.2, 4.2],
                        [1.4, 2.5, 3.6, 4.7]
                    ], dtype=np.float64
                )
        e_rad = np.array(
            [  # 0  # 1 #  2 #  3
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.1, 0.1, 0.1],
                [0.2, 0.3, 0.4, 0.5],
            ], dtype=np.float64
        )
        e = {
            'wind_speed': (np_raw_array, ta4),
            'rel_hum': (np_raw_array, ta4),
            'temperature': (np_raw_array, ta4),
            'radiation': (e_rad, ta3),
            'prepitation_acc': (e_precip_acc, ta3),
            'prepitation': (e_precip, ta3)
        }

        self.assertIsNotNone(rd)
        for k, r in rd.items():
            self.assertTrue(k in e)
            self.assertEqual(r[1], e[k][1], "expect correct time-axis")
            self.assertTrue(np.allclose(r[0], e[k][0]), "expect exact correct values")

    def test_transform_functions_variable_interval(self):
        """
        test the _transform_raw function.
        """
        return
        # TODO: add concat file in shyft-data, then implement the tests
        EPSG, bbox, bpoly = self._epsg_bbox

        # Period start
        n_hours = 30
        t0 = api.YMDhms(2015, 8, 24, 0)
        date_str = "{}{:02}{:02}_{:02}".format(t0.year, t0.month, t0.day, t0.hour)
        utc = api.Calendar()  # No offset gives Utc

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        f1 = "arome_metcoop_red_default2_5km_{}_diff_time_unit.nc".format(date_str)
        ar1 = ConcatDataRepository(EPSG, base_dir, filename=f1)
        np_raw_array = np.array(
                    [  # 0  # 1 #  2 #  3
                        [1.0, 2.0, 3.0, 4.0],
                        [1.1, 2.1, 3.1, 4.1],
                        [1.2, 2.2, 3.2, 4.2],
                        [1.4, 2.5, 3.6, 4.7]
                    ], dtype=np.float64
                )
        raw_values = {
            'wind_speed': (np_raw_array, 'wind_speed', 'm/s'),
            'rel_hum': (np_raw_array, 'relative_humidity_2m', '?'),
            'temperature': (273.15 + np_raw_array, 'air_temperature_2m', 'K'),
            'radiation': (3600.0*np_raw_array, 'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time', 'W s/m2'),
            'prepitation_acc': (np_raw_array, 'precipitation_amount_acc', 'Mg/m^2'),
            'prepitation': (np_raw_array, 'precipitation_amount', 'mm')
        }
        raw_time = np.array([0, 3600, 7200, 7200+2*3600], dtype=np.int64)  # last step is 2 hours!

        rd = ar1._transform_raw(raw_values, raw_time)
        ta3 = api.TimeAxis(api.UtcTimeVector(raw_time[:-1]), api.ConcatData(int(raw_time[-1])))
        ta4 = api.TimeAxis(api.UtcTimeVector(raw_time), api.ConcatData(int(raw_time[-1]+2*3600)))  # assume last step is also 2 hours
        e_precip_acc = np.array(
                    [  # 0  # 1 #  2 #  3
                        [100.0, 100.0, 100.0, 100.0],
                        [100.0, 100.0, 100.0, 100.0],
                        [100.0, 150.0, 200.0, 250.0],
                    ], dtype=np.float64
        )
        e_precip = np.array(
                    [  # 0  # 1 #  2 #  3
                        [1.1, 2.1, 3.1, 4.1],
                        [1.2, 2.2, 3.2, 4.2],
                        [1.4, 2.5, 3.6, 4.7]
                    ], dtype=np.float64
                )
        e_rad = np.array(
            [  # 0  # 1 #  2 #  3
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.15, 0.2, 0.25],
            ], dtype=np.float64
        )
        e = {
            'wind_speed': (np_raw_array, ta4),
            'rel_hum': (np_raw_array, ta4),
            'temperature': (np_raw_array, ta4),
            'radiation': (e_rad, ta3),
            'prepitation_acc': (e_precip_acc, ta3),
            'prepitation': (e_precip, ta3)
        }

        self.assertIsNotNone(rd)
        for k, r in rd.items():
            self.assertTrue(k in e)
            self.assertEqual(r[1], e[k][1], "expect correct time-axis")
            self.assertTrue(np.allclose(r[0], e[k][0]), "expect exact correct values")


if __name__ == "__main__":
    unittest.main()
