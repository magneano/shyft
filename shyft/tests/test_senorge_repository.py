import unittest
from os import path
from netCDF4 import Dataset
from pyproj import Proj
from pyproj import transform
import numpy as np

from shyft import shyftdata_dir
from shyft import api
from shyft.repository.netcdf.senorge_data_repository import SeNorgeDataRepository
from shyft.repository.netcdf.senorge_data_repository import SeNorgeDataRepositoryError

from shapely.geometry import box


class SeNorgeDataRepositoryTestCase(unittest.TestCase):
    def test_get_timeseries(self):
        """
        #Simple regression test of senorge data repository.
        """
        EPSG, bbox, bpoly = self.senorge_epsg_bbox

        # Period start
        n_days = 5
        t0 = api.YMDhms(1957, 1, 2)
        utc = api.Calendar()  # No offset gives Utc
        step = api.deltahours(24)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_days*24))

        base_dir = path.join(shyftdata_dir, "repository", "senorge_data_repository")
        f1 = "senorge_test.nc"

        senorge1 = SeNorgeDataRepository(EPSG, base_dir, filename=f1, allow_subset=True)
        senorge1_data_names = ("temperature",)
        sources = senorge1.get_timeseries(senorge1_data_names, period, geo_location_criteria=bpoly)
        self.assertTrue(len(sources) > 0)

        self.assertTrue(set(sources) == set(senorge1_data_names))
        p0 = sources["temperature"][0].ts
        self.assertTrue(p0.size() == n_days + 1)
        self.assertTrue(p0.time(0), period.start)

        # check boundaries of period
        t0 = utc.time(1956, 12, 31, 6)  # referring to begin of first utc period in test senorge-dataset
        t1 = utc.time(1957, 1, 10, 6)  # referring to end of last utc period in test senorge-dataset
        period = api.UtcPeriod(t0, t1)
        sources = senorge1.get_timeseries(senorge1_data_names, period, geo_location_criteria=bpoly)
        self.assertTrue(sources['temperature'][0].ts.total_period() == period)

        t0 = utc.time(1956, 12, 31, 6, 1) # should expand series to earlier start time
        t1 = utc.time(1957, 1, 10, 5, 59) # should expand series to later end time
        period = api.UtcPeriod(t0, t1)
        sources = senorge1.get_timeseries(senorge1_data_names, period, geo_location_criteria=bpoly)
        self.assertTrue(sources['temperature'][0].ts.total_period().start < period.start)
        self.assertTrue(sources['temperature'][0].ts.total_period().end > period.end)

        t0 = utc.time(1956, 12, 31, 5, 59)  # out of range
        period = api.UtcPeriod(t0, t1)
        with self.assertRaises(SeNorgeDataRepositoryError) as context:
            senorge1.get_timeseries(senorge1_data_names, period, geo_location_criteria=bpoly)
        self.assertEqual("The earliest time in repository (1956-12-31T06:00:00Z) is later than the start of the period for which data is requested (1956-12-31T05:59:00Z)", context.exception.args[0])

        t0 = utc.time(1956, 12, 31, 6, 1)
        t1 = utc.time(1957, 1, 10, 6, 1) # out of range
        period = api.UtcPeriod(t0, t1)
        with self.assertRaises(SeNorgeDataRepositoryError) as context:
            senorge1.get_timeseries(senorge1_data_names, period, geo_location_criteria=bpoly)
        self.assertEqual("The latest time in repository (1957-01-10T06:00:00Z) is earlier than the end of the period for which data is requested (1957-01-10T06:01:00Z)", context.exception.args[0])


    def test_get_dummy(self):
        """
        #Simple regression test of dummy variables from `utils.dummy_var`
        """
        EPSG, bbox, bpoly = self.senorge_epsg_bbox

        # Period start
        n_days = 5
        t0 = api.YMDhms(1957, 1)
        utc = api.Calendar()  # No offset gives Utc
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_days * 24))

        base_dir = path.join(shyftdata_dir, "repository", "senorge_data_repository")
        f1 = "senorge_test.nc"

        senorge1 = SeNorgeDataRepository(EPSG, base_dir, filename=f1, allow_subset=True)
        senorge1_data_names = ("radiation","temperature","wind_speed")
        sources = senorge1.get_timeseries(senorge1_data_names, period, geo_location_criteria=bpoly)
        self.assertTrue(set(sources) == set(senorge1_data_names))
        r0 = sources["radiation"][0].ts
        p0 = sources["temperature"][0].ts
        self.assertTrue(p0.total_period().start <= period.start and p0.total_period().end >= period.end)
        self.assertTrue(r0.time(0), period.start)

        t0 = api.YMDhms(1957, 1)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_days * 23))
        sources = senorge1.get_timeseries(senorge1_data_names, period, geo_location_criteria=bpoly)
        r0 = sources["radiation"][0].ts
        p0 = sources["temperature"][0].ts
        self.assertTrue(p0.total_period().start <= period.start and p0.total_period().end >= period.end)
        self.assertTrue(r0.time(0) == period.start)
        self.assertTrue(r0.time_axis.time(1) - r0.time_axis.time(0) == 86400)





    @property
    def senorge_epsg_bbox(self):
        """A slice of test-data located in shyft-data repository/senorge."""
        EPSG = 32633
        x0 = 374000.0  # lower left
        y0 = 6450500.0  # lower right
        nx = 102.0
        ny = 121.0
        dx = 500.0
        dy = 500.0
        #return EPSG, ([x0, x0 + nx * dx, x0 + nx * dx, x0], [y0, y0, y0 + ny * dy, y0 + ny * dy])
        return EPSG, ([x0, x0 + nx * dx, x0 + nx * dx, x0], [y0, y0, y0 + ny * dy, y0 + ny * dy]), box(x0, y0, x0 + dx * nx, y0 + dy * ny)

    def test_wrong_directory(self):
        with self.assertRaises(SeNorgeDataRepositoryError) as context:
            SeNorgeDataRepository(32632, "Foobar", filename="")
        self.assertEqual("No such directory 'Foobar'", context.exception.args[0])

    def test_wrong_file(self):
        with self.assertRaises(SeNorgeDataRepositoryError) as context:
            utc = api.Calendar()  # No offset gives Utc
            t0 = api.YMDhms(2015, 12, 25, 18)
            period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(30))
            ar1 = SeNorgeDataRepository(32632, shyftdata_dir, filename="plain_wrong.nc")
            ar1.get_timeseries(("temperature",), period, geo_location_criteria=None)
        self.assertTrue(all(x in context.exception.args[0] for x in ["File", "not found"]))

    def test_wrong_elevation_file(self):
        with self.assertRaises(SeNorgeDataRepositoryError) as context:
            SeNorgeDataRepository(32632, shyftdata_dir, filename="", elevation_file="plain_wrong.nc")
        self.assertTrue(all(x in context.exception.args[0] for x in ["Elevation file",
                                                                     "not found"]))

    # TODO: geo_location_criteria is now of type shapely.geometry.Polygon. Replace with test verifying that error is raised when no point is found within polygon.
    # def test_non_overlapping_bbox(self):
    #     EPSG, bbox, bpoly = self.wrf_epsg_bbox
    #     bbox = list(bbox)
    #     bbox[0] = [-100000.0, -90000.0, -90000.0, -100000]
    #     bpoly = box(min(bbox[0]), min(bbox[1]), max(bbox[0]), max(bbox[1]))
    #     # Period start
    #
    #     year = 1999
    #     month = 10
    #     n_hours = 30
    #     date_str = "{}-{:02}".format(year, month)
    #     utc = api.Calendar()  # No offset gives Utc
    #     t0 = api.YMDhms(year, month)
    #     period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))
    #
    #     base_dir = path.join(shyftdata_dir, "repository", "wrf_data_repository")
    #     filename = "wrfout_d03_{}".format(date_str)
    #     reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename)
    #     data_names = ("temperature", "wind_speed", "precipitation", "relative_humidity")
    #     with self.assertRaises(SeNorgeDataRepositoryError) as context:
    #         reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)
    #     self.assertEqual("Bounding box doesn't intersect with dataset.",
    #                      context.exception.args[0])

    # TODO: geo_location_criteria=None now returns all points in dataset. Replace with test verifying that all points are returned when geo_location_criteria=None.
    # def test_missing_bbox(self):
    #     EPSG, _, _ = self.wrf_epsg_bbox
    #     # Period start
    #     year = 1999
    #     month = 10
    #     n_hours = 30
    #     date_str = "{}-{:02}".format(year, month)
    #     utc = api.Calendar()  # No offset gives Utc
    #     t0 = api.YMDhms(year, month)
    #     period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))
    #
    #     base_dir = path.join(shyftdata_dir, "repository", "wrf_data_repository")
    #     filename = "wrfout_d03_{}".format(date_str)
    #     reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename)
    #     data_names = ("temperature", "wind_speed", "precipitation", "relative_humidity")
    #     with self.assertRaises(SeNorgeDataRepositoryError) as context:
    #         reader.get_timeseries(data_names, period, geo_location_criteria=None)
    #     self.assertEqual("A bounding box must be provided.", context.exception.args[0])

    def test_tiny_bbox(self):
        EPSG, _, _ = self.senorge_epsg_bbox
        x0 = 374499.5  # lower left
        y0 = 6451499.5  # lower right
        nx = 1.0
        ny = 1.0
        dx = 1.0
        dy = 1.0
        bbox = ([x0, x0 + nx * dx, x0 + nx * dx, x0], [y0, y0, y0 + ny * dy, y0 + ny * dy])
        bpoly = box(min(bbox[0]), min(bbox[1]), max(bbox[0]), max(bbox[1]))

        # Period start
        year = 1957
        month = 1
        n_hours = 30
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "senorge_data_repository")
        filename = "senorge_test.nc"
        reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename,
                                   padding=0)
        data_names = ("temperature",)

        tss = reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)

        for name, ts in tss.items():
            self.assertTrue(len(ts) == 1)


    def test_subsets(self):
        EPSG, bbox, bpoly = self.senorge_epsg_bbox
        # Period start
        year = 1957
        month = 1
        n_hours = 30
        date_str = "{}-{:02}".format(year, month)
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "senorge_data_repository")
        filename = "senorge_test.nc"

        data_names = ("temperature", "foo")
        allow_subset = False
        reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename,
                                   allow_subset=allow_subset)
        with self.assertRaises(SeNorgeDataRepositoryError) as context:
            reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)
        self.assertEqual("Input source types ['foo'] not supported", context.exception.args[0])
        allow_subset = True
        reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename,
                                   allow_subset=allow_subset)
        sources = reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)
        self.assertEqual(len(sources), len(data_names)-1)




if __name__ == "__main__":
    unittest.main()
