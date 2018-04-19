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
        Simple regression test of WRF data repository.
        """
        EPSG, bbox, bpoly = self.senorge_epsg_bbox

        # Period start
        n_hours = 60
        t0 = api.YMDhms(2015)
        date_str = "{}-{:02}".format(t0.year, t0.month)
        utc = api.Calendar()  # No offset gives Utc
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "senorge_data_repository", "senorge2")
        f1 = "seNorge2_PREC1d_grid_2015.nc"
#        f1 = "wrfout_d03_{}".format(date_str)

        senorge1 = SeNorgeDataRepository(EPSG, base_dir, filename=f1, allow_subset=True)
        senorge1_data_names = ("precipitation")
        sources = senorge1.get_timeseries(senorge1_data_names, period, geo_location_criteria=bpoly)
        self.assertTrue(len(sources) > 0)

        self.assertTrue(set(sources) == set(senorge1_data_names))
        self.assertTrue(sources["temperature"][0].ts.size() == n_hours + 1)
        p0 = sources["precipitation"][0].ts
#        temp0 = sources["temperature"][0].ts
        self.assertTrue(p0.size() == n_hours + 1)
#        self.assertTrue(p0.time(0) == temp0.time(0))
        self.assertTrue(p0.time(0), period.start)

        # Number test:
        # asserting shyft-sources time series are same as time series of corresponding location in wrf dataset.
        dset = Dataset(path.join(base_dir, f1))
        lat = dset.variables["Y"]
        lon = dset.variables["X"]

        senorge_data = {}

        senorge_data["temperature"] = dset.variables["mean_temperature"][:]
        senorge_data["precipitation"] = dset.variables["precipitaion_amount"][:]
        senorge_data["temperature"] -= 273.16

        data_cs = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
        target_cs = "+init=EPSG:32643"
        data_proj = Proj(proj=data_cs)
        target_proj = Proj(target_cs)
        x, y = transform(data_proj, target_proj, lon[0, :, :], lat[0, :, :])

        for name, senorge_d in senorge_data.items():
            srs = sources[name]
            for i, s in enumerate(srs):
                mp = s.mid_point()
                x_ts, y_ts, z_ts = mp.x, mp.y, mp.z
                ts = s.ts
                ts_values = ts.v.to_numpy()

                # find indixes in senorge-dataset
                m = (x == x_ts) & (y == y_ts)
                idxs = np.where(m > 0)
                x_idx, y_idx = idxs[0][0], idxs[1][0]  # assumung geo-location is unique in dataset
                self.assertTrue(all(ts_values == senorge_d[:n_hours + 1, x_idx, y_idx]),
                                "senorge and shyft-TS of {} are not the same.".format(name))
                # if i ==0:
                #    plt.figure()
                #    plt.plot(ts_values)
                #    plt.title([name])
                #    plt.show()

    @property
    def senorge_epsg_bbox(self):
        """A slice of test-data located in shyft-data repository/senorge."""
        EPSG = 32643
        x0 = 674085.0  # lower left
        y0 = 3476204.0  # lower right
        nx = 102
        ny = 121
        dx = 1000.0
        dy = 1000.0
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
        x0 = 726270.0  # lower left
        y0 = 3525350.0  # lower right
        nx = 1
        ny = 1
        dx = 1.0
        dy = 1.0
        bbox = ([x0, x0 + nx * dx, x0 + nx * dx, x0], [y0, y0, y0 + ny * dy, y0 + ny * dy])
        bpoly = box(min(bbox[0]), min(bbox[1]), max(bbox[0]), max(bbox[1]))

        # Period start
        year = 2015
        month = 1
        n_hours = 30
        date_str = "{}-{:02}".format(year, month)
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "senorge_data_repository", "senorge2")
        filename = "seNorge2_PREC1d_grid_2015.nc"
        reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename,
                                   padding=0)
        data_names = ("precipitation")

        tss = reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)

        for name, ts in tss.items():
            self.assertTrue(len(ts) == 1)

    def test_subsets(self):
        EPSG, bbox, bpoly = self.senorge_epsg_bbox
        # Period start
        year = 2015
        month = 1
        n_hours = 30
        date_str = "{}-{:02}".format(year, month)
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "senorge_data_repository", "senorge2")
        filename = "seNorge2_PREC1d_grid_2015.nc"


        data_names = ("precipitation")
        allow_subset = False
        reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename,
                                   allow_subset=allow_subset)
        with self.assertRaises(SeNorgeDataRepositoryError) as context:
            reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)
        self.assertEqual("Could not find all data fields", context.exception.args[0])
        allow_subset = True
        reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename,
                                   allow_subset=allow_subset)
        try:
            sources = reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)
        except SeNorgeDataRepositoryError as e:
            self.fail("AromeDataRepository.get_timeseries(data_names, period, None) "
                      "raised AromeDataRepositoryError unexpectedly.")
        self.assertEqual(len(sources), len(data_names) - 1)

    """
    def test_rel_hum_only(self):

        print("rel hum test: ")
        # relative humidity needs temperature and pressure to be calculated
        EPSG, bbox, bpoly = self.wrf_epsg_bbox
        # Period start
        year = 1999
        month = 10
        n_hours = 30
        date_str = "{}-{:02}".format(year, month)
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "wrf_data_repository")
        filename = "wrfout_d03_{}".format(date_str)

        data_names = ["relative_humidity"]
        reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename)
        sources = reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)

        self.assertTrue(list(sources.keys()) == ["relative_humidity"])

        # allow_subset = True
        # reader = SeNorgeDataRepository(EPSG, base_dir, filename=filename,
        #                             bounding_box=bbox, allow_subset=allow_subset)
        # try:
        #    sources = reader.get_timeseries(data_names, period, None)
        # except SeNorgeDataRepositoryError as e:
        #    self.fail("AromeDataRepository.get_timeseries(data_names, period, None) "
        #              "raised AromeDataRepositoryError unexpectedly.")
        # self.assertEqual(len(sources), len(data_names) - 1)
    """

if __name__ == "__main__":
    unittest.main()
