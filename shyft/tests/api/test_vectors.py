from builtins import range

from shyft import api
import numpy as np
from numpy.testing import assert_array_almost_equal
import unittest


class Vectors(unittest.TestCase):
    """
    Some basic test to ensure that the exposure of c++ vector<T> is working as expected
    """

    def test_double_vector(self):
        dv_from_list = api.DoubleVector([x for x in range(10)])
        dv_np = np.arange(10.0)
        dv_from_np = api.DoubleVector.from_numpy(dv_np)
        self.assertEqual(len(dv_from_list), 10)
        assert_array_almost_equal(dv_from_list.to_numpy(), dv_np)
        assert_array_almost_equal(dv_from_np.to_numpy(), dv_np)
        dv_from_np[5] = 8
        dv_from_np.append(11)
        dv_from_np.push_back(12)
        dv_np[5] = 8
        dv_np.resize(12)
        dv_np[10] = 11
        dv_np[11] = 12
        assert_array_almost_equal(dv_from_np.to_numpy(), dv_np)
        # this does not work yet
        # nv= api.DoubleVector(dv_np).. would be very nice!

    def test_int_vector(self):
        dv_from_list = api.IntVector([x for x in range(10)])
        dv_np = np.arange(10, dtype=np.int32)  # notice, default is int64, which does not convert automatically to int32
        dv_from_np = api.IntVector.from_numpy(dv_np)
        self.assertEqual(len(dv_from_list), 10)
        assert_array_almost_equal(dv_from_list.to_numpy(), dv_np)
        assert_array_almost_equal(dv_from_np.to_numpy(), dv_np)
        dv_from_np[5] = 8
        dv_from_np.append(11)
        dv_from_np.push_back(12)
        dv_np[5] = 8
        dv_np.resize(12)
        dv_np[10] = 11
        dv_np[11] = 12
        assert_array_almost_equal(dv_from_np.to_numpy(), dv_np)

    def test_int64_vector(self):
        dv_from_list = api.Int64Vector([x for x in range(10)])
        dv_np = np.arange(10, dtype=np.int32)  # notice, default is int64, which does not convert automatically to int32
        dv_from_np = api.Int64Vector.from_numpy(dv_np)
        self.assertEqual(len(dv_from_list), 10)
        assert_array_almost_equal(dv_from_list.to_numpy(), dv_np)
        assert_array_almost_equal(dv_from_np.to_numpy(), dv_np)
        dv_from_np[5] = 8
        dv_from_np.append(11)
        dv_from_np.append(12)
        dv_np[5] = 8
        dv_np.resize(12)
        dv_np[10] = 11
        dv_np[11] = 12
        assert_array_almost_equal(dv_from_np.to_numpy(), dv_np)

    def test_utctime_vector(self):
        dv_from_list = api.UtcTimeVector([x for x in range(10)])
        dv_np = np.arange(10, dtype=np.int64)
        dv_from_np = api.UtcTimeVector.from_numpy(dv_np)
        self.assertEqual(len(dv_from_list), 10)
        assert_array_almost_equal(dv_from_list.to_numpy(), dv_np)
        assert_array_almost_equal(dv_from_np.to_numpy(), dv_np)
        dv_from_np[5] = api.time(8)  # it should also have accepted any number here
        dv_from_np[5] = 8.5  # 8.5 seconds
        self.assertAlmostEqual(dv_from_np[5].seconds, 8.5)  # verify it stores microseconds
        dv_from_np.append(api.time(11))
        dv_from_np.push_back(12)  # this one takes any that could go as seconds
        dv_from_np.push_back(api.time(13))  # this one takes any that could go as seconds
        dv_np[5] = 8.5  # python/numpy silently ignore float -> int64
        dv_np.resize(13)
        dv_np[10] = 11
        dv_np[11] = 12
        dv_np[12] = 13
        assert_array_almost_equal(dv_from_np.to_numpy(), dv_np)
        dv2 = dv_from_np.to_numpy_double()
        self.assertAlmostEqual(dv2[5], 8.5)  # verify that to_numpy_double preserves microsecond

    def test_string_vector(self):
        # NOTE: support for string vector is very limited, e.g. numpy does not work, only lists
        #     but for now this is sufficient in Shyft
        s_list = ['abc', 'def']
        dv_from_list = api.StringVector(s_list)
        self.assertEqual(len(dv_from_list), 2)
        for i in range(len(dv_from_list)):
            self.assertEqual(s_list[i], dv_from_list[i])

    def test_geopoint_vector(self):
        gpv = api.GeoPointVector()
        for i in range(5):
            gpv.append(api.GeoPoint(i, i*2, i*3))
        self.assertEqual(len(gpv), 5)
        for i in range(5):
            self.assertEqual(gpv[i].x, i)
            self.assertEqual(gpv[i].y, i*2)
            self.assertEqual(gpv[i].z, i*3)

    def test_geo_cell_data_vector(self):
        gcdv = api.GeoCellDataVector()
        for i in range(5):
            p = api.GeoPoint(100, 200, 300)
            ltf = api.LandTypeFractions()
            ltf.set_fractions(glacier=0.1, lake=0.1, reservoir=0.1, forest=0.1)
            gcd = api.GeoCellData(p, 1000000.0, i, 0.9, ltf)
            gcdv.append(gcd)
        self.assertEqual(len(gcdv), 5)
        for i in range(5):
            self.assertEqual(gcdv[i].catchment_id(), i)
        g2 = api.GeoCellDataVector(gcdv)  # copy construct a new
        self.assertTrue(g2 == gcdv)
        g2[0].set_catchment_id(10)
        self.assertTrue(g2 != gcdv)
        # serialize
        gcdv_s= gcdv.serialize()
        self.assertGreater(len(gcdv_s),2)
        gcdv_deserialized = api.GeoCellDataVector.deserialize(gcdv_s)
        self.assertIsNotNone(gcdv_deserialized)
        self.assertTrue(gcdv_deserialized == gcdv)

