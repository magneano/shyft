import unittest
from numpy.testing import assert_array_almost_equal

from shyft.api import TimeSeries, TimeAxis, time, POINT_AVERAGE_VALUE as stair_case, TsVector, DoubleVector, pow


class VerifyTimeSeriesPow(unittest.TestCase):

    def test_pow(self):
        """ verify pow(ts,num) pow(num,ts) pow(ts,ts) """
        a = TimeSeries(TimeAxis(time('2018-01-01T00:00:00Z'), time(3600), 3), DoubleVector([1.0, 2.0, 3.0]), stair_case)
        b = TimeSeries(TimeAxis(time('2018-01-01T00:00:00Z'), time(3600), 3), DoubleVector([2.0, 2.0, 2.0]), stair_case)
        assert_array_almost_equal([1, 4, 9], a.pow(2.0).values.to_numpy())
        assert_array_almost_equal([1, 4, 9], pow(a, 2.0).values.to_numpy())
        assert_array_almost_equal([2, 4, 8], pow(2.0, a).values.to_numpy())
        assert_array_almost_equal([2, 4, 8], pow(b, a).values.to_numpy())
        assert_array_almost_equal([2, 4, 8], b.pow(a).values.to_numpy())

    def test_pow_vector(self):
        """ verify tsvector pow(tsv,num) pow(num,tsv) pow(tsv,ts) pow(tsv,tsv) """
        a = TsVector([
            TimeSeries(TimeAxis(time('2018-01-01T00:00:00Z'), time(3600), 3), DoubleVector([1.0, 2.0, 3.0]), stair_case)
        ]
        )
        b = TimeSeries(TimeAxis(time('2018-01-01T00:00:00Z'), time(3600), 3), DoubleVector([2.0, 2.0, 2.0]), stair_case)

        assert_array_almost_equal([1, 4, 9], a.pow(2.0)[0].values.to_numpy())
        assert_array_almost_equal([1, 4, 9], pow(a, 2.0)[0].values.to_numpy())
        assert_array_almost_equal([2, 4, 8], pow(2.0, a)[0].values.to_numpy())
        assert_array_almost_equal([2, 4, 8], pow(b, a)[0].values.to_numpy())
        assert_array_almost_equal([1, 4, 9], a.pow(b)[0].values.to_numpy())
        # assert_array_almost_equal([1, 4, 9], pow(a,b)[0].values.to_numpy())  # not supported
