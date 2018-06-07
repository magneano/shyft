import numpy as np
from numpy.testing import assert_array_almost_equal
import unittest
from shyft.api import TimeSeries, Calendar, TimeAxis, deltahours, DoubleVector, POINT_AVERAGE_VALUE


class InsideTimeSeries(unittest.TestCase):

    def test_simple_case(self):
        utc = Calendar()
        t0 = utc.time(2018, 1, 1)
        dt = deltahours(1)
        dv = DoubleVector()
        dv[:] = [1.0, 2.0, 2.5, 1.9, 3.0, 3.1, -1.0]
        i1_ex = [0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        ts = TimeSeries(TimeAxis(t0, dt, len(dv)), dv, POINT_AVERAGE_VALUE)
        i1 = ts.inside(2.0, 3.0)
        assert_array_almost_equal(i1.values.to_numpy(), np.array(i1_ex))

    def test_inverted_values(self):
        utc = Calendar()
        t0 = utc.time(2018, 1, 1)
        dt = deltahours(1)
        dv = DoubleVector()
        dv[:] = [1.0, 2.0, 2.5, 1.9, 3.0, 3.1, float('nan')]  # also verify nan-handling
        i1_ex = [1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
        ts = TimeSeries(TimeAxis(t0, dt, len(dv)), dv, POINT_AVERAGE_VALUE)
        i2 = ts.inside(min_v=2.0, max_v=3.0, nan_v=1.0, inside_v=0.0, outside_v=1.0)
        assert_array_almost_equal(i2.values.to_numpy(), np.array(i1_ex))
