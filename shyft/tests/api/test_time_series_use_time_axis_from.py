import unittest
from shyft.api import TimeSeries, TimeAxis, Calendar, UtcPeriod, TsVector, DoubleVector, POINT_AVERAGE_VALUE, \
    POINT_INSTANT_VALUE
import numpy as np


class TimeSeriesUseTimeAxisFrom(unittest.TestCase):

    def test_basic_case(self):
        ta = TimeAxis(0, 10, 6)
        to = TimeAxis(10, 20, 3)
        va = DoubleVector([0, 10, 20, 30, 40.0, 50.0])
        a = TimeSeries(ta, va, POINT_INSTANT_VALUE)
        o = TimeSeries(to, fill_value=0.0, point_fx=POINT_INSTANT_VALUE)  # point-ip should not matter
        r = a.use_time_axis_from(o)
        self.assertIsNotNone(r)
        self.assertEquals(r.time_axis, o.time_axis)
        ev = DoubleVector([10.0, 30.0, 50.0]).to_numpy()
        self.assertTrue(np.allclose(r.values.to_numpy(), ev))
        self.assertEquals(r.point_interpretation(), a.point_interpretation())
        tsv = TsVector([a, 2.0 * a])
        rv = tsv.use_time_axis_from(o)
        for x in rv:
            self.assertEquals(x.time_axis, o.time_axis)
        self.assertTrue(np.allclose(rv[0].values.to_numpy(), ev))
        self.assertTrue(np.allclose(rv[1].values.to_numpy(), 2.0*ev))