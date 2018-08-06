from shyft.api import TimeSeries,TimeAxis,Calendar,UtcPeriod,TsVector,DoubleVector,POINT_AVERAGE_VALUE,POINT_INSTANT_VALUE

#import numpy as np
import math
#from numpy.testing import assert_array_almost_equal
import unittest


class TimeSeriesDerivative(unittest.TestCase):

    def test_linear(self):
        ta = TimeAxis(0,10,6)
        tsv=DoubleVector([1, 1, 2, 3, -1.0, 5.0])
        f = TimeSeries(ta, tsv, POINT_INSTANT_VALUE)
        d_f = f.derivative()
        self.assertEqual(len(f),len(d_f))
        self.assertAlmostEqual(d_f.value(0), 0.0)
        self.assertAlmostEqual(d_f.value(1), 0.1)
        self.assertAlmostEqual(d_f.value(2), 0.1)
        self.assertAlmostEqual(d_f.value(3), -0.4)
        self.assertAlmostEqual(d_f(f.time(3) + 5), -0.4)
        self.assertAlmostEqual(d_f.value(4), 0.6)
        self.assertFalse(math.isfinite(d_f.value(5)))
        self.assertFalse(math.isfinite(d_f(f.time(5))))
        self.assertTrue(math.isfinite(d_f(f.time(5) - 1)))
        v = d_f.values
        self.assertAlmostEqual(len(v), len(f))
        self.assertAlmostEqual(v[0], 0.0)
        self.assertAlmostEqual(v[1], 0.1)
        self.assertAlmostEqual(v[2], 0.1)
        self.assertAlmostEqual(v[3], -0.4)
        self.assertAlmostEqual(v[4], 0.6)
        self.assertFalse(math.isfinite(v[5]))

    def test_linear_vector(self):
        ta = TimeAxis(0,10,6)
        tsv=DoubleVector([1, 1, 2, 3, -1.0, 5.0])
        fv = TsVector([TimeSeries(ta, tsv, POINT_INSTANT_VALUE)])
        d_fv = fv.derivative()
        f=fv[0]
        d_f=d_fv[0]
        self.assertEqual(len(f),len(d_f))
        self.assertAlmostEqual(d_f.value(0), 0.0)
        self.assertAlmostEqual(d_f.value(1), 0.1)
        self.assertAlmostEqual(d_f.value(2), 0.1)
        self.assertAlmostEqual(d_f.value(3), -0.4)
        self.assertAlmostEqual(d_f(f.time(3) + 5), -0.4)
        self.assertAlmostEqual(d_f.value(4), 0.6)
        self.assertFalse(math.isfinite(d_f.value(5)))
        self.assertFalse(math.isfinite(d_f(f.time(5))))
        self.assertTrue(math.isfinite(d_f(f.time(5) - 1)))
        v = d_f.values
        self.assertAlmostEqual(len(v), len(f))
        self.assertAlmostEqual(v[0], 0.0)
        self.assertAlmostEqual(v[1], 0.1)
        self.assertAlmostEqual(v[2], 0.1)
        self.assertAlmostEqual(v[3], -0.4)
        self.assertAlmostEqual(v[4], 0.6)
        self.assertFalse(math.isfinite(v[5]))
