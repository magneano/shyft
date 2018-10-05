import numpy as np
from numpy.testing import assert_array_almost_equal
import unittest
from shyft.api import TimeSeries, Calendar, TimeAxis, deltahours, DoubleVector, POINT_AVERAGE_VALUE, derivative_method


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

    def test_inside_of_derivative(self):
        """Created in response to https://github.com/statkraft/shyft/issues/352"""
        values = [1, 1, 1, 1]
        utc = Calendar()
        data = np.array(values, dtype='float64')
        data_ta = TimeAxis(utc.time(2015, 1, 1), 3600, len(data))

        orig = TimeSeries(data_ta, data, point_interpretation_policy.POINT_AVERAGE_VALUE)
        orig_derivative_inside_inf = orig.derivative(derivative_method.BACKWARD).inside(-float('inf'), float('inf'), 0,
                                                                                        1, 0)

        def check_nan(ts):
            """Method that returns 1 for all timesteps that contain nan."""
            return ts.inside(-float('inf'), float('inf'), 1, 0, 0).values.to_numpy()

        np.testing.assert_equal(check_nan(orig + orig_derivative_inside_inf), [0., 0., 0., 0.],
                                'TimeSeries.inside() should not find any NaN-values in this TimeSeries')

