import unittest
from numpy.testing import assert_array_almost_equal
from shyft.api import TimeSeries, Calendar, TimeAxis, deltahours, DoubleVector, POINT_AVERAGE_VALUE, POINT_INSTANT_VALUE


class BucketTimeSeries(unittest.TestCase):
    """
    Test/Verify the TimeSeries.bucket_to_hourly function ts.

    This function/time-series is specialized to deal with precipitation measurements
    using accumulated precipitation bucket signals.

    Due to defects in sensors, the signals fluctuate with temperature during the day.

    """

    def test_simple_case(self):
        utc = Calendar()
        t0 = utc.time(2018, 1, 1)
        dt = deltahours(1)
        dv = DoubleVector([
            0, 0, 0, 0, 0, 0,
            4, 4, 4, 4, 4, 4,
            1, 1, 1, 1, 1, 1,
            9, 9, 9, 9, 9, 9,
            9, 9, 9, 9, 9, 9,
            -6, -6, -6, -6, -6, -6,
            6, 6, 6, 6, 6, 6,
            12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12, 12
        ])
        ev = DoubleVector([
            0, 0, 0, 0, 0, 3,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 6,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 2,
            0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0
        ])
        acc_ts = TimeSeries(TimeAxis(t0, dt, len(dv)), dv, POINT_INSTANT_VALUE)
        h_ts = acc_ts.bucket_to_hourly(start_hour_utc=0, bucket_emptying_limit=-100)
        e_ts = TimeSeries(TimeAxis(t0, dt, len(ev)), ev, POINT_AVERAGE_VALUE)
        self.assertEquals(h_ts.time_axis, e_ts.time_axis)
        self.assertEquals(h_ts.point_interpretation(), e_ts.point_interpretation())
        assert_array_almost_equal(h_ts.values.to_numpy(), e_ts.values.to_numpy())
