import unittest

import numpy as np
from numpy.testing import assert_array_almost_equal

from shyft.api import TimeSeries, Calendar, TimeAxis, deltahours, DoubleVector, POINT_AVERAGE_VALUE


class DecodeTimeSeries(unittest.TestCase):
    """
    Test/Verify the TimeSeries bit decode function.
    Typical use of TimeSeries.decode(start_bit=n, n_bits=4) is for remote-sensing quality/sensing status bit
    handling.
    The .decode function help you extract the interesting part of bit-encoded status to time-series that can
    be dealt with using math & statistics
    """

    def test_simple_case(self):
        utc = Calendar()
        t0 = utc.time(2018, 1, 1)
        dt = deltahours(1)
        dv = DoubleVector()
        dv[:] = [1.2, 0.0, 2.0, 5.0, 15.0, float('nan'), -1.0]  # these are bit-encoded values, note 1.2 -> 1.0
        i0_1_e = [1.0, 0.0, 0.0, 1.0, 1.0, float('nan'), float('nan')]  # expected values for bit 0 1-bit
        i1_3_e = [0.0, 0.0, 1.0, 2.0, 7.0, float('nan'), float('nan')]  # expected values for bit 3 3-bits
        ts = TimeSeries(TimeAxis(t0, dt, len(dv)), dv, POINT_AVERAGE_VALUE)
        i0_1 = ts.decode(start_bit=0, n_bits=1)
        i1_3 = ts.decode(start_bit=1, n_bits=3)
        assert_array_almost_equal(i0_1.values.to_numpy(), np.array(i0_1_e))
        assert_array_almost_equal(i1_3.values.to_numpy(), np.array(i1_3_e))

    def test_error_handling(self):
        utc = Calendar()
        t0 = utc.time(2018, 1, 1)
        dt = deltahours(1)
        dv = DoubleVector()
        dv[:] = [1.0, 2.0, 2.5, 1.9, 3.0, 3.1, float('nan')]  # also verify nan-handling
        ts = TimeSeries(TimeAxis(t0, dt, len(dv)), dv, POINT_AVERAGE_VALUE)
        try:
            ts.decode(start_bit=0, n_bits=0)
            self.assertTrue(False, 'This should throw, n_bits >0')
        except RuntimeError as re:
            pass

        try:
            ts.decode(start_bit=41, n_bits=12)
            self.assertTrue(False, 'This should throw, start_bit + n_bits >52')
        except RuntimeError as re:
            pass

        try:
            ts.decode(start_bit=-1, n_bits=12)
            self.assertTrue(False, 'This should throw, start_bit >=0')
        except RuntimeError as re:
            pass
