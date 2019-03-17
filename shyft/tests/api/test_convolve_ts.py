import unittest
import numpy as np

from shyft.api import Calendar
from shyft.api import DoubleVector
from shyft.api import TimeAxisFixedDeltaT
from shyft.api import TimeSeries
from shyft.api import convolve_policy
from shyft.api import deltahours
from shyft.api import point_interpretation_policy as point_fx


class ConvolveTs(unittest.TestCase):
    """Verify and illustrate the ts.convolve_w(weights,policy)

     """

    def test_convolve_policy(self):
        utc = Calendar()
        ts = TimeSeries(ta=TimeAxisFixedDeltaT(utc.time(2001, 1, 1), deltahours(1), 24), fill_value=10.0, point_fx=point_fx.POINT_AVERAGE_VALUE)
        w = DoubleVector.from_numpy([0.05, 0.15, 0.6, 0.15, 0.05])
        alignments = [convolve_policy.BACKWARD, convolve_policy.CENTER, convolve_policy.FORWARD]
        for alignment in alignments:
            # policy USE_NEAREST to ensure mass-balance between source and cts
            cts = ts.convolve_w(w, convolve_policy.USE_NEAREST | alignment)
            self.assertIsNotNone(cts)
            self.assertEquals(len(cts), len(ts))
            self.assertEquals(cts.values.to_numpy().sum(), ts.values.to_numpy().sum())

    def test_bind_info(self):
        ts = TimeSeries("a")
        w = DoubleVector.from_numpy([0.05, 0.15, 0.6, 0.15, 0.05])
        tsd = ts.convolve_w(w, convolve_policy.USE_NEAREST | convolve_policy.BACKWARD)
        bi = tsd.find_ts_bind_info()
        self.assertEqual(len(bi), 1)

    def test_invalid_convolution(self):
        utc = Calendar()
        ts = TimeSeries(ta=TimeAxisFixedDeltaT(utc.time(2001, 1, 1), deltahours(1), 12),
                        fill_value=10.0,
                        point_fx=point_fx.POINT_AVERAGE_VALUE)
        even_sized_kernel = DoubleVector.from_numpy([0.5, 0.5])
        oversized_kernel = DoubleVector.from_numpy(np.ones(13))
        with self.assertRaises(RuntimeError):
            ts.convolve_w(even_sized_kernel, convolve_policy.CENTER)
        with self.assertRaises(RuntimeError):
            ts.convolve_w(oversized_kernel, convolve_policy.CENTER)
