import numpy as np
import math
import unittest
from shyft.api import Calendar
from shyft.api import TimeSeries
from shyft.api import TimeAxis
from shyft.api import TsVector
from shyft.api import TsVectorSet
from shyft.api import point_interpretation_policy as ts_point_fx
from shyft.api import deltahours
from shyft.api import DoubleVector as dv
from shyft.api import utctime_now
from shyft.api import quantile_map_forecast


class QuantileMapping(unittest.TestCase):
    def test_forecast(self):
        fx_avg = ts_point_fx.POINT_AVERAGE_VALUE
        utc = Calendar()
        ta = TimeAxis(utc.time(2017, 1, 1, 0, 0, 0), deltahours(24), 4)
        historical_data = TsVector()

        forecast_sets = TsVectorSet()
        weight_sets = dv()
        num_historical_data = 56

        # Let's make three sets, one of two elements, one of three, and one of
        # four.

        forecasts_1 = TsVector()
        forecasts_2 = TsVector()
        forecasts_3 = TsVector()

        forecasts_1.append(TimeSeries(ta, dv([13.4, 15.6, 17.1, 19.1]), fx_avg))
        forecasts_1.append(TimeSeries(ta, dv([34.1, 2.40, 43.9, 10.2]), fx_avg))
        forecast_sets.append(forecasts_1)
        weight_sets.append(5.0)

        forecasts_2.append(TimeSeries(ta, dv([83.1, -42.2, 0.4, 23.4]), fx_avg))
        forecasts_2.append(TimeSeries(ta, dv([15.1, 6.500, 4.2, 2.9]), fx_avg))
        forecasts_2.append(TimeSeries(ta, dv([53.1, 87.90, 23.8, 5.6]), fx_avg))
        forecast_sets.append(forecasts_2)
        weight_sets.append(9.0)

        forecasts_3.append(TimeSeries(ta, dv([1.5, -1.9, -17.2, -10.0]), fx_avg))
        forecasts_3.append(TimeSeries(ta, dv([4.7, 18.2, 15.3000, 8.9]), fx_avg))
        forecasts_3.append(TimeSeries(ta, dv([-45.2, -2.3, 80.2, 71.0]), fx_avg))
        forecasts_3.append(TimeSeries(ta, dv([45.1, -92.0, 34.4, 65.8]), fx_avg))
        forecast_sets.append(forecasts_3)
        weight_sets.append(3.0)

        for i in range(num_historical_data):
            historical_data.append(TimeSeries(ta, dv.from_numpy(np.random.random(ta.size()) * 50.0), fx_avg))

        # need one more exposed from core here: auto historical_order = qm::quantile_index<tsa_t>(historical_data, ta);

        interpolation_start = ta.time(2)
        interpolation_end = ta.time(3)
        # Act
        result = quantile_map_forecast(forecast_sets, weight_sets, historical_data, ta, interpolation_start, interpolation_end, False)

        self.assertIsNotNone(result)
        self.assertEqual(len(result), num_historical_data)
        for ts in result:
            self.assertEqual(ts.size(), ta.size())

    def test_speed(self):
        """
        the purpose of this test is to figure out the
        speed characteristics of qm.
        testcases of interest is
          12 (4x3)  forecast typical arome 4 times a day, use last 3 days, 1-3 hour dt, 14 days ahead
          100 historical scenarios
          time-axis wanted is
          next 14 days plus historical scenarios 3..60 weeks ahead

        expected performance:
          the first 14 days includes sorting 10 forcasts, 100 historical pr. timestep.
           the period after should be close to 'memcpy' performance.

        """
        # Arrange the inputs

        utc = Calendar()
        n_hist_ts = 100
        n_fc_days = 14
        n_hist_days = n_fc_days + 360
        n_fc_ts = 10
        t0 = utc.time(2017, 1, 1, 0, 0, 0)

        def generate_ts(ta: TimeAxis, n_fc: int) -> TsVector:
            fx_avg = ts_point_fx.POINT_AVERAGE_VALUE
            r = TsVector()
            w = 2 * 3.14 / len(ta)
            for i in range(n_fc):
                a = np.random.ranf() * 20 - 10.0
                b = np.random.ranf() * 5.0
                v = dv([a + b * math.sin(w * i) for i in range(len(ta))])
                r.append(TimeSeries(ta, v, fx_avg))
            return r

        ta_hist = TimeAxis(t0, deltahours(1), 24 * n_hist_days)
        historical_scenario_ts = generate_ts(ta_hist, n_hist_ts)
        fc_set = TsVectorSet()
        fc_weight = dv()
        n_fc_sets = 4 * 2
        fc_every_dt = deltahours(6)  # six hours between each arome fc.
        dt_fc = deltahours(1)
        for i in range(n_fc_sets):
            t0_fc = t0 + fc_every_dt * i
            fc_set.append(generate_ts(TimeAxis(t0_fc, dt_fc, 24 * n_fc_days), n_fc_ts))
            fc_weight.append(float(3 + i))

        # interpolation_start= no_utctime
        # Act
        interpolated_quantiles = False
        qm_end_idx1 = 24 * (n_fc_days - 2)
        qm_end_idx2 = 24 * (n_fc_days - 1)
        n_ts = 0
        n_v = 0
        print(r"n_days\ttime_used[s]\n")
        tot_seconds = 0.0
        for h_days in range(n_fc_days + 10, n_hist_days, 30):
            ta_qm = TimeAxis(t0 + n_fc_sets * fc_every_dt, dt_fc, 24 * h_days)
            a0 = utctime_now()
            result = quantile_map_forecast(fc_set, fc_weight, historical_scenario_ts, ta_qm, ta_qm.time(qm_end_idx1), ta_qm.time(qm_end_idx2), interpolated_quantiles)
            self.assertIsNotNone(result)
            n_ts += len(result)
            n_v += len(result[0]) * len(result)
            a1 = utctime_now()
            tot_seconds += float(a1 - a0)
            print(f' {h_days}\t{float(a1-a0)}')
        print(f'Total of {n_ts} ts, was forecasted, number values produced {n_v/1000000} Mpts, Mb/s={8.0*n_v/1000000/tot_seconds}')
