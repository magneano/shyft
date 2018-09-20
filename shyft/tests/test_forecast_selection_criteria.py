import unittest
from shyft.api import time, UtcPeriod, Calendar
from shyft.repository.interfaces import ForecastSelectionCriteria, ForecastSelectionCriteriaError


class ForecastSelectionCriteriaTestCase(unittest.TestCase):

    def test_construct_empty_raises(self):
        """ force user to supply at least one criteria """
        with self.assertRaises(ForecastSelectionCriteriaError):
            fc = ForecastSelectionCriteria()

    def test_construct_wrong_keyword_raises(self):
        """ force user to supply at least one criteria """
        with self.assertRaises(ForecastSelectionCriteriaError):
            fc = ForecastSelectionCriteria(any_key_word_that_is_wrong=UtcPeriod())

    def test_utcperiod_criteria(self):
        test_period = UtcPeriod(0, 3600)
        fc = ForecastSelectionCriteria(forecasts_created_within_period=test_period)
        self.assertTrue(fc.criterion[0] == 'forecasts_created_within_period')
        self.assertTrue(fc.criterion[1] == test_period)
        with self.assertRaises(ForecastSelectionCriteriaError):
            ForecastSelectionCriteria(forecasts_created_within_period='should throw')

    def test_latest_available_fc(self):
        fc = ForecastSelectionCriteria(latest_available_forecasts={'number_of_forecasts': 1, 'forecasts_older_than': Calendar().time(2018, 1, 20, 7)})
        self.assertTupleEqual(('latest_available_forecasts', {'forecasts_older_than': 1516431600, 'number_of_forecasts': 1}), fc.criterion)
        with self.assertRaises(RuntimeError):
            ForecastSelectionCriteria(latest_available_forecasts={'number_of_forecasts': 1, 'forecasts_older_than': 'text that is wrong'})
        with self.assertRaises(ForecastSelectionCriteriaError):
            ForecastSelectionCriteria(latest_available_forecasts={'number_of_forecasts': 'a string 1', 'forecasts_older_than': 3600.0 })

    def test_forecasts_that_cover_period(self):
        test_period = UtcPeriod(0, 3600)
        fc = ForecastSelectionCriteria(forecasts_that_cover_period=test_period)
        self.assertTrue(fc.criterion[0] == 'forecasts_that_cover_period')
        self.assertTrue(fc.criterion[1] == test_period)
        with self.assertRaises(ForecastSelectionCriteriaError):
            ForecastSelectionCriteria(forecasts_that_cover_period='hmm throw')

    def test_forecasts_that_intersect_period(self):
        test_period = UtcPeriod(0, 3600)
        fc = ForecastSelectionCriteria(forecasts_that_intersect_period=test_period)
        self.assertTrue(fc.criterion[0] == 'forecasts_that_intersect_period')
        self.assertTrue(fc.criterion[1] == test_period)
        with self.assertRaises(ForecastSelectionCriteriaError):
            ForecastSelectionCriteria(forecasts_that_intersect_period='xhmm throw')

    def test_forecasts_at_reference_times(self):
        test_ref = [0, 3600, 7200]
        fc = ForecastSelectionCriteria(forecasts_at_reference_times=test_ref)
        self.assertTrue(fc.criterion[0] == 'forecasts_at_reference_times')
        for t1, t2 in zip(test_ref, fc.criterion[1]):
            self.assertTrue(time(t1), t2)
        with self.assertRaises(ForecastSelectionCriteriaError):
            ForecastSelectionCriteria(forecasts_that_intersect_period='xhmm throw')
