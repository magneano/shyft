import unittest
from shyft.api import win_short_path,win_set_priority
import os


class TestWinShortPath(unittest.TestCase):
    def test_basic_translate(self):
        long_path = r"C:\Program Files"
        s = win_short_path(long_path)
        if os.name == 'nt':
            self.assertEqual(s, r'C:\PROGRA~1')
            s2 = win_short_path(r"C:\Program File")
            self.assertFalse(s2)
        else:
            self.assertEqual(s, long_path)

    def test_win_set_priority(self):
        win_set_priority(-1)  # below normal and no io/mem performance
        win_set_priority(0)  # normal, hopefully a sane setting
