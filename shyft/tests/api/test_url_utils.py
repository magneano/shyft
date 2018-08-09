import unittest
from shyft.api import (urlencode, urldecode, extract_shyft_url_query_parameters, extract_shyft_url_path, extract_shyft_url_container, shyft_url)


class UrlUtilsTestCase(unittest.TestCase):
    """ Verify the python exposure for url-utils """

    def test_url_encode_decode(self):
        for s in ['', 'a', ' ', 'This is æøåÆØÅ', 'Something Special %12!']:
            self.assertEqual(urldecode(urlencode(s)), s)
            self.assertEqual(urldecode(urlencode(s, False), False), s)

    def test_shyft_url(self):
        self.assertEqual(shyft_url('test', 'a/ts/path.db'), 'shyft://test/a/ts/path.db')
        self.assertEqual(shyft_url('test', 'a/ts/path.db', {'a': '123'}), 'shyft://test/a/ts/path.db?a=123')

    def test_extract_shyft_url_container(self):
        self.assertEqual('test', extract_shyft_url_container('shyft://test/something'))

    def test_extract_shyft_url_path(self):
        self.assertEqual('something/strange/here.db', extract_shyft_url_path('shyft://test/something/strange/here.db'))

    def test_extract_shyft_query_parameters(self):
        self.assertEqual({'a': '123'}, extract_shyft_url_query_parameters('shyft://test/something/strange/here.db?a=123'))
