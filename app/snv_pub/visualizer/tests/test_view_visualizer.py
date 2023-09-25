from django.test import Client, TestCase, RequestFactory
import json
from logging import DEBUG, getLogger

logger = getLogger(__name__)


class VisualizerViewTest(TestCase):
    def setUp(self):
        # Every test needs a client.
        self.client = Client()

    def test_from_validation(self):
        data = {'score_threshold': 'aaa'}

        response = self.client.post('/visualizer/network-data/', data)
        response_json = response.json()
        str_response_json = json.dumps(response_json)
        response_json_expected = {'internalStatus': {'success': False},
                                  'formValidation': {'invalid_field': 'score_threshold',
                                                     'message': 'Please input a number.'}}

        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(str_response_json, response_json_expected)
