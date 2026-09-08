import os
import tempfile
import unittest

from proteindf_bridge.functions import get_yaml, load_yaml, save_yaml


class FunctionsTest(unittest.TestCase):
    def test_save_yaml_and_load_yaml_roundtrip(self):
        data = {"a": 1, "b": [1, 2, 3], "c": "日本語"}
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "test.yaml")
            save_yaml(data, path)
            loaded = load_yaml(path)
            self.assertEqual(loaded, [data])

    def test_get_yaml_returns_str(self):
        yaml_str = get_yaml({"a": 1})
        self.assertIsInstance(yaml_str, str)


if __name__ == '__main__':
    unittest.main()
