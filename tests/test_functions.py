import os
import tempfile
import unittest

from proteindf_bridge.functions import (
    get_yaml,
    load_yaml,
    save_yaml,
    load_msgpack,
    save_msgpack,
)


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

    def test_save_msgpack_and_load_msgpack_roundtrip(self):
        force_data = [{0: 1.0, 1: 2.0, 2: 3.0}]
        data = {
            "force": force_data,
            b"bytes_key": b"bytes_val",
            100: {0: "x", 1: "y"},
        }
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "test.mpac")
            save_msgpack(data, path)
            loaded = load_msgpack(path)

            expected = {
                "force": [{0: 1.0, 1: 2.0, 2: 3.0}],
                "bytes_key": "bytes_val",
                100: {0: "x", 1: "y"},
            }
            self.assertEqual(loaded, expected)
            self.assertIn(0, loaded["force"][0])
            self.assertNotIn("0", loaded["force"][0])
            self.assertEqual(loaded["force"][0][0], 1.0)
            self.assertIn(100, loaded)
            self.assertIn(0, loaded[100])
            self.assertIn("bytes_key", loaded)


if __name__ == '__main__':
    unittest.main()
