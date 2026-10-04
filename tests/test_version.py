import os
import re
import unittest

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None

import proteindf_bridge
from proteindf_bridge._version import __version__ as version_py


class TestVersionConsistency(unittest.TestCase):
    def test_version_matches_pyproject(self):
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        pyproject_path = os.path.join(root_dir, "pyproject.toml")
        self.assertTrue(os.path.exists(pyproject_path), f"pyproject.toml not found at {pyproject_path}")

        pyproject_version = None
        with open(pyproject_path, "rb" if tomllib else "r", encoding=None if tomllib else "utf-8") as f:
            if tomllib:
                data = tomllib.load(f)
                pyproject_version = data.get("project", {}).get("version")
            else:
                content = f.read()
                match = re.search(r'(?m)^\s*version\s*=\s*"([^"]+)"', content)
                if match:
                    pyproject_version = match.group(1)

        self.assertIsNotNone(pyproject_version, "Could not extract version from pyproject.toml")
        self.assertEqual(
            version_py,
            pyproject_version,
            f"proteindf_bridge/_version.py ({version_py}) does not match pyproject.toml ({pyproject_version})"
        )
        self.assertEqual(
            proteindf_bridge.__version__,
            pyproject_version,
            f"proteindf_bridge.__version__ ({proteindf_bridge.__version__}) does not match pyproject.toml ({pyproject_version})"
        )


if __name__ == "__main__":
    unittest.main()
