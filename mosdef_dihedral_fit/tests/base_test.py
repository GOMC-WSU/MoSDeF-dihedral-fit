import numpy as np
import pytest
import pathlib


class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    def get_fn(pathname):
        """Get test file path in test/files"""
        current_path = pathlib.Path(__file__).parent.resolve()
        return str(current_path / pathname)


