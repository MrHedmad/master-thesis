import re
import unittest
import pandas as pd
from pandas.testing import assert_frame_equal
import numpy as np
from tests import TEST_DATA

from edmund.scripts import crystallize


class TestCrystallization(unittest.TestCase):
    testfiles = TEST_DATA / "crystallize"

    def test_correct_fusion(self):
        expected_result = pd.DataFrame(
            data={
                "col1": ["a", "d", "y"],
                "col2": ["b", "e", "x"],
                "col3": ["c", np.NaN, np.NaN],
                "col4": [np.NaN, "f", "z"],
            }
        )
        result = crystallize.fuse_csvs(self.testfiles)
        assert_frame_equal(
            result.sort_values(axis=0, by="col1", ignore_index=True),
            expected_result.sort_values(axis=0, by="col1", ignore_index=True),
            check_like=True,
        )

    def test_no_recursion(self):
        expected_result = pd.DataFrame(
            data={
                "col1": [
                    "a",
                    "d",
                ],
                "col2": [
                    "b",
                    "e",
                ],
                "col3": ["c", np.NaN],
                "col4": [np.NaN, "f"],
            }
        )
        result = crystallize.fuse_csvs(self.testfiles, recursive=False)
        assert_frame_equal(
            result.sort_values(axis=0, by="col1"),
            expected_result.sort_values(axis=0, by="col1"),
            check_like=True,
        )

    def test_regex_pattern(self):
        with self.assertRaises(RuntimeError):
            crystallize.fuse_csvs(self.testfiles, pattern=re.compile(".bnna$"))

    def test_max_recursion(self):
        with self.assertRaises(RecursionError):
            crystallize.fuse_csvs(self.testfiles, max_depth=0)
