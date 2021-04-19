"""Fuse together csv files into a single, larger file for better manipulation
with other programs, such as R.
"""

import logging
import math
import re
from numbers import Number
from pathlib import Path
from typing import Union

import click
import pandas as pd

from edmund.entrypoint import cli
from edmund.utils import cqdm

log = logging.getLogger(__name__)


def fuse_csvs(
    target_directory: Union[str, Path],
    pattern: re.Pattern = re.compile(".csv$"),
    recursive: bool = True,
    max_depth: Number = math.inf,
) -> pd.DataFrame:
    """Fuse .csv files in target folder(s) into a single pandas object

    Adds an additional column with the name of the file from where the entry
    was taken from. The column name for this is `Identifier`

    It assumes that the header is the first line of each file.

    Args:
        target_directory : Target directory containing .csv files to fuse together.
        pattern: A compiled regular expression (from module `re`) that the filename
            must match for it to be included in the fusion csv. Defaults to files
            that end with ".csv".
        recursive: If true, the script also looks inside folders in the target
            path recursively.
        max_depth: Raises an error if the recursive depth of searching exceeds
            this value. Defaults to no limit

    Returns:
        A pandas object with the concatenated frame.
    """
    path: Path = Path(target_directory)
    objects: list = []

    def find_objects(
        path: Path,
        pattern: re.Pattern,
        depth: int,
        max_depth: int,
        recurse: bool = True,
    ):
        if depth > max_depth:
            raise RecursionError("Exceeded max search depth")
        assert path.is_dir()

        resulting_paths: list = []
        files = cqdm(
            path.iterdir(),
            total=len(list(path.iterdir())),
            desc=f"Looking in {path.name}",
        )
        for file in files:
            if recurse and file.is_dir():
                resulting_paths.extend(
                    find_objects(file, pattern, depth + 1, max_depth, recurse)
                )
                continue
            elif file.is_dir():
                continue
            if pattern.search(file.name):
                resulting_paths.append(file)

        return resulting_paths

    file_paths = find_objects(path, pattern, max_depth, max_depth, recursive)

    if len(file_paths) == 0:
        raise RuntimeError("No files found for fusion")

    for file_path in cqdm(file_paths, desc="Reading found files"):
        frame = pd.read_csv(file_path)
        frame["Identifier"] = file_path.stem
        objects.append(frame)

    return pd.concat(
        cqdm(objects, desc="Concatenating objects"), join="outer", ignore_index=True
    )


@cli.command(name="fuse-csvs")
@click.argument("target_directory")
@click.argument("output")
@click.option("--recurse", is_flag=True, help="If set, look recursively in folder")
@click.option(
    "-d",
    "--max_depth",
    default=None,
    type=int,
    help="Specify maximum recursion depth before erroring. Defaults to unlimited",
)
def fuse_csvs_command(target_directory, output, recurse, max_depth):
    """Fuse CSV files in a target folder

    Fuse csv files (that end in .csv) in TARGET_DIRECTORY and write them to
    the OUTPUT file.
    """
    fuse_csvs(target_directory, recursive=recurse, max_depth=max_depth).to_csv(
        Path(output)
    )
