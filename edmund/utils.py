"""A collection of tools useful in the other scripts. The usual dependency
is just Python 3.7.4
"""
import logging
import math
import os
from pathlib import Path
from typing import Iterator, Mapping

from tqdm import tqdm

from edmund.config import CONFIG

log = logging.getLogger(__name__)


def get_line(file: str, number: int) -> str:
    """Get and return the i-th line in a file.

    Remember that python begins to count at 0, so the 10th line will be
    number = 9. The line is returned without the ending "\\n".

    Parameters
    ----------
    file : str
        Full path to target file.
    number : int
        Which line to return.

    Returns
    -------
    str containing the i-th line of the file.
    """
    file = Path(file)
    with file.open() as file:
        for i, line in enumerate(file):
            if i == number:
                return line.split("\n")[0]


def absolute_file_paths(directory: str) -> Iterator[str]:
    """Grab full paths for all files in a directory.

    It will return an iterable containing every not-folder full paths
    as strings.

    Parameters
    ----------
    directory : str
        Full path to target directory.

    Returns
    -------
    iterator with each full path as a string.
    """
    directory = Path(directory)
    for item in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, item)):
            yield os.path.abspath(os.path.join(directory, item))


def count_files(directory: str) -> int:
    """Count how many files are in a directory, ignores subdirectories.

    Parameters
    ----------
    directory:
        Full path to target directory.

    Returns
    -------
    int representing the number of non-directory files in the folder.
    """
    directory = Path(directory)
    i = 0
    for item in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, item)):
            i += 1
    return i


def uniques(list_in):
    """Return the unique values from a list.

    Passes through a dictionary to remove duplicate entries, then returns
    it back into a list.
    """
    return list(dict.fromkeys(list_in))


def tuple_unique_permutations(list_in):
    """Return the unique values from a list of tuples.

    Passes through a dictionary to remove duplicate entries, then returns
    it back into a list.
    """
    return list(set(tuple(sorted(l)) for l in list_in))


def drop_na(list_in):
    """Remove all na or nan from a list"""
    return [x for x in list_in if str(x).lower() not in ["nan", "na"]]


def verb(t=False, *args):
    """Small utility to write lines if verbose is True.

    Parameters
    ----------
    t : bool
        Whether to write or not to the console.
    Other args:
        Parsed to the print command.
    """
    if t is True:
        print(args[0], end=args[1])


def make_gauss(sigma=1, mu=0):
    """Returns a Gaussian function with the specified parameters.

    Precomputes some parameters in order to speed up usage.

    Parameters
    ----------
    sigma : float
        The standard deviation of the gaussian.
    mu : float
        The mean value of the gaussian.

    Returns
    -------
    A function f(x) that models the gaussian.
    """
    k = 1 / (sigma * math.sqrt(2 * math.pi))
    s = -1.0 / (2 * sigma * sigma)

    def f(x):
        return k * math.exp(s * (x - mu) * (x - mu))

    return f


def count_missing_rows(data_frame, col="any") -> int:
    """Count rows containing missing values in a Pandas DataFrame.

    Parameters
    ----------
    data_frame : Pandas DataFrame
        A Pandas DataFrame to count NA/NANs
    col : str
        The column to count NA/NANs from. If "any" is specified, considers
        any NA in any column as a hit.

    Returns
    -------
    int representing the number of missing values counted.
    """
    i = 0
    if col == "any":
        for index, row in data_frame.iterrows():
            if any(row.isna()) is True:
                i += 1
    else:
        for index, row in data_frame[col].iterrows():
            if any(row.isna()) is True:
                i += 1
    return i


# Yanked from the python discord bot
def recursive_update(original: Mapping, new: Mapping) -> Mapping:
    """Recursively update nested dictionaries

    Helper method which implements a recursive `dict.update`
    method, used for updating the original configuration with
    configuration specified by the user.

    Args:
        original: A mapping to be updated recursively.
        new: A mapping used to update the original with, recursively.

    Returns:
        The updated mapping.
    """
    for key, value in original.items():
        if key not in new:
            continue

        if isinstance(value, Mapping):
            if not any(isinstance(subvalue, Mapping) for subvalue in value.values()):
                original[key].update(new[key])
            recursive_update(original[key], new[key])
        else:
            original[key] = new[key]

    return original


def inherit_docs(original):
    """A decorator that copies the docstring from a function to the decorated one"""

    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper


@inherit_docs(tqdm)
def cqdm(*args, **kwargs):
    """Like tqdm, but obey configuration"""
    return tqdm(*args, **kwargs, disable=not CONFIG["show_loading"])


def c(*args: str) -> str:
    """Joins all strings together"""
    return "".join(args)
