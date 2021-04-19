"""White cell infiltrate may be related to tumour progression.

Data from infiltrate is contained in a standalone file. This script was
designed to fuse it with a crystal. While it fuses the crystal, it computes
the frequency table.

It outputs a new sheet containing the frequency table plus the lymphocyte
information stitched together.

This script was made ad-hoc for a single application. I'm not sure it can be
generally useful.
"""

from pathlib import Path

import pandas as pd


def stitch_lymp(file_in, reference, file_out=None):
    """Retrieve and fuse together white cell infiltrate with a crystal.

    Parameters
    ----------
    file in : str
        Full path to crystal file containing patient mutations.
    reference : str
        Full path to file containing lymphocyte infiltrate information.
    file_out : str, optional
        Full path to output file. If unspecified, does not write a file.

    Returns
    -------
    Pandas dataframe containing the frequency table together with lympocyte
    infiltrate information.

    Dependencies
    ------------
    Pandas -- v. 0.25.1
    Python -- v. 3.7.4
    """
    # Handle Paths
    file_in = Path(file_in)
    if file_out is not None:
        file_out = Path(file_out)
    reference = Path(reference)
    # Read data
    crystal = pd.read_csv(file_in)
    crystal = crystal.rename(columns={"Identifier": "ParticipantBarcode"})
    database = pd.read_excel(reference)
    # Standardize ID in the crystal
    old_ids = crystal["ParticipantBarcode"]
    new_ids = []
    for index, item in old_ids.iteritems():
        new_ids.append(item[:12])  # Hold the first 11 characters
    crystal.ParticipantBarcode = new_ids
    # Compute a frequency table
    crystal_freq = pd.crosstab(
        index=crystal["ParticipantBarcode"], columns="Mutation Freq"
    )
    # Merging sheets
    merged = crystal_freq.merge(database, on="ParticipantBarcode")
    if file_out is not None:
        merged.to_csv(file_out)
    return merged


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("file_in", help="Target file")
    parser.add_argument(
        "database", help=("Database file containing" "lymphocyte information")
    )
    parser.add_argument(
        "-o",
        "--out",
        help="Specify target output file. Defaults to writing no file.",
        nargs="?",
        default=None,
    )

    args = parser.parse_args()
    stitch_lymp(file_in=args.file_in, reference=args.database, file_out=args.out)
