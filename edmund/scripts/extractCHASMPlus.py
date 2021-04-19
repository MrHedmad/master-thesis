"""This script is used to extract and filter at the same time data contained
in a CHASM Plus output file. Please see the website for CHASM Plus for
more information on it.
"""

import logging
import os
import pickle
from io import StringIO
from itertools import zip_longest
from pathlib import Path
from shutil import rmtree
from typing import Any, Iterator, Optional, Union

import click
import pandas as pd
from statsmodels.stats.multitest import multipletests

from edmund.entrypoint import cli
from edmund.scripts.dumpClinical import call_portal
from edmund.scripts.fuse_csvs import fuse_csvs
from edmund.utils import c, count_missing_rows, cqdm, uniques

log = logging.getLogger(__name__)

# We get a warning in extractor_plus when the dataframe column is overwritten
# The warning tells about potential problems but it is a false positive
# This option gets rid of it.
pd.options.mode.chained_assignment = None  # default='warn'


def strip_unnamed(string: str) -> str:
    """When fusing headers coming from excel, the lines with NA values are
    named as `Unnamed:...`. This function filters them out.
    """
    if string.strip().startswith("Unnamed:"):
        return ""
    else:
        return string.strip()


def combine_headers(head1: str, head2: str, sep: str) -> str:
    """Combines two header strings (with several headers) into just one.

    If some headers are missing, uses the last non-zero header instead.

    Args:
        head1: SEP-separated headers in a single string.
        head2: SEP-separated headers in a single string.
        sep: Separators used to read and write the headers.
    """
    head1 = head1.split(sep)
    head2 = head2.split(sep)
    casted = []
    last_nonzero_x = ""
    last_nonzero_y = ""
    for x, y in zip_longest(head1, head2, fillvalue=""):
        last_nonzero_x = x.strip() or last_nonzero_x
        last_nonzero_y = y.strip() or last_nonzero_y
        casted.append(last_nonzero_x + "_" + last_nonzero_y)
    return sep.join(casted)


def simplify_headers(headers: list[str]) -> list[str]:
    """Simplify the headers for the data to be more readable

    The header names are fusions between the first and second header line
    in the OpenCravat input. This function removes the parts of the header
    names that are useless. In short, it keeps only the names of the tools.
    """
    to_rename = ["Variant Annotation", "Hg19", "Tag Sampler", "VCF Info"]
    new_headers = []
    for header in headers:
        split = header.split("_")
        if split[0] in to_rename:
            new_headers.append(split[1])
        else:
            new_headers.append(header)
    return new_headers


def input_caster(file_path: Path, sep: str = "auto") -> pd.DataFrame:
    """Casts an input CHASM file into a pandas dataframe

    This is needed as inputs can be written as csv, tsv, or xlsx.
    Only retains the mutation-level data. Gene level and other outputs are stripped.

    Args:
        file_path: The path to the input file
        sep: The separator of columns in the file encoding. If "xlsx" will ingore it
            and treat it as an excel file. If "auto", tries to determine it from
            filetype.

    Returns:
        A pandas dataframe with the cleaned input.
    """
    if sep == "auto":
        if file_path.suffix == ".tsv":
            sep = "\t"
        elif file_path.suffix == ".csv":
            sep = ","
        elif file_path.suffix in [".xlsx", ".xls"]:
            # We need to treat this differently
            sep = "excel"
    elif sep in ["xlsx", "xls"]:
        sep = "excel"

    if sep == "excel":
        # Get page 2 of the excel file (The one with the mutations)
        excel: pd.DataFrame = pd.read_excel(file_path.resolve(), 1)
        # Handle the horrible headers
        head1 = excel.columns
        head2 = excel.values[0]
        casted = []
        last_nonzero_x = ""
        last_nonzero_y = ""
        for x, y in zip_longest(head1, head2, fillvalue=""):
            last_nonzero_x = strip_unnamed(x) or last_nonzero_x
            last_nonzero_y = strip_unnamed(y) or last_nonzero_y
            casted.append(last_nonzero_x + "_" + last_nonzero_y)
        excel.columns = casted
        # Drop the first row that was once a header
        excel = excel.drop(0)
        df = excel
    else:
        with file_path.open("r") as file:
            # The first frame is for the variants, so we want that
            # Ignore the first lines
            first_header = ""
            for line in file:
                if line.startswith("#"):
                    continue
                first_header = line
                break

            # The first two lines are headers. We need to combine them
            header = combine_headers(first_header, file.readline(), sep)
            with StringIO(header) as io:
                io.write(header + "\n")
                for line in file:
                    if line.startswith("#"):
                        break
                    io.write(line)
                io.seek(0)
                df = pd.read_csv(io, sep=sep)
    # Rename tumor-specific columns
    cols = []
    for name in df.columns:
        if name.startswith("CHASMplus_") and len(name.split("_")) == 3:
            name = name.split("_")
            name = name[0] + "_tspec_" + name[2]
        cols.append(name)
    df.columns = cols

    return df


def extractor_plus(
    dataframe: pd.DataFrame,
    data_id: Optional[str] = None,
    rename_hugo: bool = False,
    adjust_pval: bool = True,
    drop_na_cols: Union[str, bool, list[str]] = "default",
) -> tuple[pd.DataFrame, int]:
    """Extract data from a CHASM Plus output DataFrame.

    The column names to be used are hardcoded in, so it'll not work when
    using another predictor tools inside the OpenCravat report (KeyError).

    Args:
        dataframe: The input pandas dataframe to process read by `input_caster`
        data_id: Optional ID for increased information in error logging. Will
            be included in errors related to the handling of the data.
        rename_hugo : Do HUGO symbols need to be renamed?
        adjust_pval : Whether to adjust p-values. Appends adjusted p-values to
            a new column named "CHASMplus_FDR" and "Tum-Spec FDR".
            Uses Benjamini-Hockberg corrections.
        drop_na_cols : Specify column(s) to purge from NAs as a list of
            strings with the names of the columns. If False, will not drop any lines.
            Pass "default" to filter both pancancer and tumor-specific p-values
            Warning: Calculating FDRs without filtering p-values will probably
            lead to NaN FDRs!

    Returns:
        The extracted dataframe.

    Raises:
        KeyError: If the key specified for renaming HUGOs is invalid.
        EmptyFileError: if the file is empty after dropping NAs (see drop_na_cols
            argument).
    """
    if data_id:
        log.info(f"Processing {data_id}")
    else:
        log.info("Processing data with no identifier")

    if drop_na_cols == "default":
        drop_na_cols = ["CHASMplus_P-value", "CHASMplus_tspec_P-value"]

    log.debug(
        f"Renaming and selecting columns. Initial column names: {dataframe.columns}"
    )

    dataframe.columns = simplify_headers(dataframe.columns)

    log.debug(f"Colnames after renaming: {dataframe.columns}")

    # ----- Check for missing values ----------------------------------------
    if drop_na_cols is not False:
        log.debug("Checking for missing values.")
        missing = count_missing_rows(dataframe, col=drop_na_cols)
        if missing == 0:
            log.debug("Found no missing values.")
        else:
            log.debug(
                f" Found {missing} missing value(s). Dropping lines containing them."
            )
            dataframe.dropna(inplace=True, subset=drop_na_cols)

        if len(dataframe.index) == 0:
            log.error(f"File {data_id} has no surviving rows after dropping NAs")
            raise EmptyFileError(
                f"File {data_id} has no surviving rows after dropping NAs."
            )

    # ------ Adjust P-values -------------------------------------------------
    if adjust_pval is True:
        log.debug("Adjusting P-values using BH correction.")
        if "CHASMplus_P-value" in dataframe.columns:
            p_values = dataframe["CHASMplus_P-value"].tolist()
            adj_p_values = multipletests(p_values, method="fdr_bh")
            dataframe["CHASMplus_FDR"] = adj_p_values[1]

        if "CHASMplus_tspec_P-value" in dataframe.columns:
            p_values = dataframe["CHASMplus_tspec_P-value"].tolist()
            adj_p_values = multipletests(p_values, method="fdr_bh")
            dataframe["CHASMplus_tspec_FDR"] = adj_p_values[1]

        if (
            "CHASMplus_P-value" not in dataframe.columns
            and "CHASMplus_tspec_P-value" not in dataframe.columns
        ):
            log.warn(
                "P-values need to be adjusted, but none were included."
                " Skipping P-value adjustment."
            )

    # ----- Rename HUGO symbols if specified in the input --------------------
    # Create the variable for HUGO rename, if needed.
    if isinstance(rename_hugo, bool) is False:
        try:
            hugo_names = dataframe[rename_hugo].tolist()
        except KeyError:
            raise KeyError(
                (
                    "Wrong key given to rename HUGOs. The column was"
                    " possibly not included."
                )
            )
    else:
        hugo_names = None

    if "Hugo" in list(dataframe) and rename_hugo is not False:
        log.debug("HUGO symbols are being renamed.")
        # Rename consecutive if True
        if rename_hugo is True:
            log.debug("Using consecutive method.")
            # Values are sorted as it is required to rename correctly.
            # See _rename_list
            dataframe.sort_values("Hugo")
            renamed_hugo = _rename_list(dataframe["Hugo"], method="consecutive")
            dataframe["Renamed HUGO"] = renamed_hugo
        # Rename using a column
        elif hugo_names is not None:
            log.debug(f"Fusing with {rename_hugo}.")
            renamed_hugo = _rename_list(
                dataframe["Hugo"], list2=hugo_names, method="fuse"
            )
            dataframe["Renamed HUGO"] = renamed_hugo

    elif "Hugo" not in dataframe.columns and rename_hugo is not False:
        log.warn("HUGO need to be renamed, but not included. Skipping.")

    log.debug("ExtractCHASM completed.")
    return dataframe


def _rename_list(
    list1: list[str], method: str, list2: Optional[list[str]] = None
) -> list[str]:
    """Renames multiple list items following different methods.

    Method can be "consecutive" to rename identical consecutive items
    with increasing number or "fuse" to fuse with the second list (list2)

    Args:
        list1: List of strings to be renamed
        method: Type of renaming method. Can be "consecutive" to rename consecutive
            identical items as progressive values ("tree", "tree" will become
            "tree_0", "tree_1"); "fuse" to fuse with the values from list2.
            list2 must be of the same length as list1 while fusing.

    Returns:
        list of renamed values
    """
    last = None
    new_list = []
    count = 0
    if method == "consecutive":
        for item in list1:
            if item == last:
                count += 1
                new_list.append(item + "_" + str(count))
                last = item
            else:
                new_list.append(item + "_0")
                count = 0
                last = item

    elif method == "fuse" and list2 is not None:
        if len(list1) != len(list2):
            raise ValueError(
                ("Length of lists not equal:" f" {len(list1)} vs {len(list2)}")
            )
        new_list = [str(i1) + "_" + str(i2) for i1, i2 in zip(list1, list2)]
    else:
        raise ValueError("Unrecognized method or list2 unspecified")
    return new_list


class EmptyFileError(Exception):
    pass


def absolute_file_paths(directory: Union[str, Path]) -> Iterator[Path]:
    """Grab full paths for all files in a directory.

    Args:
        directory : Path to target directory.

    Returns:
        Iterator with each full path as a `Path`
    """
    directory = Path(directory)
    item: Path
    for item in directory.iterdir():
        if item.is_file():
            yield item.absolute()


def extract_project(
    folder_in: Union[Path, str],
    folder_out: Union[Path, str],
    project_name: Optional[str] = None,
    keep_temp: bool = True,
) -> None:
    """Function that extracts all files in a folder, combines them into a
    crystal, retrieves clinical data, and creates a clinical frame.

    Args:
        folder_in: Path to folder containing the data to extract. File type
            is detected automatically at extraction time. The files
            are expected to start with the patient ID as found in the
            clinical data.
        folder_out: Path to folder with output data. Temporary files
            will also be stored here, under a `temp` folder.
        project_name: TCGA project name. If unspecified, uses the folder name.
            For convenience, the folder name is first split at any underscores,
            and the first part of the name is then used. For instance,
            `BRCA_TSV` is treated as `BRCA`.
        keep_temp: Should the temp folder be kept? If False, the temp folder
            is deleted after extraction.
    """
    log.debug("Processing input")
    folder_in = Path(folder_in)
    folder_out = Path(folder_out)
    os.makedirs(folder_out, exist_ok=True)

    temp_folder = folder_out / "temp"
    log.info(f"Temp folder set to be {temp_folder}")
    os.makedirs(temp_folder, exist_ok=True)

    if project_name is None:
        project_name = folder_in.stem.split("_")[0]

    # I need just the tumor type, so i strip the TCGA if it's present
    if not project_name.startswith("TCGA-"):
        project_name = f"TCGA-{project_name}"
    log.info(f"Retrieving clinical data for project: {project_name}")
    clinical_data = call_portal(project_name)

    log.info("Extracting files")
    for path in cqdm(
        folder_in.iterdir(),
        total=len(list(folder_in.iterdir())),
        desc="Extracting files",
    ):
        patient_id = path.name[:12]
        try:
            extracted_frame = extractor_plus(
                dataframe=input_caster(file_path=path),
                data_id=path.name,
                rename_hugo="Protein Change",
                adjust_pval=True,
            )
        except EmptyFileError:
            log.warning(f"File {path.name} was empty after filtering NAs. Ignoring it.")
            continue

        extracted_frame.to_csv(temp_folder / (patient_id + ".csv"))
    log.info("Fusing extracted files")
    crystal = fuse_csvs(temp_folder, recursive=False)

    clinical_data = clinical_data.loc[
        clinical_data["submitter_id"].isin(crystal["Identifier"])
    ]
    crystal.to_csv(folder_out / f"{project_name}_crystal.csv")
    clinical_data.to_csv(folder_out / f"{project_name}_clinical.csv")

    if keep_temp is False:
        rmtree(temp_folder, ignore_errors=True)


@cli.command(name="extract")
@click.argument("target-path")
@click.argument("output-path")
@click.option(
    "--keep-loose",
    default=False,
    is_flag=True,
    help="Force keeping intermediate files. Deletes them by default.",
)
def project_command(
    target_path: Union[Path, str],
    output_path: Union[Path, str],
    keep_loose: bool,
):
    """Process Open CRAVAT CHASMplus predictions and make R-ready csv files

    Processes all files in target folder TARGET_PATH. Looks recursively in the
    folder, finding subfolders. Expects that these folders are named with the
    TCGA tumor identifier. For instance, TCGA-BRCA tumor files should be in a
    folder named `BRCA`. Writes all outputs in similarly named folders to
    OUTPUT_PATH.

    For simplicity, strips from the subfolder names any character after `_`.
    So, for instance, `BRCA_FILES` is treated as `BRCA` when detecting the
    tumor ID.

    The single files containing patient data should be named starting with
    the patient unique ID. These are used in order to create data frames
    with clinical data relevant to the patients in the csv.
    """
    target_path = Path(target_path)
    output_path = Path(output_path)

    log.info("Finding target folders in target path.")
    names: dict = {}
    for subfolder in target_path.iterdir():
        if not subfolder.is_dir():
            continue
        tumor_id = subfolder.name.split("_")[0]
        names[f"TCGA-{tumor_id}"] = subfolder

    long_names = ", ".join(names.keys())
    log.info(f"Found {len(names)} folders: {long_names}")

    for project_id, path in cqdm(names.items(), desc="Generating Projects"):
        log.info(f"Making output path for {project_id}")
        project_output_path = output_path / project_id
        os.makedirs(project_output_path, exist_ok=True)

        log.info(f"Generating the project for {project_id}")

        extract_project(path, output_path / project_id, project_id, keep_loose)
