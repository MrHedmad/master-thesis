"""This script is used to extract and filter at the same time data contained
in a CHASM Plus output file. Please see the website for CHASM Plus for
more information on it.

The main function (extractor_plus) is feature rich but complex,
so it is wrapped in a simpler interface function (extract_chasm_plus).
extract_chasm_plus is ran from the command line.
"""

import logging
import os
import pickle
import warnings
from io import StringIO
from itertools import zip_longest
from numbers import Real
from pathlib import Path
from shutil import rmtree
from typing import Any, Iterator, Optional, Union

import click
from numpy import extract
import pandas as pd
from statsmodels.stats.multitest import multipletests

from edmund.entrypoint import cli
from edmund.scripts.dumpClinical import call_portal
from edmund.scripts.fuse_csvs import fuse_csvs
from edmund.utils import count_missing_rows, cqdm, uniques, c

log = logging.getLogger(__name__)

# We get a warning in extractor_plus when the dataframe column is overwritten
# The warning tells about potential problems but it is a false positive
# This option gets rid of it.
pd.options.mode.chained_assignment = None  # default='warn'

# The imputation for the tsv files changes the names of the columns to these
# names. This dict can be used to move between that and the ones needed
# by the extractor.
COL_NAMES_REFERENCES = {
    "Variant Annotation_UID": "UID",
    "Variant Annotation_Chrom": "Chrom",
    "Variant Annotation_Position": "Position",
    "Variant Annotation_Ref Base": "Ref Base",
    "Variant Annotation_Alt Base": "Alt Base",
    "Variant Annotation_Hugo": "Hugo",
    "Variant Annotation_Note": "Note",
    "Variant Annotation_Protein Change": "Protein Change",
    "Variant Annotation_All Mappings": "All Mappings",
    "Variant Annotation_Coding": "Coding",
    "Variant Annotation_Transcript": "Transcript",
    "Variant Annotation_Sequence Ontology": "Ontology",
    "Hg19_Hg19 Chrom": "Hg19 Chrom",
    "Hg19_Hg19 Position": "Hg19 Position",
    "Tag Sampler_Sample Count": "Sample Count",
    "Tag Sampler_Samples": "Samples",
    "Tag Sampler_Tags": "Tags",
    "CHASMplus_P-value": "P-value",
    "CHASMplus_Score": "Score",
    "CHASMplus_Transcript": "Transcript",
    "CHASMplus_All results": "All results",
    "CHASMplus_tspec_P-value": "Tum-Spec P-value",
    "CHASMplus_tspec_Score": "Tum-Spec Score",
    "CHASMplus_tspec_Transcript": "Tum-Spec Transcript",
    "CHASMplus_tspec_All results": "Tum-Spec Results",
    "VCF Info_Phred": "Phred",
    "VCF Info_VCF filter": "VCF filter",
    "VCF Info_Zygosity": "Zygosity",
    "VCF Info_Alternate reads": "Alternate reads",
    "VCF Info_Total reads": "Total reads",
    "VCF Info_Variant AF": "Variant AF",
    "VCF Info_Haplotype block ID": "Haplotype block ID",
    "VCF Info_Haplotype strand ID": "Haplotype strand ID",
}


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
    select: Union[str, dict] = "default",
    rename_hugo: bool = False,
    adjust_pval: bool = True,
    drop_na_cols: Union[str, bool, list[str]] = "default",
) -> tuple[pd.DataFrame, int]:
    """Extract and filter data from a CHASM Plus output DataFrame.

    The column names to be used are hardcoded in, so it'll not work when
    using another predictor tools inside the OpenCravat report (KeyError).

    Args:
        dataframe: The input pandas dataframe to process read by `input_caster`
        data_id: Optional ID for increased information in error logging. Will
            be included in errors related to the handling of the data.
        select : A dictionary with keys being the names of columns to extract.
            The values can be None to not filter the column,
            a string to only keep rows equal to that string, or a list such as
            ["operator", value1, value2].
            "operator" specifies the type of filtering using the value(s).
            Possible operators are:
                - ">" or "more" = Greater than "value1" rows are kept.
                - ">=" or "moreorequal" = Greater than or equal to "value1" rows
                    are kept.
                - "<" or "less" = Less than "value1" rows are kept.
                - "<=" or "lessorequal" = Less than or equal to "value1" rows
                    are kept.
                - "><" or "between" = Rows with values in between "value1" and
                    "value2" (included) are kept.
                - "<>" or "outside" = Rows with values outside of the range
                    "value1" and "value2" (included) are kept.
                - "!" or "not" = Rows that are NOT "value1" are kept.
                - "=" or "equal" = Rows that are "value1" are kept.
            (If unspecified, defaults to extract "UID", "Chrom", "Position",
            "Ref Base", "Alt Base", "Hugo", "Protein Change", "P-value",
            "Score", "Tum-Spec P-value", "Tum-Spec Score", filtering no columns.)
        rename_hugo : Do HUGO symbols need to be renamed?
            If False, doesn't rename HUGO symbols.
            If True, renames equal HUGO symbols: For example, `CROCC, CROCC` will
            become `CROCC_0, CROCC_1`, etc.
            If a string containing a column name is specified,
            renames HUGO symbols with symbol_(new col). For example, CROCC_Y240T
            Appends the renamed HUGOs as a new column named "Renamed HUGOs".
            (Defaults to False.)
        adjust_pval : Whether to adjust p-values. Appends adjusted p-values to
            a new column named "FDR" and "Tum-Spec FDR". Uses BH corrections.
            (Defaults to True)
        drop_na_cols : Specify column(s) to purge from NAs as a list of
            strings with the names of the columns. If False, will not drop any lines.
            Warning: Calculating FDRs without filtering p-values will probably
            lead to NaN FDRs!
            (Defaults to filtering only p-values, passing the string `"default"`)

    Returns:
        The extracted dataframe and the number of lines dropped by filtering.
        This is useful as you can keep track of how many mutations were removed
        as they did not meet some parameter, like passenger mutations.

    Raises:
        KeyError: if the key specified for renaming HUGOs is invalid.
        EmptyFileError: if the file is empty after dropping NAs (see drop_na_cols
            argument).
    """
    if data_id:
        log.info(f"Processing {data_id}")
    else:
        log.info("Processing data with no identifier")

    if select == "default":
        select = {
            "UID": None,
            "Chrom": None,
            "Position": None,
            "Ref Base": None,
            "Alt Base": None,
            "Hugo": None,
            "Protein Change": None,
            "P-value": None,
            "Score": None,
            "Tum-Spec P-value": None,
            "Tum-Spec Score": None,
        }

    if drop_na_cols == "default":
        if "P-value" in select.keys() and "Tum-Spec P-value" in select.keys():
            drop_na_cols = ["P-value", "Tum-Spec P-value"]
        elif "P-value" in select.keys():
            drop_na_cols = ["P-value"]
        elif "Tum-Spec P-value" in select.keys():
            drop_na_cols = ["Tum-Spec P-value"]

    log.debug(
        f"Renaming and selecting columns. Initial column names: {dataframe.columns}"
    )
    try:
        dataframe.columns = [COL_NAMES_REFERENCES[x] for x in dataframe.columns]
    except KeyError as e:
        log.error(
            f"Unmanageable column name found in dataframe ({data_id}).", exc_info=True
        )
        raise

    cols_selected = [x for x in select.keys() if x not in ["Tum-Spec FDR", "FDR"]]
    dataframe = dataframe[cols_selected]
    # ----- Check for missing values ----------------------------------------
    if drop_na_cols is not False:
        log.debug("Checking for missing values.")
        missing = count_missing_rows(dataframe, col=drop_na_cols)
        if missing == 0:
            log.debug("Found no missing values.")
        else:
            dataframe.dropna(inplace=True, subset=drop_na_cols)
            log.debug(
                f" Found {missing} missing value(s). Dropping lines containing them."
            )

        if len(dataframe.index) == 0:
            log.error(f"File {data_id} has no surviving rows after dropping NAs")
            raise EmptyFileError(
                f"File {data_id} has no surviving rows after dropping NAs."
            )

    # ------ Adjust P-values -------------------------------------------------
    if adjust_pval is True:
        log.debug("Adjusting P-values using BH correction.")
        if "P-value" in select.keys():
            p_values = dataframe["P-value"].tolist()
            adj_p_values = multipletests(p_values, method="fdr_bh")
            dataframe["FDR"] = adj_p_values[1]

        if "Tum-Spec P-value" in select.keys():
            p_values = dataframe["Tum-Spec P-value"].tolist()
            adj_p_values = multipletests(p_values, method="fdr_bh")
            dataframe["Tum-Spec FDR"] = adj_p_values[1]

        if "P-value" not in select.keys() and "Tum-Spec P-value" not in select.keys():
            log.warn(
                "P-values need to be adjusted, but none were included."
                " Skipping P-value adjustment."
            )

    # ------ Filtering -------------------------------------------------------
    log.debug("Filtering new dataframe")
    length_in = len(dataframe)
    dataframe = _filter(dataframe, select)
    dropped_mut_nr = length_in - len(dataframe)

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
    return dataframe, dropped_mut_nr


def _filter(input_data: pd.DataFrame, filter_dict: dict) -> pd.DataFrame:
    """Filter target input_data pandas DataFrame with some filter_dict.

    See extract_chasm_plus. Filter takes as filter_dict a dictionary of
    {column_name: filterKey} or a list.
    filterKey is a list like ["operator", value1, value2]. Possible operators:
    - ">" or "more" = Greater than "value1" rows are kept.
    - ">=" or "moreorequal" = Greater than or equal to
      "value1" rows are kept.
    - "<" or "less" = Less than "value1" rows are kept.
    - "<=" or "lessorequal" = Less than or equal to "value1"
      rows are kept.
    - "><" or "between" = Rows with values in between
      "value1" and "value2" (included) are kept.
    - "<>" or "outside" = Rows with values outside of the
      range "value1" and "value2" (included) are kept.
    - "!" or "not" = Rows that are NOT "value1" are kept.
    - "=" or "equal" = Rows that are "value1" are kept.
    If filterKey is a list of only strings, only rows with the string(s)
    are kept.

    Args:
        input_data: Input dataframe to be filtered
        filter_dict: Dictionary of filtering rules

    Returns:
        Filtered dataframe

    Raises:
        IndexError: Not enough numerical parameters are given to a filter rule
        ValueError: An unrecognized filter rule keyword is included
    """

    keywords = (
        ">",
        "more",
        ">=",
        "moreorequal",
        "<",
        "less",
        "<=",
        "lessorequal",
        "<>",
        "between",
        "><",
        "outside",
        "!",
        "not",
        "=",
        "equal",
    )
    for key, value in filter_dict.items():
        if value is None:  # We don't need to check if there is no key
            continue
        # If all items in the list are strings, filter_dict using strings.
        # Cannot use in since it breaks pandas, so I create a truth key and
        # filter_dict using that.
        elif all(isinstance(item, str) for item in value):
            truth = []
            for i in input_data[key]:
                if i in value:
                    truth.append(True)
                else:
                    truth.append(False)
            input_data = input_data[truth]
        # If we only have a list, sieve through it (see _sieve)
        elif value[0] in keywords:
            # The try is only to provide context for the error
            try:
                input_data = _sieve(input_data, key, value)
            except IndexError:
                log.error(
                    f"Not enough numerical parameters given as filter key `{key}`"
                )
                raise IndexError(
                    f"Not enough numerical parameters given as filter key `{key}`"
                )
        else:
            log.error(f"Unrecognized keyword {value[0]}")
            raise ValueError(f"Unrecognized keyword {value[0]}")

    return input_data


def _sieve(input_data: pd.DataFrame, key: str, value: Real) -> pd.DataFrame:
    """Subroutine of _filter. Does the actual filtering of the data.

    See the `_filter` function
    """
    if value[0] == ">" or value[0].lower() == "more":
        input_data = input_data[input_data[key] > value[1]]

    elif value[0] == ">=" or value[0].lower() == "moreorequal":
        input_data = input_data[input_data[key] >= value[1]]

    elif value[0] == "<" or value[0].lower() == "less":
        input_data = input_data[input_data[key] < value[1]]

    elif value[0] == "<=" or value[0].lower() == "lessorequal":
        input_data = input_data[input_data[key] <= value[1]]

    elif value[0] == "><" or value[0].lower() == "between":
        # The atomic & is required to not break pandas
        input_data = input_data[
            (input_data[key] >= min(value[1:])) & (input_data[key] <= max(value[1:]))
        ]

    elif value[0] == "<>" or value[0].lower() == "outside":
        # The atomic | is required to not break pandas
        input_data = input_data[
            (input_data[key] <= min(value[1:])) | (input_data[key] >= max(value[1:]))
        ]

    elif value[0] == "!" or value[0].lower() == "not":
        input_data = input_data[input_data[key] != value[1]]

    elif value[0] == "=" or value[0].lower() == "equal":
        input_data = input_data[input_data[key] == value[1]]

    return input_data


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


def extract_chasm_plus(
    file_in: Union[Path, str], variable: str, filter_tum: bool, value: Real
) -> tuple[pd.DataFrame, int]:
    """Wrapper function for extractor_plus, configuring it for the pipeline.

    Reads FILE_IN, extracts it, and filters it by VARIABLE dynamically.
    Keeps high scores, and low p-values.

    Args:
        file_in : Full path to the file to be extracted.
        variable :  Variable to filter the data by. Can be "pvalue" for CHASM
            P-value, "fdr" for CHASM FDR or "score" for CHASM score. Case-insensitive.
            Keeps low p-values/FDr values and high scores.
        filter_tum: Filter by tumor-specific values?
        value : Value to filter by. Keeps values strictly lower than this for
            p-values and FDR, or values strictly higher than this for scores.

    Returns:
        The output from the configured call to `extractor_plus`

    Raises:
        ValueError if variable is not "pvalue", "score" or "fdr".
    """
    input_dict = {
        "Chrom": None,
        "Position": None,
        "Ref Base": None,
        "Alt Base": None,
        "Hugo": None,
        "Protein Change": None,
        "P-value": None,
        "Score": None,
        "Tum-Spec P-value": None,
        "Tum-Spec Score": None,
        "FDR": None,
        "Tum-Spec FDR": None,
    }
    # Update input dictionaries
    if variable.lower() == "pvalue":
        if filter_tum is False:
            log.debug(
                f"ExtractCHASMPlus starting. Keeping P-values lower than {value}."
            )
            input_dict.update({"P-value": ["<", value]})
        else:
            log.debug(
                f"ExtractCHASMPlus starting. Keeping Tumour P-values lower than {value}."
            )
            input_dict.update({"Tum-Spec P-value": ["<", value]})

    elif variable.lower() == "fdr":
        if filter_tum is False:
            log.debug(f"ExtractCHASM starting. Keeping FDR scores lower than {value}.")
            input_dict.update({"FDR": ["<", value]})
        else:
            log.debug(
                f"ExtractCHASM starting. Keeping Tumour FDR scores lower than {value}."
            )
            input_dict.update({"Tum-Spec FDR": ["<", value]})

    elif variable.lower() == "score":
        if filter_tum is False:
            log.debug(
                f"ExtractCHASM starting. Keeping CHASM scores higher than {value}."
            )
            input_dict.update({"Score": [">", value]})
        else:
            log.debug(
                f"ExtractCHASM starting. Keeping CHASM scores higher than {value}."
            )
            input_dict.update({"Tum-Spec Score": [">", value]})
    else:
        raise ValueError(f"Unrecognized variable {variable}.")

    file_in = Path(file_in)

    return extractor_plus(
        dataframe=input_caster(file_path=file_in),
        data_id=file_in.name,
        select=input_dict,
        rename_hugo="Protein Change",
        adjust_pval=True,
    )


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


class Project:
    """Unified class to manipulate CHASM mutational data.

    Attributes:
        id (str): ID of the project. For example, "TCGA-BRCA"
        patients (dict): Dictionary of key-value pairs where each key is a
            patient identifier and each entry is the associated patient object.
        zero_mut_nr (:obj:`int`, optional): Number of patients without
            mutations as reported by iterextract. Present only if project was
            generated by extracting files.
        zero_mut_ids (list[str], optional): List of patient ids of patients
            with 0 mutations after extraction as reported by iterextract.
            Present only if project was generated by extracting files.
    """

    def __init__(
        self,
        project_name: str,
        extract_path: Union[str, Path],
        extract_options: dict = {"variable": "fdr", "filter_tum": True, "value": 0.05},
        temp_folder: Optional[Union[str, Path]] = None,
        cleanup: bool = False,
        patients_to_retrieve: int = 100_000,
    ):
        """Constructor to generate project objects.

        Generate a project object containing all data regarding both mutations
        (from either a crystal or a set of excel files) and clinical data
        (downloaded from the GDC portal directly).

        It can handle and will include in the clinical csv patients with 0 mutations
        returned from extractCHASMPlus.

        Arguments
        ---------
        project_name: Identifier of the project as reported on the GDC data
            portal. For example, TCGA-BRCA.
        extract_path: Full path to a folder containing mutational data in the
            form of files that must be extracted.
        extract_options: Options given to iterextract for specifying extraction
            method. Supports the keys `variable`, `value` and `filter_tum`
        temp_folder: Full path to temporary folded where intermediate files
            will be stored. Only used if extracting files.
            Defaults to the folder where the script is placed.
        cleanup: Whether or not to remove the temporary folder once the
            project has been generated. Defaults to False
        patients_to_retrieve:  Number of patients to download from the GDC.
            If the actual number of patients is lower than this value, all
            patients are downloaded instead. (This is implemented by the GDC
            API itself). Therefore, just keep it large enough to download all
            patients. Defaults to 100'000.

        Raises:
            ValueError if both or neither of `crystal_path` or `extract_path`
            are specified.
        """
        self.id = project_name
        self.patients = {}
        self.zero_mut_nr = None
        self.zero_mut_ids = None

        log.info(
            c(
                "Starting the generation of a new Project with",
                " the following arguments:\n",
                f"\t\tProject name = {project_name}\n",
                f"\t\tShould I cleanup after myself? {cleanup}\n",
                f"\t\tExtracting from folder {extract_path}\n",
                "\t\textract variable = {}\n".format(extract_options["variable"]),
                "\t\textract value    = {}\n".format(extract_options["value"]),
                "\t\tfilter tumor values? = {}".format(extract_options["filter_tum"]),
            )
        )
        extract_path = Path(extract_path)
        if temp_folder is None:
            log.debug("No temporary folder path was specified. Creating a new one.")
            temp_folder = Path(os.path.realpath(__file__))
            temp_folder = temp_folder.parent
            os.makedirs(name=f".proj_{project_name}_temp", exist_ok=True)
            temp_folder /= f".proj_{project_name}_temp"
            log.debug(f"Temporary folder created @ {temp_folder}")
        else:
            temp_folder = Path(temp_folder)
            temp_folder /= f".proj_{project_name}_temp"
            log.debug(f"Temporary folder set to be @ {temp_folder} by user input.")

        os.makedirs(temp_folder, exist_ok=True)
        os.makedirs(temp_folder / "extracted_files", exist_ok=True)

        # ----- Retrieve clinical data from GDC -----------------------------
        log.debug(f"Retrieving clinical data for {project_name}")
        clinical_data = call_portal(project_name, patients_to_retrieve)
        # ----- Load the crystal --------------------------------------------

        log.info("Beginning to extract files.")
        passenger_freqs = {}
        zero_mut_ids = []
        for path in cqdm(
            extract_path.iterdir(),
            total=len(list(extract_path.iterdir())),
            desc="Extracting files",
        ):
            patient_id = path.name[:12]
            try:
                extracted_frame, passenger_freq = extract_chasm_plus(
                    file_in=path,
                    variable=extract_options["variable"],
                    filter_tum=extract_options["filter_tum"],
                    value=extract_options["value"],
                )
            except EmptyFileError:
                log.warning(
                    f"File {path.name} was empty after filtering NAs. Ignoring it."
                )
                continue

            passenger_freqs[patient_id] = passenger_freq
            extracted_frame.to_csv(
                temp_folder / "extracted_files" / (path.name + ".csv")
            )
            if extracted_frame.empty:
                zero_mut_ids.append(patient_id)

        if len(zero_mut_ids) != 0:
            self.zero_mut_nr = len(zero_mut_ids)
            self.zero_mut_ids = zero_mut_ids

        log.debug(f"Extracted files are in {temp_folder}/extracted_files")
        crystal = fuse_csvs(target_directory=(temp_folder / "extracted_files"))
        # Cleanup the ids from the filenames
        crystal["Identifier"] = [x[:12] for x in crystal["Identifier"]]
        crystal_path = temp_folder / f"{self.id}_crystal.csv"
        crystal.to_csv(crystal_path)
        log.debug(f"Files have been fused. The resulting crystal is in {crystal_path}")
        log.debug(f"Finished extracting")

        # Cleaning the Identifier (remove trailing characters)
        crystal["Identifier"] = crystal["Identifier"].map(lambda x: x[0:12])
        unique_patient_ids = uniques(crystal["Identifier"])
        # Taking into account all patients in the crystal
        if self.zero_mut_ids is not None:
            unique_patient_ids += uniques(self.zero_mut_ids)
        # Putting the clinical data in the patients
        missing_ids = []
        for ident in unique_patient_ids:
            try:
                # Patient class builder just eats the clinical data
                patient = Patient(
                    clinical_data.loc[clinical_data["submitter_id"] == ident]
                )
                patient.id = ident[:12]
                self.patients.update({ident: patient})
            except ValueError:
                missing_ids.append(ident)
        if len(missing_ids) != 0:
            log.warn(
                (
                    f"There are {len(missing_ids)} Patients that are"
                    " present in the crystal that have no"
                    " match in the downloaded data."
                    ' Is "number" large enough?'
                    " Is the project name correct?\n"
                    "Perhaps the GDC doesn't contain the clinical data for"
                    " these patients.\n"
                    " The missing IDs are: {}".format(sorted(missing_ids))
                )
            )
        else:
            log.debug("No patients present in the crystal were left out.")
        # Put the mutations in the patients
        missing_muts = []
        for _, mut_data in crystal.iterrows():
            mutation = Mutation(mut_data)
            try:
                self.patients[mut_data["Identifier"]].mutations.append(mutation)
            except KeyError:
                missing_muts.append(mut_data["Identifier"])
        if len(missing_muts) != 0:
            log.warn(
                (
                    "There are some mutations that I cannot give to"
                    " a patient as they are missing from the array."
                    " It's probably due to the warning you should"
                    " have seen before. The culprits are {}".format(
                        sorted(uniques(missing_muts))
                    )
                )
            )
        # Update all frequencies ----------------------------------------
        for patient in self.patients.values():
            patient.passenger_mut_nr = passenger_freqs[patient.id]
        log.debug("Finished creating project")

        if cleanup:
            log.debug(f"Cleaning up temporary folder @ {temp_folder}")
            rmtree(temp_folder, ignore_errors=True)
        # It's a bit hard to see but __init__ finishes here -----------------

    def generate_clinical_freqs(self, path_out: Union[Path, str]) -> None:
        """Create a file with the clinical information specified + the
        mutation frequency for all patients in the project"""
        log.info(f"Creating a clinical frequency table from the project {self.id}")
        frequency_frame = {"submitter_id", "frequency"}
        frequency_frame = pd.DataFrame(data=frequency_frame)
        clinical_frame = pd.DataFrame()
        for patient in self.patients.values():
            frequency_data = {
                "submitter_id": [patient.id],
                "frequency": [patient.driver_mut_number],
                "passenger_freq": [patient.passenger_mut_nr],
            }
            frequency_data = pd.DataFrame(frequency_data)
            clinical_frame = clinical_frame.append(
                patient.clinical_data, ignore_index=True, sort=True
            )
            frequency_frame = frequency_frame.append(
                frequency_data, ignore_index=True, sort=True
            )
        merged_frame = frequency_frame.merge(clinical_frame, on="submitter_id")
        path_out = Path(path_out)
        merged_frame.to_csv(path_out)
        log.info(
            (
                "Finished creating the clinical frequency table."
                f" It can be found at {path_out}"
            )
        )
        return None

    def save_to_file(self, path_out: Union[Path, str]):
        """Pickle and store in a picklefile the project"""
        log.info(f"Writing the project to a binary file @ {path_out}")
        path_out = Path(path_out)
        with path_out.open("wb+") as output_stream:
            pickle.dump(self, output_stream)

    def generate_lymphocyte_infiltrate(self, lymph_data, path_out):
        """Write frequencies and lymphocyte infiltrate data."""
        log.info(f"Writing lymphocyte infiltrate frequencies @ {path_out}")
        log.info(f"Database of lympocyte infilitrate is @ {lymph_data}")
        patients = []
        freqs = []
        passenger_freqs = []
        for patient in self.patients.values():
            patients.append(patient.id)
            freqs.append(patient.driver_mut_number)
            passenger_freqs.append(patient.passenger_mut_nr)
        project_patients = pd.DataFrame(
            data={
                "ParticipantBarcode": patients,
                "frequency": freqs,
                "passenger_freq": passenger_freqs,
            }
        )
        database = pd.read_excel(Path(lymph_data))
        merged = project_patients.merge(database, on="ParticipantBarcode")
        merged.to_csv(Path(path_out))


class Patient:
    def __init__(self, clinical_data):
        self.id = None
        self.mutations = []
        self.clinical_data = clinical_data
        self.passenger_mut_nr = None

    @property
    def driver_mut_number(self):
        return len(self.mutations)


class Mutation:
    def __init__(self, crystal_line):
        self.chromosome = crystal_line["Chrom"]
        self.position = crystal_line["Position"]
        self.ref_base = crystal_line["Ref Base"]
        self.alt_base = crystal_line["Alt Base"]
        self.hugo_symbol = crystal_line["Hugo"]
        self.protein_change = crystal_line["Protein Change"]
        self.pvalue = crystal_line["P-value"]
        self.ts_pvalue = crystal_line["Tum-Spec P-value"]
        self.score = crystal_line["Score"]
        self.fdr = crystal_line["FDR"]
        self.ts_fdr = crystal_line["Tum-Spec FDR"]


def unpickle(path: Union[Path, str], override_response: bool = False) -> Any:
    """Load a pickled file

    WARNING!! Pickled files may pose a security threat.
    Never un pickle unknown pickle files downloaded from the internet.

    Args:
        path: Full path to the pickled file.
        override_response: Whether to override questions regarding safety
            and unpickle anyway. Defaults to False

    Returns:
        Unpickled object
    """
    if override_response is False:
        log.warn(
            "Pickled files may pose a security threat. Never unpickle"
            " unknown pickle files downloaded from the internet."
        )
        response = input("Do you wish to continue loading? [Y/N]: ")
    else:
        response = "yes"

    if response.lower().startswith("y"):
        path = Path(path)
        with path.open("rb") as file:
            obj = pickle.load(file)
        return obj
    else:
        return False


@cli.command(name="project")
@click.argument("target-path")
@click.argument("output-path")
@click.option(
    "--extract-variable",
    default="fdr",
    help="Specify extraction variable. Defaults to FDR. Can be pvalue, score or fdr",
)
@click.option(
    "--extract-value",
    default=0.05,
    type=click.FLOAT,
    help="Value to extract with. Keeps low pvalues/fdrs and high scores. Defaults to 0.05",
)
@click.option(
    "--filter-pancancer",
    default=False,
    is_flag=True,
    help="Force filtering on pancancer predictions. By default uses tumor specific ones.",
)
@click.option(
    "--keep-loose",
    default=False,
    is_flag=True,
    help="Force keeping intermediate files. Deletes them by default.",
)
@click.option(
    "--save-pickle",
    is_flag=True,
    default=False,
    help="Save as a pickle object the created project.",
)
@click.option(
    "--patients-to-download",
    default=100_000,
    type=click.INT,
    help="Specify patients to retrieve from the GDC data portal. Defaults to 100'000.",
)
def project_command(
    target_path: Union[Path, str],
    output_path: Union[Path, str],
    extract_variable: str,
    extract_value: float,
    filter_pancancer: bool,
    keep_loose: bool,
    save_pickle: bool,
    patients_to_download: int,
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
        project = Project(
            project_name=project_id,
            extract_path=path,
            extract_options={
                "variable": extract_variable,
                "filter_tum": not filter_pancancer,
                "value": extract_value,
            },
            temp_folder=project_output_path,
            cleanup=not keep_loose,
            patients_to_retrieve=patients_to_download,
        )

        project.generate_clinical_freqs(
            path_out=project_output_path / f"{project_id}_clinical.csv"
        )
        if save_pickle:
            project.save_to_file(project_output_path / f"{project_id}_project.pickle")
