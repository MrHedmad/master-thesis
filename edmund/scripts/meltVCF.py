"""VCF files containing mutations are often a mixture of different mutational
types. This script will "melt" the VCF files into several daughter files,
each containing a different type of mutation. The mutation types are
hard-coded, but all non-coded types will be still saved in a "orphan" file.
It can handle several files at a time, and uses multiprocessing for increased
speed.
The produced files are in folders inside the target one.
Ignores files starting with "." and files not ending in .vcf.

Uses tdqm to show a handy progess bar.
Uses pathlib to provide path compatibility between systems.
"""

import logging
import multiprocessing as mp
import os
from pathlib import Path

from tqdm import tqdm

from edmund.utils import absolute_file_paths, count_files, get_line

log = logging.getLogger(__name__)


class StringSorter:
    def __init__(self) -> None:
        """Utility class to sort strings in many different files"""
        self._switches = {}
        self._connections = {}

    def add_switch(self, path: Path, keyword: str, mode: str = "w+"):
        assert mode in ("w", "w+", "a", "a+"), "A switch must have writeable mode"
        assert self._connections == {}, "Cannot add switches with opened connections."
        self._switches.update({keyword: {"path": path, "mode": mode}})

    def remove_switch(self, keyword: str):
        """Remove a registered switch"""
        assert (
            self._connections == {}
        ), "Cannot remove switches with opened connections."
        try:
            del self._switches[keyword]
        except KeyError as e:
            log.error(f"Cannot remove unregistered keyword {keyword}")
            raise

    def __enter__(self):
        for switch, coordinates in self._switches.items():
            self._connections[switch] = coordinates["path"].open(coordinates["mode"])
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for stream in self._connections.values():
            stream.close()
        self._connections = {}

    def sort(self, string: str, keyword: str):
        assert (
            self._connections != {}
        ), "Connections are empty. Are you in a context manager? Did you add switches?"
        self._connections[keyword].write(string)


def melt_vcf(input_file: str, suffix: str = "default", divide_somatic: bool = True):
    """Melt a vcf file in its components.

    Takes in a .vcf file and splits it: if it contains indels, it essentially
    renames that file and puts it in an "indel" folder.
    If it contains SNPs, it renames the file into two different files in two
    different folders, snp_germline and snp_somatic. The file is split
    line-by-line into two different .vcf files containing the header of the
    original file and either only somatic or germline SNPs.

    Ignores files beginning with "." and files not ending in .vcf.
    The supported version of VCF encoding is 4.1

    Parameters
    ----------
    input_file : str
        Full path to target .vcf file.
    suffix : str, optional
        Suffix to give to the extracted files. (By default attempts to find
        the name of the tool used to generate the .vcf file, and use that.)
    divide_somatic : bool, optional
        Whether or not to divide somatic mutations further into their
        components. See _filter_somatic. (Defaults to True)

    Returns
    -------
    True if file was correctly processed.
    False if file was not processed.

    Dependencies
    ------------
    tqdm --- v. 4.31.1
    python - v. 3.7.2
    """

    input_file = Path(input_file)

    # Ignore files starting with . (unix temp files) and not ending in .vcf
    if input_file.name.startswith("."):
        return False
    elif input_file.suffix != ".vcf":
        return False

    if suffix is not None and suffix != "default":
        if not suffix.startswith("."):
            suffix = "." + suffix
    elif suffix == "default":
        # Get tool name from "softwareName" variable
        suffix = get_line(input_file, 10).split(",")[8]
        suffix = suffix.split("<")[1]
        suffix = "." + suffix

    # Detect mutation type by reading the lines:
    with input_file.open("r") as file:
        ins_check = False
        del_check = False
        snp_check = False
        for line in file:
            if line.startswith("#"):
                continue
            elif "VT=INS" in line.split(";"):
                ins_check = True
            elif "VT=DEL" in line.split(";"):
                del_check = True
            elif "VT=SNP" in line.split(";"):
                snp_check = True
        if del_check and ins_check and snp_check is True:
            raise NotImplementedError(
                "The file contains both somatic"
                " mutations and indels, which melt_vcf"
                " cannot yet melt."
            )
        elif del_check or ins_check is True and snp_check is not True:
            file_content = "indel"
        elif snp_check is True and del_check or ins_check is not True:
            file_content = "snp"
        else:
            raise ValueError(
                "I cannot detect the type of mutations contained"
                " in the file. Is the file of the correct"
                " encoding?"
            )

    output_path = input_file.parent

    # Get the 10th file line (9 counting 0) to extract the patient name
    patient_id = get_line(input_file, 9).split("=")[1]

    if "indel" in file_content:
        # No need to split this file. Just rename it in the new folder.
        output_path /= "indel"
        # Make the folder containing the output
        os.makedirs(output_path, exist_ok=True)

        if suffix is not None:
            output_path /= patient_id + ".indel" + suffix + ".vcf"
        else:
            output_path /= patient_id + ".indel.vcf"

        with input_file.open("r") as fileIN, output_path.open("w+") as fileOUT:
            for line in fileIN:
                fileOUT.write(line)

    elif "snp" in file_content:
        # We need to split this file, so check each line and put it into two
        # different files.
        out_somatic = output_path / "snp_somatic"
        os.makedirs(out_somatic, exist_ok=True)
        if suffix is not None:
            out_somatic /= patient_id + ".somaticsnp" + suffix + ".vcf"
        else:
            out_somatic /= patient_id + ".somaticsnp.vcf"

        out_germline = output_path / "snp_germline"
        os.makedirs(out_germline, exist_ok=True)
        if suffix is not None:
            out_germline /= patient_id + ".germlinesnp" + suffix + ".vcf"
        else:
            out_germline /= patient_id + ".germlinesnp.vcf"

        with input_file.open("r") as fileIN, out_somatic.open(
            "w+"
        ) as somOUT, out_germline.open("w+") as germOUT:

            for line in fileIN:
                if line.startswith("#"):
                    # "#" denotes a header line.
                    germOUT.write(line)
                    somOUT.write(line)
                    continue  # Skip to next line

                if "SS=Somatic" in line.split(";"):
                    somOUT.write(line)
                    continue
                if "SS=Germline" in line.split(";"):
                    germOUT.write(line)
                    continue

        if divide_somatic is True:
            # I give the path out
            _filter_somatic(out_somatic)

    return True


def _filter_somatic(file_path):
    """Take a somatic vcf file and divide it into its components.

    PatientID and suffix are inherited from input file, as this is a
    subroutine of melt_vcf.
    Please do not use as a standalone!
    """

    patient_id = file_path.name.split(".")[0]
    # I get the suffix used before without the patient ID
    suffix = file_path.name.split(".")[1:]
    suffix = ".".join(suffix)
    suffix = "." + suffix

    # Define folder paths for the 8 different variables
    this_folder = file_path.parent

    three_utr_path = this_folder / "threeUTR"
    five_utr_path = this_folder / "fiveUTR"
    igr_path = this_folder / "IGR"
    intron_path = this_folder / "intron"
    missense_path = this_folder / "missense_mutation"
    rna_path = this_folder / "rna"
    silent_path = this_folder / "silent"
    nonsense_path = this_folder / "nonsense"
    de_novo_path = this_folder / "deNovoStart"
    splice_site_path = this_folder / "splice_site"
    linc_rna_path = this_folder / "lincRNA"
    # I want to catch orphan lines, just in case
    orphan_path = this_folder / "orphans"

    # Make the folders if not there
    os.makedirs(three_utr_path, exist_ok=True)
    os.makedirs(five_utr_path, exist_ok=True)
    os.makedirs(igr_path, exist_ok=True)
    os.makedirs(intron_path, exist_ok=True)
    os.makedirs(missense_path, exist_ok=True)
    os.makedirs(rna_path, exist_ok=True)
    os.makedirs(silent_path, exist_ok=True)
    os.makedirs(orphan_path, exist_ok=True)
    os.makedirs(nonsense_path, exist_ok=True)
    os.makedirs(de_novo_path, exist_ok=True)
    os.makedirs(splice_site_path, exist_ok=True)
    os.makedirs(linc_rna_path, exist_ok=True)

    three_utr_path /= patient_id + ".3UTR" + suffix
    five_utr_path /= patient_id + ".5UTR" + suffix
    igr_path /= patient_id + ".IGR" + suffix
    intron_path /= patient_id + ".intron" + suffix
    missense_path /= patient_id + ".missense" + suffix
    rna_path /= patient_id + ".RNA" + suffix
    silent_path /= patient_id + ".silent" + suffix
    nonsense_path /= patient_id + ".nonsense" + suffix
    de_novo_path /= patient_id + ".deNovoStart" + suffix
    splice_site_path /= patient_id + ".spliceSite" + suffix
    linc_rna_path /= patient_id + ".lincRNA" + suffix
    # I want to catch orphan lines, just in case
    orphan_path /= patient_id + ".orphan" + suffix

    with file_path.open("r") as fileIN, three_utr_path.open(
        "w+"
    ) as threeUTR, five_utr_path.open("w+") as fiveUTR, igr_path.open(
        "w+"
    ) as igr, intron_path.open(
        "w+"
    ) as intron, missense_path.open(
        "w+"
    ) as missense, rna_path.open(
        "w+"
    ) as rna, silent_path.open(
        "w+"
    ) as silent, orphan_path.open(
        "w+"
    ) as orphan, nonsense_path.open(
        "w+"
    ) as nonsense, de_novo_path.open(
        "w+"
    ) as de_novo, splice_site_path.open(
        "w+"
    ) as splice_site, linc_rna_path.open(
        "w+"
    ) as lincRNA:

        for line in fileIN:
            if line.startswith("#"):
                # "#" denotes a header line.
                threeUTR.write(line)
                fiveUTR.write(line)
                igr.write(line)
                intron.write(line)
                missense.write(line)
                rna.write(line)
                silent.write(line)
                orphan.write(line)
                nonsense.write(line)
                de_novo.write(line)
                splice_site.write(line)
                lincRNA.write(line)
                continue  # Skip to next line

            if "VC=3'UTR" in line.split(";"):
                threeUTR.write(line)
                continue
            elif "VC=5'UTR" in line.split(";"):
                fiveUTR.write(line)
                continue
            elif "VC=IGR" in line.split(";"):
                igr.write(line)
                continue
            elif "VC=Intron" in line.split(";"):
                intron.write(line)
                continue
            elif "VC=Missense_Mutation" in line.split(";"):
                missense.write(line)
                continue
            elif "VC=RNA" in line.split(";"):
                rna.write(line)
                continue
            elif "VC=Silent" in line.split(";"):
                silent.write(line)
                continue
            elif "VC=lincRNA" in line.split(";"):
                lincRNA.write(line)
                continue
            elif "VC=Nonsense_Mutation" in line.split(";"):
                nonsense.write(line)
                continue
            elif "VC=De_novo_Start_OutOfFrame" in line.split(";"):
                de_novo.write(line)
                continue
            elif "VC=Splice_Site" in line.split(";"):
                splice_site.write(line)
                continue
            else:
                orphan.write(line)
                continue

    return True


def iterative_melt(folder, cores=None):
    """Wrapper for melt_vcf, allowing a pass in all files in a folder.

    Allows the looping over every file in the specified path.
    Supports multiprocessing, by default using all available processors.

    Parameters
    ----------
    folder : str
        Full path to target folder contining .vcf files.
    cores : int, optional
        Number of cores to use. (Defaults to the maximum number found)

    Returns
    -------
    True

    Raises
    ------
    ValueError if multiprocessing fails at dicovering the number of cores
    automatically.
    """
    folder_path = Path(folder)
    iterable = absolute_file_paths(folder_path)
    total_files = count_files(folder_path)

    if cores is None:
        try:
            cores = mp.cpu_count()
        except NotImplementedError:
            raise ValueError(
                (
                    "The number of CPUs couldn't be found"
                    " automatically. Please specify it when running"
                    " the program."
                )
            )

    print(f"Starting to melt VCF files using {cores} core(s).")
    with mp.Pool(cores) as pool:
        list(tqdm(pool.imap(melt_vcf, iterable), total=total_files))

    return True


# Run from commandline
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("folder", help="path to folder")
    parser.add_argument(
        "cores",
        help=("Number of cores to use." " Defaults to maximum available"),
        nargs="?",
        default=None,
        type=int,
    )

    args = parser.parse_args()
    iterative_melt(args.folder, args.cores)
