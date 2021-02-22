"""Ask the GDC data portal for clinical information. Unpack it (clean it) from
the JSON and paste it inside a Pandas Dataframe together with the standard
patient identifier.
"""

import json
import logging
import time
from functools import reduce
from pathlib import Path

import click
import pandas as pd
import requests

from edmund.entrypoint import cli

log = logging.getLogger(__name__)


def call_portal(project_id: str, number: int = 1_000_000):
    """Retrieves clinical data regarding a project from the GDC database

    Requires an internet connection. Downloads data regarding all patients in
    the project (up to 'number'), and will then merge it into a single pandas
    DataFrame for manipulation.

    Writes all information in a log file in the same folder as the script.

    Args:
        project_id : The project's ID, such as TCGA-BRCA
        number: Download the first "number" of patients, up to all patients.
            This should generally be kept as a large enough values. Defaults
            to 1_000_000

    Returns:
        Pandas dataframe containing clinical information
    """

    cases_endpt = "https://api.gdc.cancer.gov/cases/"
    filters = {
        "op": "in",
        "content": {"field": "project.project_id", "value": [project_id]},
    }
    data_types = ["diagnoses", "demographic", "exposures"]

    dataframes = []

    for data_type in data_types:
        params = {
            "filters": json.dumps(filters),
            "format": "JSON",
            "expand": data_type,
            "size": str(number),
        }
        log.debug(f"Calling GDC portal for {data_type} data")
        start = time.perf_counter()
        response = requests.get(cases_endpt, params=params)
        log.debug(
            f"Received {len(response.content)} bytes in {(time.perf_counter() - start):.2f} seconds"
        )
        warnings = json.loads(response.content.decode("utf-8"))["warnings"]
        if warnings:
            log.warning(f"There were some warnings when downloading data: {warnings}")
        # Decode and unpack the JSON response
        decoded_response = json.loads(response.content.decode("utf-8"))["data"]["hits"]
        # Clean up the data, and put it in a dataframe ------------------------
        missing_diagnoses = []
        cases = []
        for patient in decoded_response:
            clean_data = {}
            try:
                if data_type != "demographic":
                    # Diagnoses and exposures are in a list of 1 element, so
                    # I'm unlisting them here (the [0])
                    clean_data.update(patient[data_type][0])
                else:
                    # Demographic is just a dictionary, no need to unlist
                    clean_data.update(patient[data_type])
            except KeyError:
                missing_diagnoses.append(patient["submitter_id"])
            # Add the relevant patient ID to the cleaned data for merging
            clean_data.update({"submitter_id": patient["submitter_id"]})
            cases.append(clean_data)
        # Warn the user if something went wrong when retrieving the data
        if missing_diagnoses:
            str_missing_diagnoses = ", ".join(missing_diagnoses)
            log.warning(
                f"I found one or more missing {data_type}: {str_missing_diagnoses}"
            )
        # Finally, add the dataframe to the dataframe list
        dataframes.append(pd.DataFrame(cases))
    # Collapse all dataframes into a single one
    log.debug("Collapsing received data")
    merged_frame = reduce(
        lambda x, y: pd.merge(x, y, on="submitter_id", how="outer"), dataframes
    )
    return merged_frame


@cli.command()
@click.argument("project_id", type=str)
@click.argument("output_file", type=str)
@click.argument("number", type=int, default=1_000_000)
def get_clinical_data(project_id: str, output_file: str, number):
    """Retrieves clinical data from the GDC data portal given a TCGA ID

    Gets data from the first NUMBER patients (default to 1 Million) in the
    TCGA project with id PROJECT_ID, and save them in csv format to OUTPUT_FILE.

    The missing values are replaced with the string "not reported",
    like TCGA does with their missing variables.
    """
    # TODO: Possible that they mean "NA" as one thing, and "not reported"
    # as another, but who cares.
    with Path(output_file).open("w+") as file:
        call_portal(project_id, number).fillna("not reported").to_csv(file)
