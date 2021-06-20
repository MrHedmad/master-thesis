import logging
import os
import subprocess
from pathlib import Path
from typing import IO, Optional, Union

import click
import numpy as np
import pandas as pd
from click.decorators import command

from edmund.entrypoint import cli
from edmund.utils import cqdm

log = logging.getLogger(__name__)


def read_vcf_header(file: IO):
    """Read a vcf header. Resets pointer to 0"""
    out = []
    for line in file:
        if line.startswith("#"):
            out.append(line)
        else:
            break

    file.seek(0)
    return out


def read_vcf_mutations(file: IO):
    """Returns all mutations in a VCF. Returns pointer to 0"""
    out = []
    for line in file:
        if line.startswith("#"):
            continue
        out.append(line)
    file.seek(0)
    return out


def make_sizes(min_size, max_size, step):
    assert min_size <= max_size
    out = [min_size]

    temp = list(range((max_size - min_size) // step))
    temp = map(lambda x: x + 1, temp)
    temp = map(lambda x: x * step + min_size, temp)
    out.extend(temp)
    return out


def convert_time(timestring):
    x = timestring.count(":")
    if x == 0:
        return float(timestring)
    elif x == 1:
        minutes, sec = timestring.split(":")
        return int(minutes) * 60 + float(sec)
    elif x == 2:
        hours, minutes, sec = timestring.split(":")
        return int(hours) * 3600 + int(minutes) * 60 + float(sec)
    else:
        raise ValueError(f"Cannot convert time {timestring}, too many `:`")


def read_time_output(outstr):
    outstr = outstr.decode("utf-8")
    # The -1 is there as there is a final newline that I need to remove
    # or the split with : fails.
    outstr = outstr.split("\n")[:-1]
    outstr = map(lambda x: x.split(": ")[1], outstr)
    outstr = list(outstr)
    return {
        "kernel_time": float(outstr[2]),
        "user_time": float(outstr[1]),
        "wall_time": convert_time(outstr[4]),
        "max_resident_size_bytes": float(outstr[9]) * 1000,
    }


@cli.command(name="testtool")
@click.argument("call")
@click.argument("input")
@click.argument("output")
@click.argument("scratch")
@click.option(
    "--minsize", "-i", default=100, type=int, help="Minimum mutations in input VCF file"
)
@click.option(
    "--maxsize",
    "-a",
    default=10_000,
    type=int,
    help="Maximum mutations in input VCF file",
)
@click.option(
    "--stepsize",
    "-s",
    default=1000,
    type=int,
    help="Step size, in mutations, of input VCF file",
)
@click.option(
    "--nrsamples",
    "-n",
    default=10,
    type=int,
    help="Times the CALL is repeated before moving on to next input size",
)
@click.option(
    "--env",
    "-e",
    default=None,
    type=str,
    help="Path to custom python environment to run CALL inside",
)
def test_call(
    call: str,
    scratch: Union[str, Path],
    input: Union[str, Path],
    output: Union[str, Path],
    minsize: int,
    maxsize: int,
    stepsize: int,
    nrsamples: int,
    env: Optional[str],
):
    """Test a certain CALL, monitoring it with `time`

    The CALL is repeated multiple times by replacing the input file with
    larger and larger inputs, recording the time taken and memory used.
    CALL is a string with the exact call of the tool, with the input file
    path replaced with `{input}`. That string is replaced with the various
    input file paths.

    Input files are VCFs with more and more mutations. The full input is
    specified in the argument INPUT, and must contain more or equal mutations
    than MAXSIZE.

    At each run, the mutations in the input file are increased by STEPSIZE,
    up to or equal to MAXSIZE. The generated VCFs are saved in SCRATCH and
    deleted at the end of the run.
    Runs of a certain size are repeated SAMPLES time.

    The output is a CSV, saved in the OUTPUT path.

    `time` is called from `/bin/time`.
    """
    log.info("Starting testtool")
    scratch = Path(scratch)
    file_in = Path(input)
    output_path = Path(output)
    env = Path(env) if env else os.environ["PATH"]

    os.makedirs(scratch, exist_ok=True)
    scratch_file = scratch / "temp.vcf"

    realcall = call.format(input=scratch_file.resolve())

    log.info("Counting mutations")
    with file_in.open("r") as file:
        mutations = read_vcf_mutations(file)
        if (tot_muts := len(mutations)) < maxsize:
            raise ValueError(
                f"{maxsize} too large, as input file contains only {tot_muts} mutations"
            )

        header = read_vcf_header(file)

    steps = make_sizes(minsize, maxsize, stepsize)
    # Make empty frame with output
    frame = pd.DataFrame(
        columns=[
            "nr_muts",
            "filesize_bytes",
            "kernel_time",
            "user_time",
            "wall_time",
            "max_resident_size_bytes",
        ]
    )

    for sample_size in cqdm(steps, desc="Testing steps"):
        log.info(f"Sampling {sample_size} mutations")
        samples = np.random.choice(mutations, (nrsamples, sample_size), replace=True)
        for sample in cqdm(samples, desc="Running calls"):
            # Here, we do the single run of the tool
            # Write the temp file with the data
            with scratch_file.open("w+") as scratchfile:
                scratchfile.writelines(header)
                scratchfile.writelines(sample)

            fsize = os.path.getsize(scratch_file)
            # Run the call
            timedcall = ["/bin/time", "-v"]
            timedcall.extend(realcall.split(" "))
            output = subprocess.run(
                timedcall, capture_output=True, env={"PATH": f"{env}/bin:$PATH"}
            )
            log.debug(
                f"Call Finished. Output: {output.stdout}, Stderr: {output.stderr}"
            )
            # Write output to dictionary
            outputdict = read_time_output(output.stderr)
            outputdict.update({"nr_muts": sample_size, "filesize_bytes": fsize})
            frame = frame.append(outputdict, ignore_index=True)
            log.info("Writing partial results after call")
            frame.to_csv(output_path)
        log.info("Finished calls for sample size. Writing (partial) output")
        frame.to_csv(output_path)

    log.info("Finished")

    return None

    # Calculate how many runs we have to do
