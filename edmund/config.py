"""This file loads the configuration of the module."""
import json
import logging
import shutil
from copy import deepcopy
from enum import Enum
from pathlib import Path
from typing import Mapping

import click
import yaml

# Don't import the utils script here as it makes a circular import and
# Python gets really really angry
from edmund.entrypoint import cli


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


log = logging.getLogger(__name__)
default_config_path = Path("./edmund/default_config.yml")
user_config_path = Path("./edmund/config.yml")


def load_config():
    try:
        with default_config_path.open("r") as stream:
            default_config = yaml.safe_load(stream)
    except FileNotFoundError as err:
        log.error("Cannot find the default config file (default-config.yml)")
        raise

    try:
        with user_config_path.open("r") as stream:
            log.info("Searching and loading config file")
            config = yaml.safe_load(stream)
    except FileNotFoundError as err:
        log.info("No config file found. Using default parameters")
        config = {}

    return recursive_update(default_config, config)


CONFIG = load_config()


def update_user_config(new_parameters: dict):
    new_config = recursive_update(CONFIG, new_parameters)

    with user_config_path.open("w+") as stream:
        log.info("Dumping updated config")
        yaml.dump(new_config, stream)


class LoggingLevels(Enum):
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    SILENT = logging.CRITICAL + 1


@cli.command()
@click.option(
    "--log-path",
    type=str,
    default=None,
    help="specify logging path. Must not be a folder.",
)
@click.option(
    "--console-verbosity",
    type=str,
    default=None,
    help=(
        "Specify console logging level."
        " Possible values are, from most to least"
        " verbose: 'debug', 'info', 'warning', 'error', 'silent'."
    ),
)
@click.option(
    "--log-verbosity",
    type=str,
    default=None,
    help=(
        "Specify file logging level."
        " Possible values are, from most to least"
        " verbose: 'debug', 'info', 'warning', 'error', 'silent'."
    ),
)
@click.option(
    "--reset",
    is_flag=True,
    help="Reset all options to default. Overrides all other options passed in this call.",
)
@click.option("--show-loading/--no-show-loading", default=None)
def change_config(log_path, console_verbosity, log_verbosity, reset, show_loading):
    """Change the module configuration"""
    # I know, I know, but it works and it's easy. Just this time, ok?
    global CONFIG
    print(show_loading)
    if reset:
        # Reset all options, that is, regenerate the config file.
        log.info("Resetting configuration to defaults")
        shutil.copy(default_config_path, user_config_path)
        CONFIG = load_config()
        return

    possible_levels = ("debug", "info", "warning", "error", "silent")

    if console_verbosity and console_verbosity.lower() not in possible_levels:
        raise ValueError("Invalid console_verbosity value.")
    if log_verbosity and log_verbosity.lower() not in possible_levels:
        raise ValueError("Invalid log_verbosity value.")

    update = deepcopy(CONFIG)
    if log_path:
        if Path(log_path).is_dir():
            raise ValueError("Cannot update log path with this file. Is a folder.")
        update["logs"]["path"] = log_path
    if console_verbosity:
        update["logs"]["stdout_level"] = LoggingLevels[console_verbosity.upper()].value
    if log_verbosity:
        update["logs"]["file_level"] = LoggingLevels[log_verbosity.upper()].value
    if show_loading is None:
        pass
    elif show_loading is True or show_loading is False:
        update["show_loading"] = show_loading

    update_user_config(update)
    CONFIG = load_config()
    log.info(
        "Finished changing configuration. Changes will be applied the next time the module is started."
    )


@cli.command()
def print_config():
    """Print the contents of the current configuration"""
    print(
        "\n".join(
            [
                "Current configuration for Edmund:",
                "Logging:",
                "\tLogging file: {}".format(Path(CONFIG["logs"]["path"]).resolve()),
                "\tConsole logging level: {}".format(
                    LoggingLevels(CONFIG["logs"]["stdout_level"]).name
                ),
                "\tFile logging level: {}".format(
                    LoggingLevels(CONFIG["logs"]["file_level"]).name
                ),
                "Showing loading bars? {}".format(
                    "yes" if CONFIG["show_loading"] else "no"
                ),
            ]
        )
    )
