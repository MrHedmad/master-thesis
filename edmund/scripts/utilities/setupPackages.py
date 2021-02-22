"""This script install dependencies and can upgrade all Python packages.

It requires Python 3.7 to function.
These packages are used by several scripts in the repository. I don't check
for packages already installed as dependencies, as pip takes care of it.
"""


def main():
    import subprocess
    import sys

    import pkg_resources

    def install(package):
        subprocess.call([sys.executable, "-m", "pip", "install", "--user", package])

    def update(package):
        subprocess.call(
            [sys.executable, "-m", "pip", "install", "--upgrade", "--user", package]
        )

    if sys.version_info[0] < 3 or sys.version_info[1] < 7:
        raise Exception(
            "Script was not run with Python 3.7+.\n"
            "Try typing 'python3.7 setupPackages.py'"
        )

    required_packages = [
        "pandas",
        "numpy",
        "tqdm",
        "scipy",
        "xlwt",
        "xlrd",
        "statsmodels",
    ]

    print(
        "This utility will install all packages required by this repository, and then"
        " ask to update -ALL- packages on this machine.\n"
        "It requires a Python installation of 3.7+.\nIt will use the pip"
        " specific to the python executable used to launch this script, so"
        " it also works for virtual environments.\nThe --user flag will always"
        " be included.\n"
    )

    cont = input("Would you like to install required packages? [Y/N]: ")
    if cont.lower().startswith("y"):
        print("[[[[[ Starting package installation ]]]]]")
        for package in required_packages:
            install(package)
        print("[[[[[ Finished installing required packages ]]]]]")
    else:
        print("Aborted.")

    cont = input(
        "Would you like to update -ALL- installed python packages?"
        " (This may take a while) [Y/N]: "
    )
    if cont.lower().startswith("y"):
        print("[[[[[ Starting package update ]]]]]")
        packages = [dist.project_name for dist in pkg_resources.working_set]
        for package in packages:
            update(package)
        print("[[[[[ Finished updating packages ]]]]]")

    input("Press any key to continue...")


if __name__ == "__main__":
    main()
