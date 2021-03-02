# Master Thesis
Code written by Luca for his Master project. All of these are tested on a machine running a Linux kernel (Windows Subsystems for Linux on Ubuntu).

## Python command line tools
[![Python Version](https://img.shields.io/badge/python-3.9-blue)](https://www.python.org/)

The python command line tools are packaged under the name `Edmund`. The tools can be accessed using the command line from the module level. The scripts can be found in `./edmund`
### Installation
Install Python `3.9` and clone the repo:
```
git clone https://github.com/MrHedmad/master-thesis.git .
```
Make a virtual environment inside the repo folder:
```
cd master-thesis
python3.9 -m venv env
source ./env/bin/activate
```
Install the dependencies:
```
pip install -r requirements.txt
```
The `ed` helper script can be used to access the scripts before a full module installation is done. Try `./ed --help` (changing the permissions might be needed). Or, the full command to enter the interface is `python -m edmund`.

### Usage
The manuals are bundled in the tools themselves. All tools support a `--help` option that prints out the manual. Try `./ed --help` to begin.
### Testing
Some unit tests are available for some functions in the scripts. To run all of them, just run `./runtests` from inside the virtual environment. They are stored inside of `./tests`.
## R scripts
[![R version](https://img.shields.io/badge/R%20version-4.0.2-blue)](https://www.r-project.org/)

The R scripts are written and should probably be run in RStudio. The R environment is currently not version controlled. Install every package that is loaded with `library` and it should be fine.
