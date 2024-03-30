## Set Up Environment (to run in WSL)

### Install conda
*Note: These instructions are for installing conda in Linux.*

Run the following code in the Ubuntu terminal:

```bash
# Download the installer
# Replace <INSTALLER_VERSION> with the version you would like (I used 2024.02-1)
curl -O https://repo.anaconda.com/archive/Anaconda3-<INSTALLER_VERSION>-Linux-x86_64.sh

# Install
bash Anaconda3-<INSTALLER_VERSION>-Linux-x86_64.sh
```

### Create Environment
After installing conda, create the environment by running the following code in the terminal:

```bash
# Create the environment
# Replace <ENV_NAME> with a name of your choice (I used "heatequationenv")
# Replace <PYTHON_VERSION> with a version of your choice (I used 3.12.2)
conda create -n <ENV_NAME> python=<PYTHON_VERSION>

# Activate the environment
conda activate <ENV_NAME>
```
The code used here does not specify a python version, so the environment uses the default Python 3.6.9 version.

### Install FEniCS
*Note: This installation method is only available for conda.*

After activating the environment, run:

```bash
 conda install -c conda-forge fenics
```

### Install matplotlib
In the environment, run in the terminal:

```bash
conda install -c conda-forge matplotlib
```

### Install Jupyter Notebook
1. In the environment, run in the terminal:

```bash
conda install jupyter
```

2. Close the terminal.

### Open Jupyter Notebook
To create a directory for the project and open Jupyter Notebook for writing code, open a terminal and run the following code:

```bash
# Make directory
# Replace <DIRECTORY_NAME> with a name of your choice
conda activate <ENV_NAME>
mkdir <DIRECTORY_NAME>

# Open Jupyter Notebook
jupyter notebook
```

Running `jupyter notebook` creates one HTML file path and two URLs. Copy and paste one of the URLs or HTML file path into a browser to access Jupyter Notebook.


## Set Up Solutions for Solving the Heat Equation with the Finite Element Method

