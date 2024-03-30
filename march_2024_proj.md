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
*Note: This installation method is only available for Ubuntu. For instructions on how to install FEniCS for Mac or Windows, follow the instructions [here](https://fenicsproject.org/pub/tutorial/html/._ftut1003.html#___sec5).*

After activating the environment, run:

```bash
# Add the FEniCS package to your repository and update
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update

# Install FEniCS and install packages for dependenices of FEniCS
sudo apt-get install fenics
sudo apt-get dist-upgrade
```

### Install Jupyter Notebook
In the environment, run in the terminal:

```bash
conda install jupyter
```

### Open Jupyter Notebook
To create a directory for the project and open Jupyter Notebook for writing code, run the following code in the terminal:

```bash
# Make directory
# Replace <DIRECTORY_NAME> with a name of your choice
mkdir <DIRECTORY_NAME>

# Open Jupyter Notebook
jupyter notebook
```

Running `jupyter notebook` creates one HTML file path and two URLs. Copy and paste one of the URLs or HTML file path into a browser to access Jupyter Notebook.


## Set Up Solutions for Solving the Heat Equation with the Finite Element Method

