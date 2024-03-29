## Set Up Environment (to run in WSL)

### Install conda

Run the following code in the Ubuntu terminal:

```bash
# Download the installer
# Replace <INSTALLER_VERSION> with the version you would like
# For this project, I used 2024.02-1
curl -O https://repo.anaconda.com/archive/Anaconda3-<INSTALLER_VERSION>-Linux-x86_64.sh

# Install
bash Anaconda3-<INSTALLER_VERSION>-Linux-x86_64.sh
```

### Create Environment

After installing conda, create the environment by running the following code in the terminal:

```bash
# Create the environment
# Replace <ENV_NAME> with a name of your choice
conda create --name <ENV_NAME>

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
Run:

```bash
sudo apt install jupyter-core
```

### Open Jupyter Notebook

To create a directory for the project and open Visual Studio Code for writing code, run the following code in the terminal:

```bash
# Make directory
# Replace <DIRECTORY_NAME> with a name of your choice
mkdir <DIRECTORY_NAME>

# Open Visual Studio Code
code
```

## Set Up Solutions for Solving the Heat Equation with the Finite Element Method

