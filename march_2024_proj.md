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

### Using Visual Studio Code

To create a directory for the project and open Visual Studio Code for writing code, run the following code in the terminal:

```bash
# Make directory
# Replace <DIRECTORY_NAME> with a name of your choice
mkdir <DIRECTORY_NAME>

# Open Visual Studio Code
code
```

## Set Up Solutions for Solving the Heat Equation with the Finite Element Method

### Install and Run FEniCS
After activating the environment, run:

```bash
# Install
curl -s https://get.fenicsproject.org | bash

# Run
fenicsproject run
```
