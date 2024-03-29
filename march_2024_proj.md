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
# Replace <ENV_NAME> with a name of your choice
conda create --name <ENV_NAME>
```
