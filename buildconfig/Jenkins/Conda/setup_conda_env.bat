:: This script will setup a working conda environment given the arguements
::
:: Expected args:
::   1. WORKSPACE: path to the workspace
::   2. EXPECTED_CONDA_PATH: path to the conda executable
::   3. CONDA_ENV_NAME: name of the conda environment

set WORKSPACE=%1
set EXPECTED_CONDA_PATH=%2
set CONDA_ENV_NAME=%3

:: In case environment is already present deactivate so it can be removed.
CALL %EXPECTED_CONDA_PATH% deactivate
CALL %EXPECTED_CONDA_PATH% env remove -n %CONDA_ENV_NAME%

CALL %EXPECTED_CONDA_PATH% env create -f %WORKSPACE%\mantid-developer-win.yml