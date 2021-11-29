:: This script will download and install mambaforge if it's not present where expected
::
:: Expected args:
::   1. EXPECTED_MAMBAFORGE_PATH: path to where mambaforge should be installed
::   2. EXPECTED_CONDA_PATH: path to the conda executable
::   3. CLEAN_BUILD: whether or not to force mambaforge to be removed before attempting to install it again

set EXPECTED_MAMBAFORGE_PATH=%1
set EXPECTED_CONDA_PATH=%2
set CLEAN_BUILD=%3

set MAMBAFORGE_SCRIPT_NAME=Mambaforge-Windows-x86_64.exe
set URL=https://github.com/conda-forge/miniforge/releases/latest/download/%MAMBAFORGE_SCRIPT_NAME%
set DOWNLOAD_LOCATION=%~dp0\%MAMBAFORGE_SCRIPT_NAME%

:: If clean build delete mambaforge installation
IF %CLEAN_BUILD%==true (
    IF EXIST %EXPECTED_MAMBAFORGE_PATH% (
       RMDIR /S /Q %EXPECTED_MAMBAFORGE_PATH%
    )
)

:: Install conda if not here
IF NOT EXIST %EXPECTED_CONDA_PATH% (
    IF NOT EXIST %DOWNLOAD_LOCATION% (
       CALL bitsadmin.exe /transfer "Download Mambaforge" %URL% %DOWNLOAD_LOCATION%
    )
    START /wait "" %DOWNLOAD_LOCATION% /InstallationType=JustMe /RegisterPython=0 /S /D=%EXPECTED_MAMBAFORGE_PATH%
    DEL %DOWNLOAD_LOCATION%
)