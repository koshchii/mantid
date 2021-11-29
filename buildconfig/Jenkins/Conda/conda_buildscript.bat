@echo on
::  This script expects to be in a Windows environment, it will run our CI workflow on a Windows environment depending
::  on the flags and args passed. This script will always compile the code, and run the unit tests if a change was
::  made to something other than a .rst file.
::
::  Script usage:
::  conda-buildscript <workspace> <cmake-preset-name> [options]
::
::  Example command to run a PR build on ubuntu:
::  conda-buildscript /jenkins/workspace_dir/ linux-ci --enable-systemtests --enable-coverage --enable-docs
::
::  Expected parameters:
::    WORKSPACE: path to the workspace/source code that this should run inside
::    CMAKE_PRESET: the CMake preset that should be ran to generate the cmake files for this CI job
::
::  Possible flags:
::    /enable-systemtests: Runs the system tests from being compiled or ran
::    /enable-package: Runs a package being produced
::    /enable-docs: Runs the docs from being built
::    /enable-dev-docs: Runs the developer docs from being built
::    /enable-doctests: Runs the documentation tests
::    /clean-build: Clears the build folder and builds from scratch
::    /clean-external-projects: Clear the external projects from the build folder
::
::  Possible parameters:
::    --extra-cmake-flags: Extra flags to pass directly to cmake, enclose in "", defaults to nothing
::    --build-threads: pass the number of threads that can be used to build with, default is 1 thread per logical core

:: Check args for issues
set NUMBER_OF_ARGS=0
for %%x in (%*) do Set /A NUMBER_OF_ARGS+=1
set ARGS_OK=true
if %NUMBER_OF_ARGS% LSS 2 (
   set ARGS_OK=false
)
if %ARGS_OK% == false (
    echo "Pass 2 arguements followed by optional flags usage"
    exit /b 1
)

:: Set variables
set WORKSPACE=%1
shift
set CMAKE_PRESET=%1
shift

set EXPECTED_MAMBAFORGE_PATH=%WORKSPACE%\mambaforge
set EXPECTED_CONDA_PATH=%EXPECTED_MAMBAFORGE_PATH%\condabin\mamba.bat
set SCRIPT_DIR=%WORKSPACE%\buildconfig\Jenkins\Conda
set CONDA_ENV_NAME=mantid-developer
set BUILD_THREADS=%NUMBER_OF_PROCESSORS%
set BUILD_DIR=%WORKSPACE%\build
set ENABLE_BUILD_CODE=true
set ENABLE_UNITTEST=true
set ENABLE_SYSTEMTEST=false
set ENABLE_PACKAGE=false
set ENABLE_DOCS=false
set ENABLE_DEV_DOCS=false
set ENABLE_DOCTESTS=false
set CLEAN_BUILD=false
set CLEAN_EXTERNAL_PROJECTS=false

:: Handle the optional flags
:loop
if "%~1" neq "" (
    IF "%1"=="/enable-systemtests" set ENABLE_SYSTEMTEST=true
    IF "%1"=="/enable-package" set ENABLE_PACKAGE=true
    IF "%1"=="/enable-docs" set ENABLE_DOCS=true
    IF "%1"=="/enable-dev-docs" set ENABLE_DEV_DOCS=true
    IF "%1"=="/enable-doctests" set ENABLE_DOCTESTS=true
    IF "%1"=="/clean-build" set CLEAN_BUILD=true
    IF "%1"=="/clean-external-projects" set CLEAN_EXTERNAL_PROJECTS=true
    IF "%1"=="/extra-cmake-flags" (
       set EXTRA_CMAKE_FLAGS=%2
       shift
    )
    IF "%1"=="/build-threads" (
       set BUILD_THREADS=%2
       shift
    )
    shift
    goto :loop
)

echo "Workspace is %WORKSPACE%"
echo "Cmake preset is %CMAKE_PRESET%"
echo "Enable Systemtests is %ENABLE_SYSTEMTEST%"
echo "Enable package is %ENABLE_PACKAGE%"
echo "Enable docs is %ENABLE_DOCS%"
echo "Enable dev docs is %ENABLE_DEV_DOCS%"
echo "Enable doctests is %ENABLE_DOCTESTS%"
echo "Enable clean build is %CLEAN_BUILD%"
echo "Enable clean external projects is %CLEAN_EXTERNAL_PROJECTS%"
echo "Extra cmake flags are %EXTRA_CMAKE_FLAGS%"
echo "Build threads total are %BUILD_THREADS%"

:: Setup mambaforge
CALL %SCRIPT_DIR%\download_and_install_mambaforge.bat %EXPECTED_MAMBAFORGE_PATH% %EXPECTED_CONDA_PATH% %CLEAN_BUILD%

:: Setup Conda environment
CALL %SCRIPT_DIR%\setup_conda_env.bat %WORKSPACE% %EXPECTED_CONDA_PATH% %CONDA_ENV_NAME%

:: Activate Conda environment
CALL %EXPECTED_CONDA_PATH% activate %CONDA_ENV_NAME%

:: Clean the source tree to remove stale configured files ignoring the build directory
git clean -d -x --force --exclude=%BUILD_DIR% --exclude=mambaforge --exclude=".Xauthority-*" --exclude=".vscode"

:: Clean up build folder
IF EXIST %BUILD_DIR%\CMakeCache.txt (
  call rg CMAKE_GENERATOR:INTERNAL %BUILD_DIR%\CMakeCache.txt > %BUILD_DIR%\cmake_generator.log
  call rg -q "%CM_GENERATOR%" %BUILD_DIR%\cmake_generator.log
  if ERRORLEVEL 1 (
    set CLEAN_BUILD=yes
    echo Previous build used a different compiler. Performing a clean build.
  ) else (
    set CLEAN_BUILD=no
    echo Previous build used the same compiler. No need to clean.
  )
)

IF %CLEAN_BUILD%==true (
    IF EXIST %BUILD_DIR% (
          RMDIR /S /Q %BUILD_DIR%
    )
) 

IF EXIST %BUILD_DIR% (
  RMDIR /S /Q %BUILD_DIR%\bin %BUILD_DIR%\ExternalData %BUILD_DIR%\Testing
  PUSHD %BUILD_DIR%
  FOR /f %%F in ('dir /b /a-d /S "TEST-*.xml"') do DEL /Q %%F >/nul
  POPD
  if "!CLEAN_EXTERNAL_PROJECTS!" == "true" (
    RMDIR /S /Q %BUILD_DIR%\eigen-prefix
    RMDIR /S /Q %BUILD_DIR%\googletest-download %BUILD_DIR%\googletest-src
  )
) else (
  MD %BUILD_DIR%
)

del /Q %BUILD_DIR%\*.exe

:: Build

CALL %SCRIPT_DIR%\build.bat %WORKSPACE% %CMAKE_PRESET% %ENABLE_DOCS% %ENABLE_DEV_DOCS% %ENABLE_BUILD_CODE% %ENABLE_UNITTEST% %ENABLE_SYSTEMTEST% %EXTRA_CMAKE_FLAGS% %BUILD_THREADS%

exit /b 1
:: Tests
CALL %SCRIPT_DIR%\run_tests.bat %WORKSPACE% %ENABLE_SYSTEMTEST% %ENABLE_UNITTEST% %ENABLE_DOCS% %ENABLE_DOCTESTS% %BUILD_THREADS%
