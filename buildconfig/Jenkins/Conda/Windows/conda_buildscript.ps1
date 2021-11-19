#  This script expects to be in a Windows environment, it will run our CI workflow on a Windows environment depending
#  on the flags and args passed. This script will always compile the code, and run the unit tests if a change was
#  made to something other than a .rst file.
#
#  Script usage:
#  conda-buildscript <path-to-workspace> <cmake-preset-name> [options]
#
#  Example command to run a PR build on ubuntu:
#  conda-buildscript /jenkins/workspace_dir/ linux-ci --enable-systemtests --enable-coverage --enable-docs
#
#  Expected args:
#    1. WORKSPACE: path to the workspace/source code that this should run inside
#    2. CMAKE_PRESET: the CMake preset that should be ran to generate the cmake files for this CI job
#
#  Possible flags:
#    --enable-systemtests: Runs the system tests from being compiled or ran
#    --enable-package: Runs a package being produced
#    --enable-docs: Runs the docs from being built
#    --enable-dev-docs: Runs the developer docs from being built
#    --enable-doctests: Runs the documentation tests
#    --clean-build: Clears the build folder and builds from scratch
#    --clean-external-projects: Clear the external projects from the build folder
#
#  Possible parameters:
#    --extra-cmake-flags: Extra flags to pass directly to cmake, enclose in "", defaults to nothing
#    --build-threads: pass the number of threads that can be used to build with, default is 1 thread per logical core

param (
    [switch] $enable_system_tests=$false,
    [switch] $enable_package=$false,
    [switch] $enable_docs=$false,
    [switch] $enable_doc_tests=$false,
    [switch] $enable_dev_docs=$false,
    [switch] $clean_build = $false,
    [switch] $clean_build_external_projects=$false,
    [string] $extra_cmake_flags=""
)

$workspace=$Args[0]
$cmake_preset=$Args[1]
$expected_mambaforge_path="$workspace\mambaforge"
$expected_conda_path="$expected_mambaforge_path\condabin\conda.bat"
$conda_env_name="mantid-developer"
$build_threads=
$build_dir="$workspace\build"
$enable_build_code=$true
$enable_unit_tests=$true

# Setup Mambaforge
. $PSScriptRoot/download_and_install_mambaforge.ps1 $expected_mambaforge_path $expected_conda_path $clean_build

# Setup Conda environment
. $PSScriptRoot/setup_conda_env.ps1 $workspace $expected_conda_path $conda_env_name

# Activate Conda environment
. $expected_conda_path activate $conda_env_name

# Clean the source tree to remove stale configured files ignoring the build directory
. git clean -d -x --force --exclude=build --exclude=mambaforge --exclude=".Xauthority-*"

# Clean up build folder
if ($clean_build -And (Test-Path -Path $build_dir)) {
    Remove-Item -Path $build_dir -Recurse
}

New-Item -Path $workspace -Name "build" -ItemType "directory"

# Clean up items that will cause issues in build
if (Test-Path -Path $build_dir/bin) {
    Remove-Item -Path $build_dir/bin -Recurse
}
if (Test-Path -Path $build_dir/ExternalData) {
    Remove-Item -Path $build_dir/ExternalData -Recurse
}
if (Test-Path -Path $build_dir/Testing) {
    Remove-Item -Path $build_dir/Testing -Recurse
}
Get-ChildItem $workspace -Include TEST-*.xml -Recurse | Remove-Item
if ($clean_build_external_projects -And (Test-Path -Path $build_dir/eigen-*)) {
    Remove-Item -Path $build_dir/eigen-* -Recurse
}
if ($clean_build_external_projects -And (Test-Path -Path $build_dir/eigen-*)) {
    Remove-Item -Path $build_dir/eigen-* -Recurse
}

# TODO implement check_for_changes.sh for dev and user docs only
# If only docs changes:
#    $enable_build_code=$false
#    $enable_unit_tests=$false
#    $enable_system_tests=$false

. $PSScriptRoot/build.ps1 $workspace $cmake_preset $enable_docs $enable_dev_docs $enable_build_code $enable_unit_tests $enable_system_tests $extra_cmake_flags

. $PSScriptRoot/run-tests.ps1 $workspace $enable_system_tests $enable_unit_tests $enable_docs $enable_doc_tests $build_threads