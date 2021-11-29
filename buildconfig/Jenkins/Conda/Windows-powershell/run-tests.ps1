#!/bin/bash -ex

# This script expects a setup environment with all dependencies availiable in
# the environment
#
# Expected args:
#   1. WORKSPACE: path to the workspace/source code that this should run inside
#   2. ENABLE_SYSTEM_TESTS: whether or not system tests are being ran/built
#   3. ENABLE_UNIT_TESTS: whether or not unit tests are being ran/built
#   4. ENABLE_DOCS: Whether or not the docs have been built
#   5. ENABLE_DOC_TESTS: Whether or not the doc tests should be ran
#   6. BUILD_THREADS: The number of threads to use during testing

$workspace=$Args[0]
$enable_system_tests=$Args[1]
$enable_unit_tests=$Args[2]
$enable_docs=$Args[3]
$enable_doc_tests=$Args[4]
$build_threads=$Args[5]

# Clean up prior to testing
# Prevent race conditions when creating the user config directory
$userconfig_dir="~/.mantid"
if (Test-Path -Path $userconfig_dir ){
    Remove-Item -Path $userconfig_dir -Recurse
}

New-Item -Path "~" -Name ".mantid" -ItemType "directory"
New-Item -Path $userconfig_dir -Name "Mantid.user.properties" -Itemtype "file" -Value "MultiThreaded.MaxCores=2"

if ($enable_unit_tests) {
    . ctest -C Release --no-compress-output -T Test -j$build_threads --schedule-random --output-on-failure
}

if ($enable_docs -and $enable_doc_tests) {
    . cmake --build . --target docs-doctest -j$build_threads
}

if ($enable_system_tests) {
    $userconfig_dir="~/.mantid"
    if (Test-Path -Path $userconfig_dir ){
        Remove-Item -Path $userconfig_dir -Recurse
    }

    New-Item -Path "~" -Name ".mantid" -ItemType "directory"
    New-Item -Path $userconfig_dir -Name Mantid.user.properties -Itemtype "file" -Value "UpdateInstrumentDefinitions.OnStartup = 0"
    Add-Content -Path $userconfig_dir/Mantid.user.properties -Value "usagereports.enabled = 0"
    Add-Content -Path $userconfig_dir/Mantid.user.properties -Value "CheckMantidVersion.OnStartup = 0"

    . $workspace/build/systemtest -j$build_threads
}
