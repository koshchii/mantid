# This script expects a setup environment with all dependencies availiable in
# the environment
#
# Expected args:
#   1. WORKSPACE: path to the workspace/source code that this should run inside
#   2. CMAKE_PRESET: the CMake preset that should be ran to generate the cmake files for this CI job
#   3. ENABLE_DOCS: build the user docs
#   4. ENABLE_DEV_DOCS: build the user docs
#   5. ENABLE_BUILD_CODE: whether or not to build the main code target
#   6. ENABLE_UNIT_TESTS: whether or not to build the unit tests target
#   7. ENABLE_SYSTEM_TESTS: whether or not to build the system test data
#   8. EXTRA_CMAKE_FLAGS: Extra flags to pass directly to cmake, enclose in "".
#   9. BUILD_THREADS: Pass the number of threads that can be used to build with

$workspace=$Args[0]
$cmake_preset=$Args[1]
$enable_docs=$Args[2]
$enable_dev_docs=$Args[3]
$enable_build_code=$Args[4]
$enable_unit_tests=$Args[5]
$enable_system_tests=$Args[6]
$extra_cmake_flags=$Args[7]
$build_threads=$Args[8]
$mantid_data_store_cmake=""

if ($null -eq $env:MANTID_DATA_STORE) {
    $mantid_data_store_cmake = "-DMANTID_DATA_STORE=$env:MANTID_DATA_STORE"
}

. cmake -preset=$cmake_preset $mantid_data_store_cmake $extra_cmake_flags $workspace

cd $workspace/build
if ($enable_build_code) {
    . cmake --build . -j$build_threads
}

if ($enable_unit_tests) {
    . cmake --build . --target AllTests -j$build_threads
}

if ($enable_system_tests) {
    . cmake --build . --target StandardTestData -j$build_threads
    . cmake --build . --target SystemTestData -j$build_threads
}

if ($enable_docs) {
    # Remove doctrees directory so it forces a full reparse. It seems that
    # without this newly added doctests are not executed
    if (Test-Path -Path $workspace/build/docs/doctrees) {
        Remove-Item -Path $workspace/build/docs/doctrees -Recurse
    }
    . cmake --build . --target docs-html -j$build_threads
}

if ($enable_dev_docs) {
    if (Test-Path -Path $workspace/build/dev-docs/doctrees) {
        Remove-Item -Path $workspace/build/dev-docs/doctrees -Recurse
    }
    if (Test-Path -Path $workspace/build/dev-docs/dev_docs_warnings.txt) {
        Remove-Item -Path $workspace/build/dev-docs/dev_docs_warnings.txt
    }
    . cmake --build . --target dev-docs-html -j$build_threads
}
