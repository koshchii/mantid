:: This script expects a setup environment with all dependencies availiable in
:: the environment
::
:: Expected args:
::   1. WORKSPACE: path to the workspace/source code that this should run inside
::   2. CMAKE_PRESET: the CMake preset that should be ran to generate the cmake files for this CI job
::   3. ENABLE_DOCS: build the user docs
::   4. ENABLE_DEV_DOCS: build the user docs
::   5. ENABLE_BUILD_CODE: whether or not to build the main code target
::   6. ENABLE_UNIT_TESTS: whether or not to build the unit tests target
::   7. ENABLE_SYSTEM_TESTS: whether or not to build the system test data
::   8. EXTRA_CMAKE_FLAGS: Extra flags to pass directly to cmake, enclose in "".
::   9. BUILD_THREADS: Pass the number of threads that can be used to build with

set WORKSPACE=%1
set CMAKE_PRESET=%2
set ENABLE_DOCS=%3
set ENABLE_DEV_DOCS=%4
set ENABLE_BUILD_CODE=%5
set ENABLE_UNIT_TESTS=%6
set ENABLE_SYSTEM_TESTS=%7
set EXTRA_CMAKE_FLAGS=%8
set BUILD_THREADS=%9

:: Set up the location for local object store outside of the build and source
:: tree, which can be shared by multiple builds.
if NOT DEFINED MANTID_DATA_STORE (
  set MANTID_DATA_STORE=%USERPROFILE%\MantidExternalData
)

:: Set CMake generator
CALL %SCRIPT_DIR%\..\cmakegenerator.bat

:: Call CMake
cmake --preset=%CMAKE_PRESET% %MANTID_DATA_STORE% %EXTRA_CMAKE_FLAGS% %WORKSPACE%

cd %WORKSPACE%\build
:: CMake build
if %ENABLE_BUILD_CODE%==true (
    cmake --build . -j%BUILD_THREADS%
)

if %ENABLE_UNIT_TESTS%==true (
    cmake --build . --target AllTests -j%BUILD_THREADS%
)

if %ENABLE_SYSTEM_TESTS%==true (
    cmake --build . --target StandardTestData -j%BUILD_THREADS%
    cmake --build . --target SystemTestData -j%BUILD_THREADS%
)

if %ENABLE_DOCS%==true (
    if %WORKSPACE%\build\docs\doctrees (
       RMDIR /S /Q %WORKSPACE%\build\docs\doctrees
    )
    cmake --build . --target docs-html -j%BUILD_THREADS%
)

if %ENABLE_DEV_DOCS%==true (
    if %WORKSPACE%\build\dev-docs\doctrees (
       RMDIR /S /Q %WORKSPACE%\build\dev-docs\doctrees
    )
    if %WORKSPACE%\build\dev-docs\dev_docs_warnings.txt (
       RM %WORKSPACE%\build\dev-docs\dev_docs_warnings.txt
    )

    cmake --build . --target dev-docs-html -j%BUILD_THREADS%
)