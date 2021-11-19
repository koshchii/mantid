# This script will download and install mambaforge if it's not present where expected
#
# Expected args:
#   1. EXPECTED_MAMBAFORGE_PATH: path to where mambaforge should be installed
#   2. EXPECTED_CONDA_PATH: path to the conda executable
#   3. CLEAN_BUILD: whether or not to force mambaforge to be removed before attempting to install it again

$expected_mambaforge_path=$Args[0]
$expected_conda_path=$Args[1]
$clean_build=$Args[2]

$mambaforge_script_name="Mambaforge-Windows-x86_64.exe"
$url="https://github.com/conda-forge/miniforge/releases/latest/download/$mambaforge_script_name"
$download_location="$PSScriptRoot/$mambaforge_script_name"

if ($clean_build -And (Test-Path -Path $expected_mambaforge_path)) {
   Remove-Item -Path $expected_mambaforge_path -Recurse
}

if (!(Test-Path $expected_conda_path)) {
   if (!(Test-Path $download_location)) {
      # Turn off progress bar because it's horribly inefficient
      $ProgressPreference = 'SilentlyContinue'
      Invoke-WebRequest $url -OutFile $download_location
   }
   Start-Process -NoNewWindow -Wait "$PSScriptRoot/execute_mambaforge.bat" -ArgumentList $download_location, $expected_mambaforge_path
   Remove-Item -Path $download_location
}
