# This script will setup a working conda environment given the arguements
#
# Expected args:
#   1. WORKSPACE: path to the workspace
#   2. EXPECTED_CONDA_PATH: path to the conda executable
#   3. CONDA_ENV_NAME: name of the conda environment

$workspace=$Args[0]
$expected_conda_path=$Args[1]
$conda_env_name=$Args[2]

Start-Process -NoNewWindow -Wait "$expected_conda_path" -ArgumentList "env remove -n $conda_env_name"

Start-Process -NoNewWindow -Wait "$expected_conda_path" -ArgumentList "env create -f $workspace/mantid-developer-win.yml"
