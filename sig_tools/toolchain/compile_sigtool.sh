#!/bin/bash
set -euxo pipefail
self_path=$(dirname $0)
sig_tool=${1?"sig_tool not specified"}
tool_dir=${self_path}/../tools

sig_tool_path="$tool_dir/${sig_tool}.m"
[[ -f $sig_tool_path ]] || echo "$sig_tool_path not found"

#***Change to point to the correct path(s) ***
export MSCRIPT_PATH=${self_path}/..

# Path to Mortar library
export MORTARPATH=$MSCRIPT_PATH
# Update Mortar before compilation
export UPDATEMORTAR=-1
# use version control for binary files
export USE_SUBVERSION=-1

export MATLAB_ROOT=/opt/matlab_2018b
#export MCR_ROOT=/opt/mcr/v95
export MCR_VERSION=v95

cd $tool_dir
../toolchain/makemcc $sig_tool
