#!/bin/sh
set -eu
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
./Herwig read -L ../lib -L .libs --repo HerwigDefaults.rpo "$src_dir/testDecayPolicy.in"
