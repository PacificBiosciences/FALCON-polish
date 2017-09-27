#!/bin/bash
THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load smrttools/incremental
set -vex
cd ${THISDIR}

echo export PATH=/pbi/dept/secondary/builds/develop/current_smrttools-release_installdir/private/otherbins/all/bin:$PATH
echo $PATH
export SMRT_PYTHON_PYTHONPATH_PREPEND=$(pwd)

ls -larth
cd testkit/hgap5_fake_synth5k
ls -larth
make test-local
