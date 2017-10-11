#!/bin/bash
mkdir -p test-reports
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load anaconda/4.x.x
module load htslib/1.3.1
export PYTHONUSERBASE=$PWD/deployment

cd falcon_polish
nosetests -v \
    --verbose --with-xunit \
    --xunit-file=${bamboo_build_working_directory}/test-reports/fcpolish_xunit.xml \
    utest
chmod +w -R .
