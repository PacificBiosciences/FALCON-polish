#!/bin/bash -e
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load git/2.8.3
module load gcc/6.4.0
module load ccache
#module load python/2.7.13-UCS4
#module load htslib/1.3.1

# The following is missing some stuff.
##export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.mobs.ccachedir
##module load boost/1.60
##if [[ $BOOST_ROOT =~ /include ]]; then
##  set -x
##  BOOST_ROOT=$(dirname $BOOST_ROOT)
##  set +x
##fi
##module load anaconda/4.x.x
##
##if [ ! -e .distfiles/gtest/release-1.7.0.tar.gz ]; then
##  mkdir -p .distfiles/gtest
##  curl -sL http://nexus/repository/unsupported/distfiles/googletest/release-1.7.0.tar.gz \
##    -o .distfiles/gtest/release-1.7.0.tar.gz
##fi
##tar zxf .distfiles/gtest/release-1.7.0.tar.gz
##ln -sfn googletest-release-1.7.0 gtest
##
##rm -rf deployment && mkdir deployment
##export PYTHONUSERBASE=$PWD/deployment
##
##GenomicConsensus_VERSION=$(cd GenomicConsensus && python setup.py --version 2> /dev/null)
##ConsensusCore_VERSION=$(cd ConsensusCore && python setup.py --version 2> /dev/null)
##ConsensusCore2_VERSION=$(cd unanimity && python setup.py --version 2> /dev/null)
##
##pip install --user --no-index http://nexus/repository/unsupported/pitchfork/gcc-4.9.2/pythonpkgs/networkx-1.10-py2.py3-none-any.whl
##pip install --user --no-index http://nexus/repository/unsupported/pitchfork/gcc-4.9.2/pythonpkgs/avro-1.7.7-cp27-none-any.whl
##pip install --user --no-index http://nexus/repository/unsupported/pitchfork/gcc-4.9.2/pythonpkgs/iso8601-0.1.12-py2.py3-none-any.whl
##pip install --user --no-index http://nexus/repository/unsupported/gcc-4.9.2/pythonpkgs/ConsensusCore-${ConsensusCore_VERSION}-cp27-cp27mu-linux_x86_64.whl
##pip install --user --no-index http://nexus/repository/unsupported/gcc-4.9.2/pythonpkgs/ConsensusCore2-${ConsensusCore2_VERSION}-cp27-cp27mu-linux_x86_64.whl
##pip install --user --no-index http://nexus/repository/unsupported/pitchfork/gcc-4.9.2/pythonpkgs/pysam-0.9.1.4-cp27-cp27mu-linux_x86_64.whl
##pip install --user --no-index http://nexus/repository/unsupported/gcc-4.9.2/pythonpkgs/GenomicConsensus-${GenomicConsensus_VERSION}-cp27-cp27mu-linux_x86_64.whl

module load smrttools/incremental
module load python/2.7.13-UCS4

set -vex

which python
which pip

rm -rf $(pwd)/../LOCAL

export PYTHONUSERBASE=$(pwd)/../LOCAL
export PATH=${PYTHONUSERBASE}/bin:${PATH}

# For speed, use cdunn wheelhouse.
WHEELHOUSE=/home/UNIXHOME/cdunn/wheelhouse/gcc-6/
#pip install --user --no-index --find-links=${WHEELHOUSE} pip pytest pytest-cov pylint

cd ..
ppath=$(pwd)/falcon-polish:$(pwd)/falcon:$(pwd)/pypeFLOW:$(pwd)/pbcommand:$(pwd)/pbcore:$(pwd)/pbcoretools
export SMRT_PYTHON_PYTHONPATH_PREPEND=${ppath}
export PYTHONPATH=${ppath}:${PYTHONPATH}
# Really, we only need to build falcon, for the C extensions.
#pip install --user --no-index ./pbcommand ./pbcore ./pbcoretools
#pip install --user --no-index ./pypeFLOW
mods="./pbcommand ./pbcore ./pbcoretools ./pypeFLOW ./FALCON"
pip install --user --find-links=${WHEELHOUSE} ${mods}
cd -
pip install --user .

pip install --user --find-links=${WHEELHOUSE} pytest pytest-cov pylint

python -c 'import falcon_kit; print falcon_kit'
python -c 'import pbcore; print pbcore'
python -c 'import pbcoretools; print pbcoretools'
python -c 'import pysam; print pysam.__version__; print pysam.faidx'

export MY_TEST_FLAGS="-v -s --durations=0 --cov=./falcon_polish/ --cov-report=term-missing --cov-report=xml:coverage.xml --cov-branch"
make pytest
#cat coverage.xml
sed -i -e 's@filename="@filename="./falcon_polish/@g' coverage.xml

make pylint

pwd
ls -larth
find . -name '*.pyc' | xargs rm
#chmod +w -R .
