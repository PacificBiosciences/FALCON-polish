#!/bin/bash -xe
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load git/2.8.3
module load gcc/4.9.2
module load ccache/3.2.3
export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.mobs.ccachedir
module load boost/1.60
if [[ $BOOST_ROOT =~ /include ]]; then
  set -x
  BOOST_ROOT=$(dirname $BOOST_ROOT)
  set +x
fi
module load anaconda/4.x.x

if [ ! -e .distfiles/gtest/release-1.7.0.tar.gz ]; then
  mkdir -p .distfiles/gtest
  curl -sL http://nexus/repository/unsupported/distfiles/googletest/release-1.7.0.tar.gz \
    -o .distfiles/gtest/release-1.7.0.tar.gz
fi
tar zxf .distfiles/gtest/release-1.7.0.tar.gz
ln -sfn googletest-release-1.7.0 gtest

rm -rf deployment && mkdir deployment
export PYTHONUSERBASE=$PWD/deployment

GenomicConsensus_VERSION=$(cd GenomicConsensus && python setup.py --version 2> /dev/null)
ConsensusCore_VERSION=$(cd ConsensusCore && python setup.py --version 2> /dev/null)
ConsensusCore2_VERSION=$(cd unanimity && python setup.py --version 2> /dev/null)

pip install --user --no-index http://nexus/repository/unsupported/pitchfork/gcc-4.9.2/pythonpkgs/networkx-1.10-py2.py3-none-any.whl
pip install --user --no-index ./falcon
pip install --user --no-index http://nexus/repository/unsupported/pitchfork/gcc-4.9.2/pythonpkgs/avro-1.7.7-cp27-none-any.whl
pip install --user --no-index http://nexus/repository/unsupported/pitchfork/gcc-4.9.2/pythonpkgs/iso8601-0.1.12-py2.py3-none-any.whl
pip install --user --no-index ./pbcommand
pip install --user --no-index http://nexus/repository/unsupported/gcc-4.9.2/pythonpkgs/ConsensusCore-${ConsensusCore_VERSION}-cp27-cp27mu-linux_x86_64.whl
pip install --user --no-index http://nexus/repository/unsupported/gcc-4.9.2/pythonpkgs/ConsensusCore2-${ConsensusCore2_VERSION}-cp27-cp27mu-linux_x86_64.whl
pip install --user --no-index http://nexus/repository/unsupported/pitchfork/gcc-4.9.2/pythonpkgs/pysam-0.9.1.4-cp27-cp27mu-linux_x86_64.whl
pip install --user --no-index ./pbcore
pip install --user --no-index http://nexus/repository/unsupported/gcc-4.9.2/pythonpkgs/GenomicConsensus-${GenomicConsensus_VERSION}-cp27-cp27mu-linux_x86_64.whl
pip install --user --no-index ./falcon_polish



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
