#!/bin/bash -xe
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load git/2.8.3
module load gcc/4.9.2
module load ccache/3.2.3

mkdir -p prebuilts/libbzip2-1.0.6
curl -sL http://ossnexus/repository/unsupported/pitchfork/gcc-4.9.2/libbzip2-1.0.6.tgz \
| tar zxf - -C prebuilts/libbzip2-1.0.6
if [ ! -e .distfiles/gtest/release-1.7.0.tar.gz ]; then
  mkdir -p .distfiles/gtest
  curl -sL http://ossnexus/repository/unsupported/distfiles/googletest/release-1.7.0.tar.gz \
    -o .distfiles/gtest/release-1.7.0.tar.gz
fi

cat > pitchfork/settings.mk << EOF
CCACHE_DIR            = /mnt/secondary/Share/tmp/bamboo.mobs.ccachedir
DISTFILES             = $PWD/.distfiles
# from Herb
HAVE_OPENSSL      = /mnt/software/o/openssl/1.0.2a
HAVE_PYTHON       = /mnt/software/p/python/2.7.9/bin/python
HAVE_BOOST        = /mnt/software/b/boost/1.58.0
HAVE_ZLIB         = /mnt/software/z/zlib/1.2.8
HAVE_SAMTOOLS     = /mnt/software/s/samtools/1.3.1mobs
HAVE_NCURSES      = /mnt/software/n/ncurses/5.9
# from MJ
HAVE_HDF5         = /mnt/software/a/anaconda2/4.2.0
HAVE_OPENBLAS     = /mnt/software/o/openblas/0.2.14
HAVE_CMAKE        = /mnt/software/c/cmake/3.2.2/bin/cmake
HAVE_LIBBZIP2     = $PWD/prebuilts/libbzip2-1.0.6
#
htslib_REPO           = $PWD/htslib
pbbam_REPO            = $PWD/pbbam
blasr_REPO            = $PWD/blasr
blasr_libcpp_REPO     = $PWD/blasr_libcpp
pbcore_REPO           = $PWD/pbcore
pbalign_REPO          = $PWD/pbalign
pbcopper_REPO         = $PWD/pbcopper
seqan_REPO            = $PWD/seqan
unanimity_REPO        = $PWD/unanimity
ConsensusCore_REPO    = $PWD/ConsensusCore
GenomicConsensus_REPO = $PWD/GenomicConsensus
pbcommand_REPO        = $PWD/pbcommand
pbcoretools_REPO      = $PWD/pbcoretools
pypeFLOW_REPO         = $PWD/pypeFLOW
daligner_REPO         = $PWD/daligner
dextractor_REPO       = $PWD/dextractor
dmasker_REPO          = $PWD/dmasker
dazzdb_REPO           = $PWD/dazzdb
pbdagcon_REPO         = $PWD/pbdagcon
bam2fastx_REPO        = $PWD/bam2fastx
falcon_polish_REPO    = $PWD/falcon_polish
falcon_kit_REPO       = $PWD/falcon
EOF
cd pitchfork
echo y|make _startover
make -j10 falcon_polish nose
