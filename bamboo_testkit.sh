#!/bin/bash
set -vex

export PATH=/pbi/dept/secondary/builds/develop/current_smrttools-release_installdir/private/otherbins/all/bin:$PATH

ls -larth
cd testkit/hgap5_fake_synth5k
ls -larth
make test-local
