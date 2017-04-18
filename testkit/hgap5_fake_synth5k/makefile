PB_TKR_DEBUG:=--log-level=DEBUG
#P:=/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/private/pacbio/dazzdb/bin:/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/private/pacbio/daligner/bin:/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/private/otherbins/all/bin:/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/private/pacbio/pythonpkgs/falcon/binwrap:/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/private/pacbio/pythonpkgs/falcon/bin:/usr/bin:/bin
#PATH:=${PATH}:$P
#export PATH

try: clean
	pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.hgap_lean --debug -e eid_subread:/pbi/dept/secondary/siv/testdata/hgap/synth5k/synth5k.xml --preset-xml=./preset.xml --force-distributed --output-dir=./job_output
try-old: clean
	pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.polished_falcon_lean --debug -e eid_subread:/pbi/dept/secondary/siv/testdata/hgap/synth5k/synth5k.xml --preset-xml=./preset.xml --force-distributed --output-dir=./job_output
test-local: clean
	pbtestkit-runner ${PB_TKR_DEBUG} testkit.cfg --local-only
test-local-chunk: clean
	pbtestkit-runner ${PB_TKR_DEBUG} testkit.cfg --local-only --force-chunk-mode
test: clean
	pbtestkit-runner ${PB_TKR_DEBUG} testkit.cfg
test-chunk: clean
	pbtestkit-runner ${PB_TKR_DEBUG} testkit.cfg --force-distributed --force-chunk-mode
sanity:
	python -c 'import pbcore; print pbcore'
import:
	pbservice import-dataset --host smrtlink-beta --port 8081 /pbi/dept/secondary/siv/testdata/hgap/synth5k/synth5k.xml
clean:
	rm -rf job_output/ 0.std*
