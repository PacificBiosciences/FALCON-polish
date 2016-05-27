from __future__ import absolute_import
from .. import sys
from .. import stats_preassembly

from pypeflow.pwatcher_bridge import PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase
from pypeflow.controller import PypeThreadWorkflow
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase

import contextlib
import gzip
import json
import logging
import os
import pprint
import re
import cStringIO

log = logging.getLogger(__name__)

#PypeWorkflow = PypeThreadWorkflow
#PypeTaskBase = PypeThreadTaskBase
PypeWorkflow = PypeProcWatcherWorkflow
PypeTaskBase = MyFakePypeThreadTaskBase

def get_length_cutoff(length_cutoff, fn):
    if length_cutoff < 0:
        try:
            length_cutoff = int(open(fn).read().strip())
            log.info('length_cutoff=%d from %r' %(length_cutoff, fn))
        except Exception:
            log.exception('Unable to read length_cutoff from "%s".' %fn)
    return length_cutoff # possibly updated

def updated_cfg(options_dict):
    opts = dict()
    for key, val in options_dict.iteritems():
        # Drop comments (keys w/ leading ~).
        if key.startswith('~'):
            continue
        # Strip leading and trailing ws, just in case.
        opts[key] = val
    opts['job_type'] = 'local' # OVERRIDE FOR NOW
    def add(key, val):
        if not key in opts:
            opts[key] = val
    add('input_fofn', 'NA') # actually, we do not need this anymore
    add('target', 'assembly')
    add('sge_option_da', 'NA')
    add('sge_option_la', 'NA')
    add('sge_option_pda', 'NA')
    add('sge_option_pla', 'NA')
    add('sge_option_fc', 'NA')
    add('sge_option_cns', 'NA')
    return opts
def dict2ini(ofs, options_dict):
    ofs.write('[General]\n')
    for key, val in sorted(options_dict.items()):
        ofs.write('{} = {}\n'.format(key, val))
def dict2json(ofs, options_dict):
    content = json.dumps(options_dict, sort_keys=True, indent=4, separators=(',', ': '))
    ofs.write(content + '\n')

@contextlib.contextmanager
def ContentUpdater(fn):
    """Write new content only if differs from old.
    """
    if os.path.exists(fn):
        with open(fn) as f:
            old_content = f.read()
    else:
        old_content = None
    new_writer = cStringIO.StringIO()
    yield new_writer
    new_content = new_writer.getvalue()
    if new_content != old_content:
        with open(fn, 'w') as f:
            f.write(new_content)
def run_prepare_falcon(falcon_parameters, i_fasta_fn, fc_cfg_fn, fc_json_config_fn, input_fofn_fn):
    wdir = os.path.dirname(fc_json_config_fn)
    mkdirs(wdir)
    with ContentUpdater(input_fofn_fn) as f:
        f.write('{}\n'.format(i_fasta_fn))
    config_falcon = updated_cfg(dict(falcon_parameters))
    config_falcon['input_fofn'] = input_fofn_fn
    with ContentUpdater(fc_cfg_fn) as f:
        dict2ini(f, config_falcon)
    with ContentUpdater(fc_json_config_fn) as f:
        dict2json(f, config_falcon)

def task_bam2fasta(self):
    """
        {
            "other_filters": "rq >= 0.7",
            "read_length": 0
            OR
            "filters": "rq>=.7, length gte 1000, length &lt;= 50000"
        }
    """
    i_dataset_fn = fn(self.dataset)
    o_fasta_fn = fn(self.fasta)
    f_dataset_fn = 'filtered.subreadset.xml'
    config = self.parameters.get('pbcoretools.tasks.filterdataset', None) # TODO: Drop this.
    if not config:
        config = self.parameters['pbcoretools']
    other_filters = config.get('other_filters', 'rq >= 0.7')
    read_length = config.get('read_length', 0)
    filters = config.get('filters', None)
    if not filters:
        filters = other_filters + ', length gte {:d}'.format(int(read_length))
    wdir, o_fasta_fn = os.path.split(o_fasta_fn)
    bash = """
python -m falcon_polish.mains.run_filterbam {i_dataset_fn} {f_dataset_fn} '{filters}'
python -m falcon_polish.mains.run_bam2fasta {f_dataset_fn} {o_fasta_fn}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_bam2fasta.sh')
    mkdirs(wdir)
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_prepare_falcon(self):
    """Pre-process FALCON cfg.
    This is super-fast, so it can always run locally.
    """
    i_fasta_fn = fn(self.fasta)
    input_fofn_fn = fn(self.input_fofn)
    fc_cfg_fn = fn(self.fc_cfg)
    fc_json_config_fn = fn(self.fc_json_config)
    config_falcon = self.parameters['falcon']
    run_prepare_falcon(config_falcon, i_fasta_fn, fc_cfg_fn, fc_json_config_fn, input_fofn_fn)
def task_falcon(self):
    fc_cfg_fn = fn(self.fc_cfg)
    o_fasta_fn = fn(self.asm_fasta)
    o_preads_fofn_fn = fn(self.preads_fofn)
    wdir, o_fasta_fn = os.path.split(o_fasta_fn)
    odir, o_preads_fofn_fn = os.path.split(o_preads_fofn_fn)
    assert odir == wdir
    o_preads_fofn_fn = os.path.basename(o_preads_fofn_fn)
    bash = """
#rm -f {o_fasta_fn} {o_preads_fofn_fn} # preassembly report depends on this, so we must not pre-delete
#TODO: Let falcon use logging.json?
fc_run1 {fc_cfg_fn}
ln -sf 2-asm-falcon/p_ctg.fa {o_fasta_fn}
ln -sf 1-preads_ovl/input_preads.fofn {o_preads_fofn_fn}
ls -ltr
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_falcon.sh')
    mkdirs(wdir)
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
    # TODO: Optionally run this on local machine.
    # By running distributed, this would cause qsub within qsub. This might
    # be impossible on some systems, so we need to make this configurable.
def task_report_pre_assembly(self):
    # TODO(CD): Bashify this, in case it is slow.
    i_raw_reads_fofn_fn = fn(self.raw_reads_fofn)
    i_preads_fofn_fn = fn(self.preads_fofn)
    i_length_cutoff_fn = fn(self.length_cutoff)
    o_json_fn = fn(self.pre_assembly_report)
    cfg = self.parameters['falcon']
    genome_length = int(cfg.get('genome_size', 0)) # different name in falcon
    length_cutoff = int(cfg['length_cutoff'])
    length_cutoff = get_length_cutoff(length_cutoff, i_length_cutoff_fn)
    kwds = {
        'i_raw_reads_fofn_fn': i_raw_reads_fofn_fn,
        'i_preads_fofn_fn': i_preads_fofn_fn,
        'genome_length': genome_length,
        'length_cutoff': length_cutoff,
    }
    log.info('Report inputs: {}'.format(repr(kwds)))
    report_dict = stats_preassembly.make_dict(**kwds)
    content = json.dumps(report_dict, sort_keys=True, indent=4, separators=(',', ': '))
    open(o_json_fn, 'w').write(content)
def task_fasta2referenceset(self):
    """Copied from pbsmrtpipe/pb_tasks/pacbio.py:run_fasta_to_referenceset()
    """
    input_file_name = fn(self.fasta)
    output_file_name = fn(self.referenceset)
    wdir, output_file_name = os.path.split(output_file_name)
    bash = """
rm -f {output_file_name} {input_file_name}.fai
dataset create --type ReferenceSet --generateIndices {output_file_name} {input_file_name}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_falcon.sh')
    mkdirs(wdir)
    open(bash_fn, 'w').write(bash)
    cmd = 'cd {wdir} && /bin/bash {bash_fn}'.format(**locals())
    sys.system(cmd)
def task_pbalign_scatter(self):
    """This might have problems if run in /tmp.
    """
    reads_fn = fn(self.dataset)
    referenceset_fn = fn(self.referenceset)
    out_json_fn = fn(self.out_json)
    wdir, out_json_fn = os.path.split(out_json_fn)
    mkdirs(wdir)
    bash = r"""
python -m pbcoretools.tasks.scatter_subread_reference -v --max_nchunks=5 \
        {reads_fn} \
        {referenceset_fn} \
        {out_json_fn}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_pbalign_scatter.sh')
    mkdirs(wdir)
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_pbalign_gather(self):
    ds_out_fn = fn(self.ds_out)
    dos = self.inputDataObjs
    dset_fns = [fn(v) for k,v in dos.items() if k.startswith('alignmentset')]
    wdir = os.path.dirname(ds_out_fn)
    mkdirs(wdir)
    ds_fofn_fn = os.path.join(wdir, 'gathered.alignmentsets.fofn')
    open(ds_fofn_fn, 'w').write('\n'.join(dset_fns) + '\n')
    bash = r"""
python -m falcon_polish.mains.run_pbalign_gather \
        {ds_fofn_fn} \
        {ds_out_fn}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_pbalign_gather.sh')
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_pbalign(self):
    """pbalign will eventually call blasr, like this:
 BlasrService: Align reads to references using blasr.
 BlasrService: Call "blasr /pbi/dept/secondary/siv/testdata/SA3-DS/lambda/2372215/0007_tiny/Analysis_Results/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.all.subreadset.xml /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/falcon_ns.tasks.task_falcon2_run_asm-0/file.fasta -out /scratch/tmpTbV4Ec/wLCUdL.bam  -bam  -bestn 10 -minMatch 12  -nproc 16  -minSubreadLength 50 -minAlnLength 50  -minPctSimilarity 70 -minPctAccuracy 70 -hitPolicy randombest  -randomSeed 1  -minPctSimilarity 70.0 "
 FilterService: Filter alignments using samFilter.
 FilterService: Call "rm -f /scratch/tmpTbV4Ec/aM1Mor.bam && ln -s /scratch/tmpTbV4Ec/wLCUdL.bam /scratch/tmpTbV4Ec/aM1Mor.bam"
 BamPostService: Sort and build index for a bam file.
 BamPostService: Call "samtools sort -m 4G /scratch/tmpTbV4Ec/aM1Mor.bam /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/pbalign.tasks.pbalign-0/aligned.subreads.alignmentset"
 BamPostService: Call "samtools index /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/pbalign.tasks.pbalign-0/aligned.subreads.alignmentset.bam /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/pbalign.tasks.pbalign-0/aligned.subreads.alignmentset.bam.bai"
 BamPostService: Call "pbindex /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/pbalign.tasks.pbalign-0/aligned.subreads.alignmentset.bam"
] OutputService: Generating the output XML file
    """
    reads_fn = fn(self.dataset)
    referenceset_fn = fn(self.referenceset)
    o_alignmentset_fn = fn(self.alignmentset)
    task_opts = self.parameters['pbalign']
    options = task_opts.get('options', '')
    algorithmOptions = task_opts.get('algorithmOptions', '')
    wdir, o_alignmentset_fn = os.path.split(o_alignmentset_fn)
    #'--debug', # requires 'ipdb'
    #'--profile', # kinda interesting, but maybe slow?
    #'--algorithmOptions "-minMatch 12 -bestn 10 -minPctSimilarity 70.0"',
    #'--concordant',
    #'--hitPolicy randombest',
    #'--minAccuracy 70.0',
    #'--minLength 50',
    bash = """
pbalign --verbose --nproc 16 {options} --algorithmOptions "{algorithmOptions}" {reads_fn} {referenceset_fn} {o_alignmentset_fn}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_pbalign.sh')
    mkdirs(wdir)
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_gc_scatter(self):
    alignmentset_fn = fn(self.alignmentset)
    referenceset_fn = fn(self.referenceset)
    chunks_fofn_fn = fn(self.out_fofn)
    wdir, chunks_fofn_fn = os.path.split(chunks_fofn_fn)
    bash = """
python -m falcon_polish.mains.run_gc_scatter \
        {alignmentset_fn} \
        {referenceset_fn} \
        {chunks_fofn_fn}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_gc_scatter.sh')
    mkdirs(wdir)
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_gc_gather(self):
    dos = self.inputDataObjs
    ds_out_fn = fn(self.ds_out)
    fastq_out_fn = fn(self.fastq_out)
    wdir, ds_out_fn = os.path.split(ds_out_fn)
    odir, fastq_out_fn = os.path.split(fastq_out_fn)
    assert odir == wdir
    mkdirs(wdir)

    fasta_ds_fofn_fn = os.path.join(wdir, 'fasta.contigset.fofn')
    dset_fns = [fn(v) for k,v in dos.items() if k.startswith('contigset_')]
    open(fasta_ds_fofn_fn, 'w').write('\n'.join(dset_fns) + '\n')

    fastq_fofn_fn = os.path.join(wdir, 'fastq.fofn')
    dset_fns = [fn(v) for k,v in dos.items() if k.startswith('fastq_')]
    open(fastq_fofn_fn, 'w').write('\n'.join(dset_fns) + '\n')

    bash = r"""
python -m falcon_polish.mains.run_gc_gather \
        {fasta_ds_fofn_fn} \
        {fastq_fofn_fn} \
        {ds_out_fn} \
        {fastq_out_fn}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_gc_gather.sh')
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_genomic_consensus(self):
    alignmentset_fn = fn(self.alignmentset)
    referenceset_fn = fn(self.referenceset)
    polished_fastq_fn = fn(self.polished_fastq)
    variants_gff_fn = fn(self.variants_gff)
    consensus_contigset_fn = fn(self.consensus_contigset)
    task_opts = self.parameters['variantCaller']
    options = task_opts.get('options', '')
    if '--alignmentSetRefWindows' not in options:
        options += ' --alignmentSetRefWindows'
    fasta_fn = re.sub(".contigset.xml", ".fasta", consensus_contigset_fn)
    wdir, fasta_fn = os.path.split(fasta_fn)
    odir, consensus_contigset_fn = os.path.split(consensus_contigset_fn)
    assert odir == wdir
    odir, variants_gff_fn = os.path.split(variants_gff_fn)
    assert odir == wdir
    wdir, polished_fastq_fn = os.path.split(polished_fastq_fn)
    assert odir == wdir
    # Possibly we should escape '{options}'
    bash = """
python -m falcon_polish.mains.run_variantCaller --log-level DEBUG --options '{options}' \
        {alignmentset_fn} \
        {referenceset_fn} \
        {polished_fastq_fn} \
        {variants_gff_fn} \
        {fasta_fn} \
        {consensus_contigset_fn}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_gc.sh')
    mkdirs(wdir)
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_polished_assembly_report(self):
    task_opts = self.parameters['pbreports.tasks.summarize_coverage']
    options = task_opts.get('options', '')
    referenceset_fn = fn(self.referenceset)
    gathered_alignmentset_fn = fn(self.gathered_alignmentset)
    polished_fastq_fn = fn(self.polished_fastq)
    report_fn = fn(self.report_json)
    wdir, report_fn = os.path.split(report_fn)
    mkdirs(wdir)
    alignment_summary_gff_fn = 'alignment.summary.gff'
    """
    If necessary, we could call this:
    from pbreports.report.summarize_coverage.summarize_coverage import summarize_coverage
    summarize_coverage(args.aln_set, args.aln_summ_gff, args.ref_set,
                       args.num_regions, args.region_size,
                       args.force_num_regions)
    """
    bash = r"""
python -m pbreports.report.summarize_coverage.summarize_coverage \
        {options} \
        {gathered_alignmentset_fn} \
        {referenceset_fn} \
        {alignment_summary_gff_fn}
python -m pbreports.report.polished_assembly \
        {alignment_summary_gff_fn} \
        {polished_fastq_fn} \
        {report_fn}
""".format(**locals())
    bash_fn = os.path.join(wdir, 'run_report.sh')
    mkdirs(wdir)
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_foo(self):
    log.debug('WARNME1 {!r}'.format(__name__))
    #print repr(self.parameters), repr(self.URL), repr(self.foo1)
    sys.system('touch {}'.format(fn(self.foo2)))

def yield_pipeline_chunk_names_from_json(ifs, key):
    d = json.loads(ifs.read())
    for cs in d['chunks']:
        #chunk_id = cs['chunk_id']
        chunk_datum = cs['chunk']
        yield chunk_datum[key]
def create_tasks_pbalign(chunk_json_pfn, referenceset_pfn, parameters):
    """Create a pbalign task for each chunk, plus a gathering task.
    """
    tasks = list()
    alignmentsets = dict()
    chunk_dir = os.path.dirname(fn(chunk_json_pfn))
    for i, subreadset_fn in enumerate(sorted(yield_pipeline_chunk_names_from_json(open(fn(chunk_json_pfn)), '$chunk.subreadset_id'))):
        wdir = 'run-pbalign-{:02d}'.format(i)
        subreadset_fn = os.path.join(chunk_dir, os.path.basename(subreadset_fn))
        subreadset_pfn = makePypeLocalFile(subreadset_fn)
        alignmentset_pfn = makePypeLocalFile('{wdir}/align.subreads.{i:02d}.alignmentset.xml'.format(**locals()))
        alignmentsets['alignmentsets_{:02d}'.format(i)] = alignmentset_pfn
        """Also produces:
        aligned.subreads.i.alignmentset.bam
        aligned.subreads.i.alignmentset.bam.bai
        aligned.subreads.i.alignmentset.bam.pbi
        """
        make_task = PypeTask(
                inputs = {"chunk_json": chunk_json_pfn,
                          "dataset": subreadset_pfn,
                          "referenceset": referenceset_pfn,},
                outputs = {"alignmentset": alignmentset_pfn,},
                parameters = parameters,
                TaskType = PypeTaskBase,
                URL = "task://localhost/pbalign/{}".format(os.path.basename(subreadset_fn)))
        task = make_task(task_pbalign)
        tasks.append(task)
    alignmentset_pfn = makePypeLocalFile('run-pbalign_gather/aligned.subreads.alignmentset.xml')
    make_task = PypeTask(
            inputs = alignmentsets,
            outputs = {"ds_out": alignmentset_pfn,},
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/pbalign_gather")
    task = make_task(task_pbalign_gather)
    tasks.append(task)
    return tasks, alignmentset_pfn
def mkdirs(d):
    log.debug('mkdir -p {}'.format(d))
    if not os.path.isdir(d):
        os.makedirs(d)
def create_tasks_gc(fofn_pfn, referenceset_pfn, parameters):
    """Create a gc task for each chunk, plus a gathering task.
    Here is the convoluted workflow:
    1. For each gc instance "chunk":
      A. variantCaller writes .fasta
      B. We create a contigset for the .fasta
    2. We keep the contigset output filenames in a FOFN (from run_gc_scatter)
       and pass that to run_gc_gather().
    3. We read each contigset and add them to a gathered ContigSet.
    4. We "consolidate" their underlying .fasta "resources",
       assuming their filenames match except extenion.
    5. Finally, we write the gathered contigset.
    Whew!
    We also gather fastq here, for convenience.
    """
    tasks = list()
    contigsets = dict()
    fastqs = dict()
    for i, alignmentset_fn in enumerate(open(fn(fofn_pfn)).read().split()):
        wdir = 'run-gc-{:02}'.format(i)
        mkdirs(wdir) # Assume CWD is correct.
        alignmentset_pfn = makePypeLocalFile(alignmentset_fn) # New pfn cuz it was not pfn before.
        polished_fastq_pfn = makePypeLocalFile(os.path.join(wdir, 'consensus.fastq'))
        variants_gff_pfn = makePypeLocalFile(os.path.join(wdir, 'variants.gff'))
        consensus_contigset_pfn = makePypeLocalFile(os.path.join(wdir, 'consensus.contigset.xml'))
        """Also produces:
        consensus.fasta
        consensus.fasta.fai

        And note that these files names are important, as pbcoretools gathering expects
        a particular pattern.
        """
        contigsets['contigset_{:02d}'.format(i)] = consensus_contigset_pfn
        fastqs['fastq_{:02d}'.format(i)] = polished_fastq_pfn
        make_task = PypeTask(
                inputs = {"alignmentset": alignmentset_pfn,
                          "referenceset": referenceset_pfn,},
                outputs = {
                    "polished_fastq": polished_fastq_pfn,
                    "variants_gff": variants_gff_pfn,
                    "consensus_contigset": consensus_contigset_pfn,
                },
                parameters = parameters,
                TaskType = PypeTaskBase,
                URL = "task://localhost/genomic_consensus/{}".format(os.path.basename(alignmentset_fn)))
        task = make_task(task_genomic_consensus)
        tasks.append(task)
    contigset_pfn = makePypeLocalFile('run-gc-gather/contigset.xml')
    gathered_fastq_pfn = makePypeLocalFile('run-gc-gather/gathered.fastq')
    inputs = dict(contigsets)
    inputs.update(fastqs)
    log.debug('inputs to gc_gather:{}'.format(pprint.pformat(contigsets)))
    make_task = PypeTask(
            inputs = inputs,
            outputs = {"ds_out": contigset_pfn,
                       "fastq_out": gathered_fastq_pfn,
            },
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/gc_gather")
    task = make_task(task_gc_gather)
    tasks.append(task)
    return tasks, contigset_pfn, gathered_fastq_pfn

def flow(config):
    #import pdb; pdb.set_trace()
    parameters = config
    #exitOnFailure=config['stop_all_jobs_on_failure'] # only matter for parallel jobs
    #wf.refreshTargets(exitOnFailure=exitOnFailure)
    #concurrent_jobs = config["pa_concurrent_jobs"]
    #PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    #wf = PypeThreadWorkflow()
    #wf = PypeWorkflow()
    #wf = PypeWorkflow(job_type='local')
    log.debug('config=\n{}'.format(pprint.pformat(config)))
    wf = PypeWorkflow(job_type=config['hgap']['job_type'])

    dataset_pfn = makePypeLocalFile(config['pbsmrtpipe']['input_files'][0])
    filtered_fasta_pfn = makePypeLocalFile('run-bam2fasta/input.fasta')
    make_task = PypeTask(
            inputs = {"dataset": dataset_pfn, },
            outputs =  {"fasta": filtered_fasta_pfn, },
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/bam2fasta")
    task = make_task(task_bam2fasta)
    wf.addTask(task)
    wf.refreshTargets()

    input_fofn_pfn = makePypeLocalFile('run-falcon/raw_reads.fofn')
    fc_cfg_pfn = makePypeLocalFile('run-falcon/fc.cfg')
    fc_json_config_pfn = makePypeLocalFile("run-falcon/fc.json")
    make_task = PypeTask(
            inputs = {"fasta": filtered_fasta_pfn, },
            outputs = {"fc_cfg": fc_cfg_pfn,
                       "fc_json_config": fc_json_config_pfn,
                       "input_fofn": input_fofn_pfn,
            },
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/prepare_falcon")
    task = make_task(task_prepare_falcon)
    wf.addTask(task)
    wf.refreshTargets()

    # We could integrate the FALCON workflow here, but for now we will just execute it,
    # so we can repeat the sub-flow easily.
    asm_fasta_pfn = makePypeLocalFile('run-falcon/asm.fasta')
    preads_fofn_pfn = makePypeLocalFile('run-falcon/preads.fofn') # for the preassembly report
    length_cutoff_pfn = makePypeLocalFile('run-falcon/0-rawreads/length_cutoff')
    make_task = PypeTask(
            inputs = {"fc_cfg": fc_cfg_pfn,},
            outputs = {"asm_fasta": asm_fasta_pfn,
                       "length_cutoff": length_cutoff_pfn,
                       "preads_fofn": preads_fofn_pfn,
            },
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/falcon")
    task = make_task(task_falcon)
    wf.addTask(task)

    pre_assembly_report_pfn = makePypeLocalFile("pre_assembly_stats.json")
    make_task = PypeTask(
            inputs = {"length_cutoff": length_cutoff_pfn,
                      "raw_reads_fofn": input_fofn_pfn,
                      "preads_fofn": preads_fofn_pfn, },
            outputs = {"pre_assembly_report": pre_assembly_report_pfn, },
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/report_pre_assembly")
    task = make_task(task_report_pre_assembly)
    wf.addTask(task)
    wf.refreshTargets()

    # The reset of the workflow will operate on datasets, not fasta directly.
    referenceset_pfn = makePypeLocalFile('run-fasta2referenceset/asm.referenceset.xml')
    make_task = PypeTask(
            inputs =  {"fasta": asm_fasta_pfn,},
            outputs = {"referenceset": referenceset_pfn,},
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/fasta2referenceset")
    task = make_task(task_fasta2referenceset)
    wf.addTask(task)
    wf.refreshTargets()

    # scatter the subreads for pbalign
    """Produces:
    pbalign_chunk.json
    chunk_subreadset_*.subreadset.xml
    """
    pbalign_chunk_json_pfn = makePypeLocalFile('run-pbalign-scatter/pbalign_chunk.json')
    make_task = PypeTask(
            inputs = {"dataset": dataset_pfn,
                      "referenceset": referenceset_pfn,},
            outputs = {"out_json": pbalign_chunk_json_pfn,},
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/pbalign_scatter")
    task = make_task(task_pbalign_scatter)
    wf.addTask(task)
    wf.refreshTargets()

    # After scattering, we can specify the pbalign jobs.
    tasks, alignmentset_pfn = create_tasks_pbalign(pbalign_chunk_json_pfn, referenceset_pfn, parameters)
    wf.addTasks(tasks)
    wf.refreshTargets()

    # scatter the alignmentset for genomic_consensus (variantCaller)
    """Produces:
    gc.chunks.fofn
    ???*.congitset.xml ???
    """
    gc_chunks_fofn_pfn = makePypeLocalFile('run-gc_scatter/gc.chunks.fofn')
    make_task = PypeTask(
            inputs = {"alignmentset": alignmentset_pfn,
                      "referenceset": referenceset_pfn,},
            outputs = {"out_fofn": gc_chunks_fofn_pfn,},
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/gc_scatter")
    task = make_task(task_gc_scatter)
    wf.addTask(task)
    wf.refreshTargets()

    tasks, contigset_pfn, gathered_fastq_pfn = create_tasks_gc(gc_chunks_fofn_pfn, referenceset_pfn, parameters)
    wf.addTasks(tasks)
    wf.refreshTargets()


    # Final report

    polished_assembly_report_json_pfn = makePypeLocalFile('run-polished-assembly-report/polished_assembly_report.json')
    make_task = PypeTask(
            inputs = {"referenceset": referenceset_pfn,
                      "gathered_alignmentset": alignmentset_pfn,
                      "polished_fastq": gathered_fastq_pfn,},
            outputs = {"report_json": polished_assembly_report_json_pfn,},
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/polished_assembly_report")
    task = make_task(task_polished_assembly_report)
    wf.addTask(task)

    wf.refreshTargets()
    #return
    ##############

    if not os.path.exists('foo.bar1'):
        sys.system('touch foo.bar1')
    foo_fn1 = makePypeLocalFile('foo.bar1')
    foo_fn2 = makePypeLocalFile('foo.bar2')
    make_task = PypeTask(
            inputs = {"foo1": foo_fn1,},
            outputs =  {"foo2": foo_fn2,},
            parameters = parameters,
            TaskType = PypeTaskBase,
            URL = "task://localhost/foo")
    task = make_task(task_foo)
    wf.addTask(task)
    wf.refreshTargets()

"""
also:
genomic_consensus.tasks.gff2vcf-0
genomic_consensus.tasks.gff2bed-0
pbcoretools.tasks.gather_gff-1

/pbi/dept/secondary/siv/smrtlink/smrtlink-alpha/smrtsuite_170220/userdata/jobs_root/000/000114/tasks
"""
