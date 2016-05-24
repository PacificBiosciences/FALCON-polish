"""A simple script to scatter an alignmentset prior to Genomic Consensus.
"""
from pbcore.io import (AlignmentSet, ReferenceSet)
import argparse
import logging
import os
import sys
log = logging.getLogger(__name__)

def run(alignmentset, referenceset, fofn):
    #'python -m pbcoretools.tasks.scatter_alignments_reference alignment_ds ds_reference json_out'
    dir_name = os.getcwd()
    dset = AlignmentSet(alignmentset, strict=True)
    maxChunks = 2
    dset_chunks = dset.split(contigs=True, maxChunks=maxChunks, breakContigs=True)

    # referenceset is used only for sanity checking.
    ReferenceSet(referenceset, strict=True)

    chunk_fns = []
    for i, dset in enumerate(dset_chunks):
        chunk_name = 'chunk_alignmentset_{}.alignmentset.xml'.format(i)
        chunk_fn = os.path.join(dir_name, chunk_name)
        dset.write(chunk_fn)
        chunk_fns.append(chunk_fn)
    with open(fofn, 'w') as ofs:
        for fn in chunk_fns:
            ofs.write('{}\n'.format(fn))
    log.info('Wrote {} chunks into "{}"'.format(len(dset_chunks), fofn))

def main(argv=sys.argv):
    description = """Scatter alignmentsets prior to parallel GenomicConsensus.
"""
    epilog = """
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('--logging',
    #        help='.ini or .json config file for Python logging module')
    parser.add_argument('alignmentset',
        help='Input alignmentset XML filename.')
    parser.add_argument('referenceset',
        help='Input referenceset XML filename.')
    parser.add_argument('fofn',
        help='File of File Names of scattered alignmentsets. Relative paths.')
    args = vars(parser.parse_args(argv[1:]))
    run(**args)

if __name__ == "__main__":
    logging.basicConfig()
    log.setLevel(logging.INFO)
    main(sys.argv)
