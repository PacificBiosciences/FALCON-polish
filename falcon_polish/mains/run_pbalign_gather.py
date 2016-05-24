"""A simple wrapper for gather_alignmentset
from pbcoretools.chunking.gather
"""
from pbcoretools.chunking.gather import gather_alignmentset
import argparse
import logging
import os
import sys
log = logging.getLogger(__name__)

def run(idatasets, odataset):
    odir = os.path.dirname(odataset)
    if odir and not os.path.isdir(odir):
        # I think pbcoretools tend to expect the output dir to exist.
        os.makedirs(odir)
    dset_out_fn = odataset
    dset_fns = open(idatasets).read().strip().split()
    gather_alignmentset(dset_fns, dset_out_fn)

def main(argv=sys.argv):
    description = """Gather results of parallel pbalign (blasr).
"""
    epilog = """
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('idatasets',
        help='File of input alignmentset XML filenames, whitespace delimited.')
    parser.add_argument('odataset',
        help='Output alignmentset XML filename.')
    args = vars(parser.parse_args(argv[1:]))
    run(**args)

if __name__ == "__main__":
    logging.basicConfig()
    log.setLevel(logging.INFO)
    main(sys.argv)
