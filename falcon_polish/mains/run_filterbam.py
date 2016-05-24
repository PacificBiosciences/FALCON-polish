"""A simple wrapper for run_filter_dataset
in pbcoretools.tasks.filters
"""
#from pbcoretools.tasks.filters import run_filter_dataset
from pbcoretools.DataSetEntryPoints import parse_filter_list
from pbcore.io import openDataSet
import argparse
import logging
import sys

def run(in_file, out_file, filterstr):
    dataSet = openDataSet(in_file)
    dataSet.updateCounts() # just in case
    filters = parse_filter_list(filterstr.split(','))
    dataSet.filters.addFilter(**filters)
    #log.info("{i} other filters added".format(i=len(filters)))
    dataSet.updateCounts()
    dataSet.write(out_file)

def main(argv=sys.argv):
    description = """Set up filters in the dataset for a BAM resource.
"""
    epilog = """
If you have only '.bam', you will first need to run 'dataset create --type ReferenceSet --generateIndices outdataset inbam*'
to generate the dataset.
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('--logging',
    #        help='.ini or .json config file for Python logging module')
    parser.add_argument('in_file',
        help='Input dataset XML filename. (Must have .bam resource.)')
    parser.add_argument('out_file',
        help='Output dataset XML filename. (Will have .bam resource.)')
    parser.add_argument('filterstr', type=str,
        help='All filters to apply, comma-separated, logically ANDed. E.g. "rq>=.7, length gte 1000, length &lt;= 50000"')
    args = vars(parser.parse_args(argv[1:]))
    #logging.basicConfig()
    #log.warning('RUNNING filterbam: {}'.format(repr(args)))
    run(**args)

if __name__ == "__main__":
    main(sys.argv)
