import falcon_kit.stats_preassembly as func
import StringIO
from nose.tools import assert_equal


def test_cutoff_reads():
    read_lens = [1,2,2,3]
    def check(n, expected):
        got = func.cutoff_reads(read_lens, n)
        assert_equal(expected, got)
    for n, expected in (
            (0, [1,2,2,3]),
            (1, [1,2,2,3]),
            (2, [2,2,3]),
            (3, [3]),
            (4, []),
        ):
        yield check, n, expected

def test_stats_from_sorted_readlengths():
    read_lens = [1,2,2,3]
    stats = func.stats_from_sorted_readlengths(read_lens)
    assert_equal(stats.nreads, 4)
    assert_equal(stats.total, 8)
    assert_equal(stats.n50, 2)
    read_lens = [1,2,3]
    stats = func.stats_from_sorted_readlengths(read_lens)
    assert_equal(stats.nreads, 3)
    assert_equal(stats.total, 6)
    assert_equal(stats.n50, 3)
