import falcon_kit.stats_preassembly as func
import StringIO


def test_cutoff_reads():
    read_lens = [1,2,2,3]
    def check(n, expected):
        got = func.cutoff_reads(read_lens, n)
        assert(expected == got)
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
    assert(stats.nreads == 4)
    assert(stats.total == 8)
    assert(stats.n50 == 2)
    read_lens = [1,2,3]
    stats = func.stats_from_sorted_readlengths(read_lens)
    assert(stats.nreads == 3)
    assert(stats.total == 6)
    assert(stats.n50 == 3)
