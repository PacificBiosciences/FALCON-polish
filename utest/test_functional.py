import falcon_polish.functional as func
import StringIO
from nose.tools import assert_equal, assert_raises


def test_calc_cutoff():
    pairs = {0: 1, 1:2, 2:3, 3:4}.items()
    total = 20
    def check(n, expected):
        got = func.calc_cutoff(n, pairs)
        assert_equal(expected, got)
    for n, expected in ((0, 3), (12, 3), (13, 2), (20, 1)):
        yield check, n, expected
    assert_raises(Exception, func.calc_cutoff, 21, pairs)

def test_total_length():
    pairs = {0: 1, 1:2, 2:3, 3:4}.items()
    total = func.total_length(pairs)
    assert_equal(20, total)

def test_fns_from_fofn():
    data = """
a-b
  
c.d
"""
    fofn = StringIO.StringIO(data)
    expected = ['a-b', 'c.d']
    got = list(func.fns_from_fofn(fofn))
    assert_equal(expected, got)

def test_joined_strs():
    def verify(args, expected):
        got = func.joined_strs(*args)
        assert_equal(list(got), list(expected))
    yield verify, ([], 1), []
    yield verify, (['a'], 1), ['a']
    yield verify, (['a'], 2), ['a']
    yield verify, (['a', 'b'], 1), ['ab']
    yield verify, (['a', 'b'], 2), ['a', 'b']
    yield verify, (['a', 'b'], 3), ['a', 'b']
    yield verify, (['a', 'b', 'c'], 1), ['abc']
    yield verify, (['a', 'b', 'c'], 2), ['ab', 'c']
    yield verify, (['a', 'b', 'c'], 3), ['a', 'b', 'c']
    yield verify, (['a', 'b', 'c', 'd'], 1), ['abcd']
    yield verify, (['a', 'b', 'c', 'd'], 2), ['ab', 'cd']
    yield verify, (['a', 'b', 'c', 'd'], 3), ['ab', 'c', 'd']
    yield verify, (['a', 'b', 'c', 'd'], 4), ['a', 'b', 'c', 'd']
