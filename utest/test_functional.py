import falcon_polish.functional as func
import StringIO
import pytest


def test_calc_cutoff():
    pairs = {0: 1, 1:2, 2:3, 3:4}.items()
    total = 20
    def check(n, expected):
        got = func.calc_cutoff(n, pairs)
        assert(expected == got)
    for n, expected in ((0, 3), (12, 3), (13, 2), (20, 1)):
        yield check, n, expected
    with pytest.raises(Exception) as excinfo:
        func.calc_cutoff(21, pairs)

def test_total_length():
    pairs = {0: 1, 1:2, 2:3, 3:4}.items()
    total = func.total_length(pairs)
    assert(20 == total)

def test_fns_from_fofn():
    data = """
a-b
  
c.d
"""
    fofn = StringIO.StringIO(data)
    expected = ['a-b', 'c.d']
    got = list(func.fns_from_fofn(fofn))
    assert(expected == got)

def test_joined_strs():
    def verify(args, expected):
        got = func.joined_strs(*args)
        assert(list(got) == list(expected))
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
