from pyslic.features import featinfo
def test_unique_names():
    A = featinfo.get_slf_names('+')
    assert len(A) == len(set(A))

    B = featinfo.get_names('+')
    assert len(B) == len(set(B))


def test_plus():
    assert len(featinfo.get_slf_names('+')) == len(featinfo.featinfo)

