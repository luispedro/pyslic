import pyslic.image.io.dirtransversal
def test_detect_dirtransversal():
    assert pyslic.image.io.dirtransversal.detect_dirtransversal('tests/data/hela-2d-shell/')
    assert len(pyslic.image.io.dirtransversal.dirtransversal('tests/data/hela-2d-shell/')) == 10
    assert type(pyslic.image.io.dirtransversal.dirtransversal('tests/data/hela-2d-shell/')) is list

