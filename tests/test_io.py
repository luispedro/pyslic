import pyslic
def tests_emptydir():
    assert not pyslic.image.io.detect_ksr_dir('tests/data/emptydir')
    assert not pyslic.image.io.detect_ic100dir('tests/data/emptydir')

