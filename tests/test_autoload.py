import pyslic.image.io.autoload
import pyslic.image.io.dirtransversal
def test_autoload():
    imgsauto = pyslic.image.io.autoload.auto_detect_load('tests/data/hela-2d-shell/')
    imgsdir = pyslic.image.io.dirtransversal.dirtransversal('tests/data/hela-2d-shell/')
    assert len(imgsauto) == len(imgsdir)
