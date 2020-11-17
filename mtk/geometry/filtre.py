from mtk.geometry.vol import BasicVolume
from numpy.linalg import norm
from numpy import zeros, sum

def ConnollyFiltre(vol, radius, value=1):
    """
    Fill Connolly filtre in a volume.

    >>> V = BasicVolume((-5,-5,-5),(1,1,1),zeros((10,10,10)))
    >>> V = ConnollyFiltre(V, 4)
    >>> sum(V._data) == 134
    True
    """
    center = vol.center
    boxmin = center - radius
    boxmax = center + radius
    for x, i, v in vol.boxiterator(boxmin, boxmax):
        d = norm(x - center)
        if radius >= d and d > radius - vol.delta[0]:
            vol[i] = value
    return vol

def CrochaFiltre(vol, radius, value=1):
    """
    Fill Crocha filtre in a volume.

    >>> V = BasicVolume((-5,-5,-5),(1,1,1),zeros((10,10,10)))
    >>> V = CrochaFiltre(V, 4)
    >>> sum(V._data) == 257
    True
    """
    center = vol.center
    boxmin = center - radius
    boxmax = center + radius
    for x, i, v in vol.boxiterator(boxmin, boxmax):
        d = norm(x - center)
        if radius >= d:
            vol[i] = value
    return vol

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

