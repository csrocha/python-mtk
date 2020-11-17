# -*- coding: utf-8 -*-
"""
Este código está generado a partir de lo descripto en: http://en.wikipedia.org/wiki/Kd-tree
"""
from math import sqrt
from numpy.linalg import norm
from numpy import array

class Node:pass

neginf = float('-inf')
inf = float('inf')

def mod(T, a, v):
    """
    Set a tuple in the a position with the value v.

    >>> T = (10, 20, 30)
    >>> mod(T, 1, 5)
    (10, 5, 30)
    """
    L = list(T)
    L[a] = v
    return tuple(L)

def kdtree(pointList, depth=0, limitMin = None, limitMax = None):
    """
    Generate a kd-tree

    >>> pointList = [(2,3), (5,4), (9,6), (4,7), (8,1), (7,2)]
    >>> tree = kdtree(pointList)
    >>> tree.location
    (7, 2)
    >>> tree.min == (neginf, neginf)
    True
    >>> tree.max == (inf, inf)
    True
    >>> tree.leftChild.location
    (5, 4)
    >>> tree.rightChild.location
    (9, 6)
    >>> tree.leftChild.min
    (-inf, -inf)
    >>> tree.leftChild.max
    (7, inf)
    >>> tree.rightChild.min
    (7, -inf)
    >>> tree.rightChild.max
    (inf, inf)
    >>> tree.leftChild.leftChild.location
    (2, 3)
    >>> tree.leftChild.leftChild.min
    (-inf, -inf)
    >>> tree.leftChild.leftChild.max
    (7, 4)
    >>> tree.leftChild.rightChild.location
    (4, 7)
    >>> tree.leftChild.rightChild.min
    (-inf, 4)
    >>> tree.leftChild.rightChild.max
    (7, inf)
    """
    if not pointList:
        return

    # Select axis based on depth so that axis cycles through all valid values
    k = len(pointList[0]) # assumes all points have the same dimension
    axis = depth % k

    # Set limits of the root
    if limitMin == None:
        limitMin = (neginf,) * k
    if limitMax == None:
        limitMax = (inf,) * k

    # Sort point list and choose median as pivot element
    pointList.sort(key=lambda point: point[axis])
    median = len(pointList)/2 # choose median

    # Create node and construct subtrees
    node = Node()
    node.min = limitMin
    node.max = limitMax
    node.location = pointList[median]
    node.leftChild = kdtree(pointList[0:median], depth+1, limitMin,
                            mod(limitMax, axis, node.location[axis]))
    node.rightChild = kdtree(pointList[median+1:], depth+1,
                             mod(limitMin, axis, node.location[axis]), limitMax)
    return node

def distNode(x, node):
    """
    Calculate distance beetween a coordinate and a node.
    If node is None distance is infinite.
    """
    if not node == None:
        return (node.location, norm(array(x)-array(node.location)))
    else:
        return (None, inf)

def search(x, node, depth=0, best=None):
    """
    Busca el nodo más cerca de la coordenada x

    >>> pointList = [(2,3), (5,4), (9,6), (4,7), (8,1), (7,2)]
    >>> tree = kdtree(pointList)
    >>> search((7,2), tree)
    ((7, 2), 0.0)
    >>> search((7,1), tree)
    ((7, 2), 1.0)
    >>> search((6,2), tree)
    ((7, 2), 1.0)
    >>> search((2,3), tree)
    ((2, 3), 0.0)
    >>> search((8,0), tree)
    ((8, 1), 1.0)
    """
    if not node:
        return best

    k = len(node.location)
    axis = depth % k

    if x[axis] == node.location[axis]:
        return distNode(x, node)
    if x[axis] < node.location[axis]:
        return min([distNode(x, node),
                    search(x, node.leftChild, depth+1,
                           best=distNode(x, node))],
                   key=lambda (c, d): d)
    if x[axis] > node.location[axis]:
        return min([distNode(x, node),
                    search(x, node.rightChild, depth+1,
                           best=distNode(x, node))],
                   key=lambda (c, d): d)

def neighbourgs(x, distance, node, depth=0, l=[]):
    """
    Busca los nodos en un rango de distancia d a x

    >>> pointList = [(2,3), (5,4), (9,6), (4,7), (8,1), (7,2)]
    >>> tree = kdtree(pointList)
    >>> neighbourgs((6,7), 5, tree)
    [(5, 4), (4, 7), (9, 6)]
    >>> neighbourgs((2,8), 3, tree)
    [(4, 7)]
    """
    if not node:
        return l

    d = distNode(x, node)
    if d[1] < distance:
        l.append(d[0])

    k = len(node.location)
    axis = depth % k

    if node.leftChild and \
       node.leftChild.min[axis] < x[axis] + distance and \
       x[axis] - distance < node.leftChild.max[axis]:
        l = l + neighbourgs(x, distance, node.leftChild, depth=depth+1, l=[])
    if node.rightChild and \
       node.rightChild.min[axis] < x[axis] + distance and \
       x[axis] - distance < node.rightChild.max[axis]:
        l = l + neighbourgs(x, distance, node.rightChild, depth=depth+1, l=[])

    return l

def test_suite():
    import doctest
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

