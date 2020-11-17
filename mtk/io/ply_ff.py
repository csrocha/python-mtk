from numpy import array, zeros

class StandfordPlyLoader:
    def __init__(self):
        self.cols = {}
        self.case = {
            '': lambda tokens : True,
            'ply': lambda tokens : self._begin(),
            'format': lambda tokens : self._format(*tokens),
            'comment': lambda tokens : True,
            'element': lambda tokens : self._element(*tokens),
            'property': lambda tokens : self._property(tokens),
            'end_header':lambda tokens : self._end_header(),
        }

    def _begin(self):
        self.state = 'begin'
        return True

    def _format(self, encoding, version):
        assert(self.state == 'begin')
        self.encoding = encoding
        self.version = version
        return True

    def _element(self, elementType, count):
        self.state = elementType
        if elementType == 'vertex':
            self.vertexRows = int(count)
        elif elementType == 'face':
            self.faceRows = int(count)
        self.cols[elementType] = []
        return True

    def _property(self, args):
        assert(self.state != 'begin')
        dataType = args[0]
        name = args[-1]
        if dataType == 'float':
            assert(len(args) == 2)
            self.cols[self.state].append(name)
        elif dataType == 'list':
            assert(len(args) == 4)
            self.cols[self.state].append(name)
        return True

    def _end_header(self):
        assert(self.state != 'begin')
        return False

    def read(self, line):
        tokens = line.split(' ')
        return self.case[tokens[0]](tokens[1:])

def load_ply_header(infile):
    """
    >>> loader = load_ply_header(test_headers[0])
    >>> loader.version == '1.0'
    True
    """
    loader = StandfordPlyLoader()
    for line in infile:
        line = line[:-1]
        if not loader.read(line):
            break
    return loader

def load_ply_vertexs(loader, infile):
    """
    >>> filename = 'data/bunny/reconstruction/bun_zipper.ply'
    >>> fileobj = open(filename)
    >>> loader = load_ply_header(fileobj)
    >>> vertexs = load_ply_vertexs(loader, fileobj)
    """
    vertexs = zeros((loader.vertexRows, len(loader.cols['vertex'])))
    for rowno in xrange(vertexs.shape[0]):
        line = infile.next()[:-1] # Remove eol
        data = [ item for item in line.split(' ') if item != '' ]
        assert(len(data) == vertexs.shape[1])
        vertexs[rowno] = map(float, data)
    return vertexs

def load_ply_faces(loader, infile):
    """
    >>> filename = 'data/bunny/reconstruction/bun_zipper.ply'
    >>> fileobj = open(filename)
    >>> loader = load_ply_header(fileobj)
    >>> vertexs = load_ply_vertexs(loader, fileobj)
    >>> faces = load_ply_faces(loader, fileobj)
    """
    faces = []
    for rowno in xrange(loader.faceRows):
        line = infile.next()[:-1] # Remove eol
        data = [ item for item in line.split(' ') if item != '' ]
        assert(int(data[0])+1 == len(data))
        faces.append(map(int, data[1:]))
    return faces

def load_ply(fileobj):
    """
    >>> filename = 'data/bunny/reconstruction/bun_zipper.ply'
    >>> vertexs, faces = load_ply(open(filename))
    >>> map(min, [vertexs[:,0], vertexs[:,1], vertexs[:,2]])
    >>> map(max, [vertexs[:,0], vertexs[:,1], vertexs[:,2]])
    """
    loader = load_ply_header(fileobj)
    vertexs = load_ply_vertexs(loader, fileobj)
    faces = load_ply_faces(loader, fileobj)
    return vertexs, faces

def test_suite():
    import doctest
    import StringIO
    from numpy import min
    global test_headers
    test_headers = map(StringIO.StringIO, ["""
ply
format ascii 1.0
comment zipper output
element vertex 35947
property float x
property float y
property float z
property float confidence
property float intensity
element face 69451
property list uchar int vertex_indices
end_header
"""])
    return doctest.DocTestSuite()

if __name__ == "__main__":
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())

