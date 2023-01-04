import sys
import numpy as np
import h5py

def prefix_sum(x):
    prefix = [0]
    for i in range(len(x)):
        back = prefix[-1]
        prefix.append(back + x[i])
    return prefix

def to_csrlist(x):
    data = []
    size = []
    for piece in x:
        data += piece if isinstance(piece, list) else piece.tolist()
        size.append(len(piece))
    offset = prefix_sum(size)
 
    return data, offset 


def remove_duplicates(x):
    return list(set(x))


def topologic_dim(e):
    m = {'vertex': 0, 'line': 1, 'triangle': 2, 'tetra': 3}
    return m[e]

class Mesh:
    def __init__(self, fname, generator='gmsh'):
        if generator == 'gmsh':
            key = 'gmsh:physical'
        elif generator == 'ansa':
            key = 'ansa:geometrical'
        else:
            raise ValueError(f'Unknown generator: {generator}')

        import meshio
        major, minor, _ = meshio.__version__.split('.')
        m = meshio.read(fname)
        self.node = m.points
        dim = self.node.shape[1]
        cells = m.cells
        cell_data = m.cell_data[key]

        self.cell = []
        self.facet = []
        self.cell_data = []
        self.facet_data = []

        if len(cells) != len(cell_data):
            raise ValueError('Unmatched cell and cell data')
        for itype in range(len(cells)):
            cell = cells[itype]
            data = cell_data[itype] 

            if int(major) <= 4:
                cell_dim = topologic_dim(cell[0])
            else:
                cell_dim = cell.dim
            if cell_dim == dim:
                self.cell += cell.data.tolist()
                self.cell_data += data.tolist()
            elif cell_dim == dim - 1:
                self.facet += cell.data.tolist()
                self.facet_data += data.tolist()
        
        self.cell_data = np.array(self.cell_data, dtype=np.int64)
        self.facet_data = np.array(self.facet_data, dtype=np.int64)


    def match(self):
        cells = self.cell #[element for celltype in self.cell for element in celltype] 
        facets = self.facet #[element for celltype in self.facet for element in celltype]

        def collide(f, e):
            elem = set(e)
            return all([i in elem for i in f])
        def orientation(f, e):
            if len(e) != 4:
                raise RuntimeError('Support tetrahedron only')
            unfound = sum(e) - sum(f)
            return np.argwhere(np.array(e) == unfound).squeeze()

        #facet_to_element = [[ielem for ielem in range(len(cells)) if collide(f, cells[ielem])] for f in facets]
        #reverse cell
        v2c = [[] for _ in range(len(cells))]
        for ielem, elem in enumerate(cells):
            for v in elem:
                v2c[v].append(ielem) 
        v2c = [remove_duplicates(m) for m in v2c]

        facet_to_element = []
        for ii, f in enumerate(facets):
            #print(f'{ii}/{len(facets)}')
            connected_cells = []
            for v in f:
                connected_cells += v2c[v]
            values, counts = np.unique(connected_cells, return_counts=True)
            #facet_to_element.append([ielem for ielem in range(len(cells)) if collide(f, cells[ielem])])
            argmax_all = np.argwhere(counts == np.max(counts)).squeeze().tolist()
            if isinstance(argmax_all, int):
                facet_to_element.append([values[argmax_all]])
            elif isinstance(argmax_all, list):
                facet_to_element.append(values[argmax_all])
            else:
                raise RuntimeError(f"Unsupported type: {type(argmax_all)}")
        facet_orientation = [[orientation(facets[ifacet], cells[ielem]) for ielem in c] for ifacet, c in enumerate(facet_to_element)]
        return facet_to_element, facet_to_element

    def adjacency_list(self):
        adjacency = [[] for _ in range(self.node.shape[0])]
        for elem in self.cell:
            for vtx in elem:
                adjacency[vtx] += elem

        return [remove_duplicates(conn) for conn in adjacency]


class MeshPart:
    def __init__(self, mesh, n):
        self.mesh = mesh
        self.epart, self.npart = self.partition(n)
        self.noffset = self.reorder()
    def partition(self, n):
        # concat elements
        #cells = [element for celltype in self.mesh.cell for element in celltype.data] 
        import pymetis
        _, epart, npart = pymetis.part_mesh(n, self.mesh.cell, None, None, gtype=pymetis.GType.DUAL) 
        return epart, npart

    def reorder(self):
        numproc = max(self.npart) + 1
        def mapping():
            npart = np.array(self.npart)

            num_nodes_each_part = [np.argwhere(npart == i).squeeze().shape[0]  for i in range(numproc)]
            num_nodes_offset = [0]
            for nn in num_nodes_each_part:
                back = num_nodes_offset[-1]
                num_nodes_offset.append(back + nn)
            
            new2old = []
            for i in range(numproc):
                new2old += np.argwhere(npart == i).squeeze().tolist()
            old2new = np.zeros((len(new2old),), dtype=np.int32)
            for i, v in enumerate(new2old):
                old2new[v] = i

            assert min(new2old) == 0 and max(new2old) + 1 == len(new2old), 'Wrong map'
            assert min(old2new) == 0 and max(old2new) + 1 == len(old2new), 'Wrong map'
            return new2old, old2new, num_nodes_offset 
        
        new2old, old2new, offset = mapping()
        # set alias
        m = self.mesh
        # reorder nodes
        m.node = m.node[new2old]
        # update elements
        m.cell = [[old2new[n] for n in e] for e in m.cell]
        m.facet = [[old2new[n] for n in e] for e in m.facet]
        return offset 

class HDF5File:
    def __init__(self, fname, mode='w'):
        self.name = fname
        self.handler = h5py.File(fname, mode)

    def _write_rawdata(self, path, data, dtype=None):
        if dtype is None:
            dtype = 'f'

        if isinstance(data, list):
            if len(data) == 0 or isinstance(data[0], list):
                data, offset = to_csrlist(data)
                self.handler.create_dataset(path + '/0', data=data, dtype='i8')
                self.handler.create_dataset(path + '/1', data=offset, dtype='i8')
            elif isinstance(data[0], (int, float)):
                self.handler.create_dataset(path, data=np.array(data))
        elif isinstance(data, np.ndarray):
            self.handler.create_dataset(path, data=data, dtype=dtype)
        else:
            raise TypeError(f'Unsupported datatype: {type(data)}')
           
    def write(self, path, meshpart):
        
        n = max(meshpart.epart) + 1
        f = self.handler
        
        dpath = path[:]
        if path[-1] != '/':
            dpath = path + '/'
        
        write_data = self._write_rawdata
        m = meshpart.mesh
        # write public data (reordered): nodal coordinates, connectivity, adjacency list, element ID
        write_data(dpath + 'node/x', m.node.flatten())
        
        write_data(dpath + 'cell', m.cell, 'i8')

        write_data(dpath + 'facet', m.facet, 'i8')
        _f2e, _for = m.match()
        write_data(dpath + 'facet/f2e', _f2e, 'i8')
        write_data(dpath + 'facet/for', _for, 'i8')

        write_data(dpath + 'node/adjacency', m.adjacency_list(), 'i8')

        write_data(dpath + 'cell/id', m.cell_data, 'i8')
        write_data(dpath + 'facet/id', m.facet_data, 'i8')

        # write private data: epart, node offset
        c = [[] for _ in range(max(meshpart.epart) + 1)]
        for i, e in enumerate(meshpart.epart):
            c[e].append(i)
        write_data(dpath + 'epart', c, 'i8')
        write_data(dpath + 'node/offset', meshpart.noffset, 'i8')

def main(argv):
    input_mesh = argv[1]
    output_mesh = argv[2]
    nparts = int(argv[3])
    print("Reading mesh: ", end='', flush=True)
    mesh = Mesh(input_mesh)
    print("Done", flush=True)
    print("Paritioning: ", end='', flush=True)
    meshpart = MeshPart(mesh, nparts)
    print("Done")

    print("Writing mesh: ", end='', flush=True)
    f = HDF5File(output_mesh, 'w')
    f.write('mesh', meshpart)
    print("Done")

if __name__ == '__main__':
    if len(sys.argv) == 4:
        sys.exit(main(sys.argv))
    else:
        print("python3 main.py [INPUT] [OUTPUT] [NPART]")

