import sys
import numpy as np
import meshio
import h5py

class Mesh:
    def __init__(self, fname, generator='gmsh'):
        if generator == 'gmsh':
            key = 'gmsh:physical'
        elif generator == 'ansa':
            key = 'ansa:geometrical'
        else:
            raise ValueError(f'Unknown generator: {generator}')

        m = meshio.read(fname)
        self.node = m.points
        dim = self.node.shape[1]
        cells = m.cells
        cell_data = m.cell_data[key]

        self.cell = []
        self.facet = []
        self.cell_data = []
        self.facet_data = []

        assert len(cells) == len(cell_data), 'Unmatched cell and cell data'
        for itype in range(len(cells)):
            cell = cells[itype]
            data = cell_data[itype] 

            if cell.dim == dim:
                self.cell.append(cell)
                self.cell_data.append(data)
            elif c.dim == dim - 1:
                self.facet.append(cell)
                self.facet_data.append(data)
        

    def match(self):
        cells = [element for celltype in self.cell for element in celltype] 
        facets = [element for celltype in self.facet for element in celltype]

        def collide(f, e):
            return all([any(e == i) for i in f])
        def orientation(f, e):
            assert len(e) == 4, 'Support tetrahedron only'
            unfound = sum(e) - sum(f)
            return np.argwhere(e == unfound).squeeze()

        facet_to_element = [[ielem for ielem in range(cells) if collide(f, cells[ielem])] for f in facets]
        facet_orientation = [[orientation(facets[ifacet], cells[ielem]) for ielem in c] for ifacet, c in enumerate(facet_to_element)]
        return facet_to_element, facet_to_element

    def adjacency_list(self):
        adjacency = [[] for _ in range(self.node.shape[0])]
        for elem in self.cell:
            for vtx in elem:
                adjacency[vtx].append(elem)

        return [list(set(conn)) for conn in adjacency]


class MeshPart:
    def __init__(self, mesh, n):
        self.mesh = mesh
        self.epart, self.npart = self.partition(n)

    def partition(self, n):
        assert self.epart is None and self.npart is None, 'Try to partition the mesh twice'

        from pymetis import part_mesh
        # concat elements
        cells = [element for celltype in self.mesh.cell for element in celltype] 
        _, epart, npart = part_mesh(n, cells) 
        return epart, npart

    def part(self, rank):   
        nodes = []
        n_owned, n_ghoested, n_all = 0, 1, 2
        nodes.append(np.argwhere(self.npart == rank).squeeze().tolist())
        allnodes = list(set([n for n in c for c in self.mesh.cells if self.epart[c] == rank]))
        nodes.append([n for n in allnodes if self.npart[n] != rank])
        nodes.append(node[0] + node[1])
        num = [len(nodes[n_owned]), len(nodes[n_ghosted]), len(allnodes)]
        assert num[0] + num[1] == num[2], "Miss/Repeat nodes"
         
        x = self.mesh.node[node[n_all]]

        element = np.argwhere(self.epart == rank).squeeze().tolist()
        
        #find facet
        
def write(fname, meshpart):
    n = max(meshpart.epart) + 1
    
    f = h5py.File(fname, 'w')
    # write public data (reordered): nodal coordinates, connectivity, adjacency list
    # write private data: 

def main(argv):
    input_mesh = argv[1]
    output_mesh = argv[2]
    nparts = int(argv[3])
    mesh = Mesh(input_mesh)
    meshpart = MeshPart(mesh, nparts)

    write(meshpart)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
