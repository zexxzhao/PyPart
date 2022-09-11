import meshio

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



