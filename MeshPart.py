from pymetis import part_mesh

class MeshPart:
    def __init__(self, mesh, n):
        self.mesh = mesh
        self.epart, self.npart = self.partition(n)

    def partition(self, n):
        assert self.epart is None and self.npart is None, 'Try to partition the mesh twice'

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
        # unfinished 
