import sys
import numpy as np
import h5py

       
def write(fname, meshpart):
    n = max(meshpart.epart) + 1
    
    f = h5py.File(fname, 'w')
    # write public data (reordered): nodal coordinates, connectivity, adjacency list, element ID
    # write private data: epart, node offset

def main(argv):
    input_mesh = argv[1]
    output_mesh = argv[2]
    nparts = int(argv[3])
    mesh = Mesh(input_mesh)
    meshpart = MeshPart(mesh, nparts)

    write(meshpart)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
