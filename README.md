# PyPart
A Python3 tool to partition meshes for XXX
## External dependencies
0. Python (>= 3.8)
1. numpy
2. meshio
3. h5py
4. PyMetis
## Installation
Install wheel first if you don't have it
```bash
pip3 install wheel
```
Then,
```bash
pip3 install . 
```
# Usage
Use `pypart --help` to print the usage information. It typically goes like
```bash
pypart cylinder.msh cylinder.h5 16
```
where `cylinder.msh` is the raw mesh file and `cylinder.h5` is the partitioned one, which is divided into 16 submeshes.
## Mesh file hierarchy
(TODO: Add more descriptions)
```
mesh/
├── cell
│   ├── 0
│   ├── 1
│   └── id
├── epart
│   ├── 0
│   └── 1
├── facet
│   ├── 0
│   ├── 1
│   ├── f2e
│   │   ├── 0
│   │   └── 1
│   ├── for
│   │   ├── 0
│   │   └── 1
│   └── id
└── node
    ├── adjacency
    │   ├── 0
    │   └── 1
    ├── offset
    └── x
```
## Known issues
### limited element types
#### 3D
- [x] tetrahedron
  - [x] triangle
- [x] prism
  - [x] triangle
  - [ ] quadrilateral
- [ ] hexahedron
  - [ ] quadrilateral
- [ ] B-spline
- [ ] Nurbs
#### 2D
- [ ] triangle
  - [ ] line
- [ ] quadrilateral
  - [ ] line
#### 1D
- [ ] line
  - [ ] point


