This is an implementation of the paper "Memory-Efficient Bijective Parameterizations of Very-Large-Scale Models".

### Dependencies
* [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/download/)
* [Eigen](http://eigen.tuxfamily.org/)
* [Triangle](http://www.cs.cmu.edu/~quake/triangle.html)
* [CGAL](https://www.cgal.org/download.html)
* [PARDISO](https://pardiso-project.org/)
* [AMGCL](https://amgcl.readthedocs.io/en/stable/)
* [glog](https://github.com/google/glog) optional for Output information
* [Qt](http://download.qt.io/archive/qt/) optional for GUI

### Usage

```
exe <mesh_file: *.obj> <method_type: our/outofcore/amgcl> <init_type: tutte/partition>
```
Note that: the input mesh should be a disk-topology tri-mesh. The parameter `num_procs` should be set according to your CPU cores number.