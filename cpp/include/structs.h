//
// Here we define all our Data Structures to keep it organized
//

#ifndef ASSIGNMENT3_STRUCTS_H
#define ASSIGNMENT3_STRUCTS_H

#include <vector>
#include <queue>
#include <string>
#include <CGAL/Simple_cartesian.h>

// -------------------------------------------------
// Define Data Structures
// -------------------------------------------------


// Kernel and geometric types
typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  Kernel;

// reading the obj
struct Shell {
  std::vector<Kernel::Triangle_3> triangles;
};

struct Object {
  std::string id;
  std::vector<Shell> shells;
};

// creating the voxel grid
struct VoxelGrid {
  std::vector<unsigned int> voxels;
  unsigned int max_x, max_y, max_z;

  VoxelGrid(unsigned int x, unsigned int y, unsigned int z) {
    max_x = x;
    max_y = y;
    max_z = z;
    unsigned int total_voxels = x*y*z;
    voxels.reserve(total_voxels);
    for (unsigned int i = 0; i < total_voxels; ++i) voxels.push_back(0);
  }

  unsigned int &operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
    assert(x < max_x);
    assert(y < max_y);
    assert(z < max_z);
    return voxels[x + y*max_x + z*max_x*max_y];
  }

  unsigned int operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z) const {
    assert(x < max_x);
    assert(y < max_y);
    assert(z < max_z);
    return voxels[x + y*max_x + z*max_x*max_y];
  }
};

#endif //ASSIGNMENT3_STRUCTS_H