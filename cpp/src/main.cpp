#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include <queue>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>

#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

using Kernel      = CGAL::Simple_cartesian<double>;
using FT          = Kernel::FT;
using Point3      = Kernel::Point_3;
using Vector3     = Kernel::Vector_3;
using Grid        = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values      = CGAL::Isosurfacing::Value_function_3<Grid>;
using Gradients   = CGAL::Isosurfacing::Gradient_function_3<Grid>;
// using Domain      = CGAL::Isosurfacing::Dual_contouring_domain_3<Grid, Values, Gradients>;
using Domain = CGAL::Isosurfacing::Marching_cubes_domain_3<Grid, Values>;
using SurfaceMesh = CGAL::Surface_mesh<Point3>;

namespace IS  = CGAL::Isosurfacing;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef double FT;
typedef CGAL::Simple_cartesian<FT> K;
typedef K::Point_3 Point_3;
typedef K::Triangle_3 Triangle;


struct Shell {
    std::vector<Triangle> triangles;
};

struct Object {
    std::string id;
    std::vector<Shell> shells;
};

// shamelessly copied from the assignment page
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

Point_3 voxel_to_world(
    int vx, int vy, int vz,
    double min_x, double min_y, double min_z,
    double voxel_size
) {
    double x = min_x + vx * voxel_size;
    double y = min_y + vy * voxel_size;
    double z = min_z + vz * voxel_size;
    return Point_3(x, y, z);
}


std::map<std::string, Object> objects;

const std::string input_file = "../../data/IfcOpenHouse_IFC2x3_simplified.obj";

// Flood fill from all boundary voxels outward — marks exterior as OUTSIDE_ID (2).
void flood_fill_exterior(VoxelGrid& grid, unsigned int outside_id) {
    const unsigned int X = grid.max_x;
    const unsigned int Y = grid.max_y;
    const unsigned int Z = grid.max_z;

    const int dx[] = {1,-1, 0, 0, 0, 0};
    const int dy[] = {0, 0, 1,-1, 0, 0};
    const int dz[] = {0, 0, 0, 0, 1,-1};

    std::queue<std::tuple<int,int,int>> q;

    q.emplace(0, 0, 0);

    while (!q.empty()) {
        auto [cx, cy, cz] = q.front();
        q.pop();
        for (int d = 0; d < 6; ++d) {
            int nx = cx + dx[d];
            int ny = cy + dy[d];
            int nz = cz + dz[d];
            if (nx < 0 || nx >= (int)X) continue;
            if (ny < 0 || ny >= (int)Y) continue;
            if (nz < 0 || nz >= (int)Z) continue;
            if (grid(nx, ny, nz) == 0) {
                grid(nx, ny, nz) = outside_id;
                q.emplace(nx, ny, nz);
            }
        }
    }
}

SurfaceMesh voxelsToMesh(
    const VoxelGrid& vg,
    double voxel_size = 1.0)
{
    const int nx = (int)vg.max_x;
    const int ny = (int)vg.max_y;
    const int nz = (int)vg.max_z;

    // -- Grid ----------------------------------------------------------------
    const CGAL::Bbox_3 bbox(0.0, 0.0, 0.0,
                            nx * voxel_size, ny * voxel_size, nz * voxel_size);
    Grid grid { bbox, CGAL::make_array<std::size_t>(nx, ny, nz) };

    // -- Value function: trilinear interpolation of the 0/1 grid ------------
    // Iso-surface extracted at 0.5 — exactly at the solid/void boundary.
    auto value_fn = [&](const Point3& p) -> FT
    {
        auto clamp = [](int v, int lo, int hi){ return std::max(lo, std::min(hi, v)); };

        double x = p.x() / voxel_size - 0.5;
        double y = p.y() / voxel_size - 0.5;
        double z = p.z() / voxel_size - 0.5;

        int i0 = (int)std::floor(x);
        int j0 = (int)std::floor(y);
        int k0 = (int)std::floor(z);

        double tx = x - i0, ty = y - j0, tz = z - k0;

        double v = 0.0;
        for (int di = 0; di <= 1; ++di)
        for (int dj = 0; dj <= 1; ++dj)
        for (int dk = 0; dk <= 1; ++dk)
        {
            int ci = clamp(i0+di, 0, nx-1);
            int cj = clamp(j0+dj, 0, ny-1);
            int ck = clamp(k0+dk, 0, nz-1);

            double wx = (di == 0) ? 1.0-tx : tx;
            double wy = (dj == 0) ? 1.0-ty : ty;
            double wz = (dk == 0) ? 1.0-tz : tz;
            v += wx * wy * wz * vg(ci, cj, ck);
        }
        return FT(1.0 - v);
    };

    // // -- Gradient: central finite differences of the value function ----------
    // auto gradient_fn = [&](const Point3& p) -> Vector3
    // {
    //     const double h = voxel_size * 0.5;
    //     double gx = (value_fn({p.x()+h, p.y(),   p.z()  })
    //                - value_fn({p.x()-h, p.y(),   p.z()  })) / (2*h);
    //     double gy = (value_fn({p.x(),   p.y()+h, p.z()  })
    //                - value_fn({p.x(),   p.y()-h, p.z()  })) / (2*h);
    //     double gz = (value_fn({p.x(),   p.y(),   p.z()+h})
    //                - value_fn({p.x(),   p.y(),   p.z()-h})) / (2*h);
    //     return Vector3(gx, gy, gz);
    // };

    // -- Domain --------------------------------------------------------------
    Values    values    { value_fn,    grid };
    // Gradients gradients { gradient_fn, grid };
    Domain domain = IS::create_marching_cubes_domain_3(grid, values);

    // -- Run dual contouring -------------------------------------------------
    std::vector<Point3>                  points;
    std::vector<std::vector<std::size_t>> faces;

    IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, FT(0.5), points, faces);

    // -- Convert polygon soup to Surface_mesh --------------------------------
    SurfaceMesh mesh;
    if (!PMP::is_polygon_soup_a_polygon_mesh(faces))
        std::cerr << "Warning: output soup is not a 2-manifold surface.\n";
    else
        PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);

    return mesh;
}

// For every voxel still at 0 (interior air not yet labelled), flood fill
// its connected region with the next available ID (starting at first_id).
// Returns the number of interior regions found.
unsigned int flood_fill_interior(VoxelGrid& grid, unsigned int first_id) {
    const unsigned int X = grid.max_x;
    const unsigned int Y = grid.max_y;
    const unsigned int Z = grid.max_z;

    const int dx[] = {1,-1, 0, 0, 0, 0};
    const int dy[] = {0, 0, 1,-1, 0, 0};
    const int dz[] = {0, 0, 0, 0, 1,-1};

    unsigned int current_id = first_id;

    for (unsigned int x = 0; x < X; ++x) {
        for (unsigned int y = 0; y < Y; ++y) {
            for (unsigned int z = 0; z < Z; ++z) {
                // Only start a new fill from an unlabelled interior voxel
                if (grid(x, y, z) != 0) continue;

                // BFS over this connected pocket
                grid(x, y, z) = current_id;
                std::queue<std::tuple<int,int,int>> q;
                q.emplace(x, y, z);

                while (!q.empty()) {
                    auto [cx, cy, cz] = q.front();
                    q.pop();
                    for (int d = 0; d < 6; ++d) {
                        int nx = cx + dx[d];
                        int ny = cy + dy[d];
                        int nz = cz + dz[d];
                        if (nx < 0 || nx >= (int)X) continue;
                        if (ny < 0 || ny >= (int)Y) continue;
                        if (nz < 0 || nz >= (int)Z) continue;
                        if (grid(nx, ny, nz) == 0) {
                            grid(nx, ny, nz) = current_id;
                            q.emplace(nx, ny, nz);
                        }
                    }
                }

                ++current_id;
            }
        }
    }

    return current_id - first_id; // number of regions found
}

int main() {
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;

    std::vector<Point_3> vertices;

    std::cout << "Reading file..." << std::endl;
    std::ifstream input_stream(input_file);

    std::string current_object_id = "";

    if (!input_stream.is_open()) {
        std::cerr << "Could not open file: " << input_file << std::endl;
        return 1;
    }

    std::string line;
    while (getline(input_stream, line)) {
        // takes line contents
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;

        // Vertex
        if (tag == "v") {
            double x, y, z;
            // assigns x,y,z according to line contents
            ss >> x >> y >> z;
            vertices.emplace_back(x, y, z);
        }

        // Object group
        else if (tag == "g") {
            ss >> current_object_id;
            objects[current_object_id] = Object();
            objects[current_object_id].id = current_object_id;
        }

        // Shell
        else if (tag == "usemtl") {
            if (!current_object_id.empty()) {
                objects[current_object_id].shells.emplace_back();
            }
        }

        // Face (f v//n v//n v//n)
        else if (tag == "f") {
            std::string a, b, c;
            // takes elements as  v//n
            // we dont really need the normals so we filter those out
            ss >> a >> b >> c;

            auto parse_index = [](const std::string& s) {
                size_t slash = s.find('/');
                return std::stoi(s.substr(0, slash)) - 1;
            };

            // processed elements to leave out normals
            int ia = parse_index(a);
            int ib = parse_index(b);
            int ic = parse_index(c);

            Triangle tri(vertices[ia], vertices[ib], vertices[ic]);

            if (!current_object_id.empty() &&
                !objects[current_object_id].shells.empty()) {
                objects[current_object_id].shells.back().triangles.push_back(tri);
            }
        }
    }

    // Summary output
    size_t total_faces = 0;
    for (auto &p : objects) {
        for (auto &shell : p.second.shells) {
            total_faces += shell.triangles.size();
        }
    }

    std::cout << "Finished reading OBJ." << std::endl;
    std::cout << "Stored vertices: " << vertices.size() << std::endl;
    std::cout << "Stored faces:    " << total_faces << std::endl;
    std::cout << "Stored objects:  " << objects.size() << std::endl;


    // bounding box extremes
    CGAL::Bbox_3 bbox = CGAL::bbox_3(vertices.begin(), vertices.end());
    double min_x = bbox.xmin();
    double min_y = bbox.ymin();
    double min_z = bbox.zmin();
    double max_x = bbox.xmax();
    double max_y = bbox.ymax();
    double max_z = bbox.zmax();

    // rows, columns, height
    double n = 0.5;
    unsigned int rows_x = (unsigned int)std::ceil((max_x - min_x) / n) + 2;
    unsigned int rows_y = (unsigned int)std::ceil((max_y - min_y) / n) + 2;
    unsigned int rows_z = (unsigned int)std::ceil((max_z - min_z) / n) + 2;

    // suboptimal turning object map into a list of objects but works better for us
    std::vector<Object> object_list;
    object_list.reserve(objects.size());

    for (auto& [id, obj] : objects) {
        object_list.push_back(obj);
    }

    // array of triangle elements
    std::vector<Triangle> triangles;
    triangles.reserve(total_faces);

    // loops through objects
    int num_objects = object_list.size();
    for (int i = 0; i < num_objects; ++i) {
        Object& obj = object_list[i];
        int num_shells = obj.shells.size();
        // loops through shells
        for (int j = 0; j < num_shells; ++j) {
            Shell& shell = obj.shells[j];
            int num_tris = shell.triangles.size();
            // loops through faces
            for (int k = 0; k < num_tris; ++k) {
                Triangle& tri = shell.triangles[k];
                triangles.push_back(tri);
            }
        }
    }

    double origin_x = min_x - n;
    double origin_y = min_y - n;
    double origin_z = min_z - n;

    VoxelGrid grid(rows_x, rows_y, rows_z);

    int num_tris = triangles.size();

    // voxelise surfaces
    // loops through all triangles
    for (int t = 0; t < num_tris; ++t) {
        const Triangle& tri = triangles[t];
        CGAL::Bbox_3 tb = tri.bbox();

        // set range of voxels within triangle bounding box
        int min_vx = std::floor((tb.xmin() - origin_x) / n);
        int max_vx = std::ceil((tb.xmax() - origin_x) / n);
        int min_vy = std::floor((tb.ymin() - origin_y) / n);
        int max_vy = std::ceil((tb.ymax() - origin_y) / n);
        int min_vz = std::floor((tb.zmin() - origin_z) / n);
        int max_vz = std::ceil((tb.zmax() - origin_z) / n);

        min_vx = std::max(0, min_vx);
        min_vy = std::max(0, min_vy);
        min_vz = std::max(0, min_vz);
        max_vx = std::min((int)grid.max_x - 1, max_vx);
        max_vy = std::min((int)grid.max_y - 1, max_vy);
        max_vz = std::min((int)grid.max_z - 1, max_vz);

        // loop through x,y,z coordinates of voxels
        for (int vx = min_vx; vx <= max_vx; ++vx) {
            for (int vy = min_vy; vy <= max_vy; ++vy) {
                for (int vz = min_vz; vz <= max_vz; ++vz) {

                    // The 8 corners of this voxel
                    double x0 = origin_x + vx * n,       x1 = x0 + n;
                    double y0 = origin_y + vy * n,       y1 = y0 + n;
                    double z0 = origin_z + vz * n,       z1 = z0 + n;

                    CGAL::Bbox_3 bbox(x0, y0, z0, x1, y1, z1);

                    if (CGAL::do_intersect(bbox, tri)) {
                        grid(vx, vy, vz) = 1;
                    }
                }
            }
        }
    }

    SurfaceMesh mesh = voxelsToMesh(grid, n);

    std::cout << mesh.num_vertices() << " vertices, "
              << mesh.num_faces()    << " faces\n";

    CGAL::IO::write_polygon_mesh("output.off", mesh);

    return 0;
}