#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include <queue>

#include <CGAL/Simple_cartesian.h>

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

    const unsigned int OUTSIDE_ID  = 2;

    const unsigned int INTERIOR_FIRST_ID = 3;

    flood_fill_exterior(grid, OUTSIDE_ID);

    std::map<unsigned int, unsigned int> value_counts;
    for (unsigned int x = 0; x < rows_x; ++x)
        for (unsigned int y = 0; y < rows_y; ++y)
            for (unsigned int z = 0; z < rows_z; ++z)
                value_counts[grid(x, y, z)]++;

    std::cout << "All values in grid after exterior fill:" << std::endl;
    for (auto& [val, count] : value_counts)
        std::cout << "  value " << val << ": " << count << " voxels" << std::endl;

    unsigned int num_rooms = flood_fill_interior(grid, INTERIOR_FIRST_ID);

    // After both fills:
    //   1                          = surface voxel
    //   2                          = exterior air
    //   3, 4, 5, ... (first_id+N)  = individual interior air pockets

    std::cout << "Flood fill complete. Interior regions found: " << num_rooms << std::endl;
}