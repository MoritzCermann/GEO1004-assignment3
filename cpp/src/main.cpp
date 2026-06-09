#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

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

const std::string input_file = "../data/OpenHouse2x3-keep-relevant-only.obj";

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

    // loops through all triangles
    for (int t = 0; t < num_tris; ++t) {
        const Triangle& tri = triangles[t];
        CGAL::Bbox_3 tb = tri.bbox();

        // set range of pixels within triangle bounding box
        int min_vx = std::floor((tb.xmin() - origin_x) / n);
        int max_vx = std::floor((tb.xmax() - origin_x) / n);
        int min_vy = std::floor((tb.ymin() - origin_y) / n);
        int max_vy = std::floor((tb.ymax() - origin_y) / n);
        int min_vz = std::floor((tb.zmin() - origin_z) / n);
        int max_vz = std::floor((tb.zmax() - origin_z) / n);

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

                    // voxel center point
                    Point_3 C(
                        origin_x + (vx + 0.5) * n,
                        origin_y + (vy + 0.5) * n,
                        origin_z + (vz + 0.5) * n
                    );

                    // x,y,z rays to check intersection
                    K::Segment_3 segX(C, C + K::Vector_3(n, 0, 0));
                    K::Segment_3 segY(C, C + K::Vector_3(0, n, 0));
                    K::Segment_3 segZ(C, C + K::Vector_3(0, 0, n));

                    // if triangle intersects with any of these segments, fill pixel
                    if (CGAL::do_intersect(tri, segX) ||
                        CGAL::do_intersect(tri, segY) ||
                        CGAL::do_intersect(tri, segZ))
                    {
                        grid(vx, vy, vz) = 1;
                        std::cout << "pixel colored...yay!" << std::endl;
                    }
                }
            }
        }
    }


}

