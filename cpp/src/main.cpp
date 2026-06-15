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
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>

#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>
#include "json.hpp" // nlohmann::json

using json = nlohmann::json;


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

const unsigned int SOLID_PART_ID = 1;
const unsigned int EXTERIOR_ID = 2;
const unsigned int FIRST_ROOM_ID = 3;


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

//const std::string input_file = "../../data/IfcOpenHouse_IFC2x3_simplified.obj";
const std::string input_file = "../../data/WellnessCenterSama_simplified.obj";

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
    double voxel_size = 0.35)
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

struct ExtractedBoundaries {
    SurfaceMesh exterior;
    std::vector<SurfaceMesh> interiors;
};

ExtractedBoundaries extract_boundaries(SurfaceMesh& mesh)
{
    ExtractedBoundaries result;

    SurfaceMesh::Property_map<SurfaceMesh::Face_index, std::size_t> fccmap =
        mesh.add_property_map<SurfaceMesh::Face_index, std::size_t>("f:CC").first;

    std::size_t num_components = PMP::connected_components(mesh, fccmap);
    std::cout << "Connected components: " << num_components << "\n";

    std::vector<SurfaceMesh> components(num_components);
    PMP::split_connected_components(mesh, components);

    if (components.empty()) return result;

    // Find the outermost component by largest bounding box volume
    auto bbox_volume = [](const SurfaceMesh& m) {
        CGAL::Bbox_3 bb;
        for (auto v : m.vertices()) bb += m.point(v).bbox();
        return (bb.xmax()-bb.xmin()) * (bb.ymax()-bb.ymin()) * (bb.zmax()-bb.zmin());
    };

    std::size_t outer_idx = 0;
    double max_vol = 0;
    for (std::size_t i = 0; i < components.size(); ++i) {
        double vol = bbox_volume(components[i]);
        if (vol > max_vol) { max_vol = vol; outer_idx = i; }
    }

    // Split into exterior and interiors
    for (std::size_t i = 0; i < components.size(); ++i) {
        if (components[i].num_vertices() == 0) continue;
        if (i == outer_idx)
            result.exterior = std::move(components[i]);
        else
            result.interiors.push_back(std::move(components[i]));
    }

    return result;
}

// Extract the outer envelope and per-room surfaces directly from the
// flood-filled voxel grid. For every SOLID (ID = 1) voxel, each of its 6 faces is
// emitted into the mesh of whichever region (exterior or room) lies on
// the other side of that face.
ExtractedBoundaries extract_surfaces(
    const VoxelGrid& grid,
    double voxel_size,
    double origin_x, double origin_y, double origin_z,
    unsigned int num_rooms)
{
    ExtractedBoundaries result;
    result.interiors.resize(num_rooms);

    // one vertex-dedup map per output mesh, keyed by grid-corner index
    std::map<std::tuple<int,int,int>, SurfaceMesh::Vertex_index> ext_vmap;
    std::vector<std::map<std::tuple<int,int,int>, SurfaceMesh::Vertex_index>> room_vmaps(num_rooms);

    auto get_or_add_vertex = [&](SurfaceMesh& mesh,
                                  std::map<std::tuple<int,int,int>, SurfaceMesh::Vertex_index>& vmap,
                                  int cx, int cy, int cz) -> SurfaceMesh::Vertex_index
    {
        auto key = std::make_tuple(cx, cy, cz);
        auto it = vmap.find(key);
        if (it != vmap.end()) return it->second;
        Point_3 p = voxel_to_world(cx, cy, cz, origin_x, origin_y, origin_z, voxel_size);
        auto vd = mesh.add_vertex(p);
        vmap[key] = vd;
        return vd;
    };

    // 4 corners of each of the 6 unit-cube faces, in CCW order such that
    // the resulting normal points OUTWARD in direction d (i.e. into the
    // neighboring voxel). Order: -X, +X, -Y, +Y, -Z, +Z
    static const int face_corners[6][4][3] = {
        {{0,0,0},{0,0,1},{0,1,1},{0,1,0}}, // -X
        {{1,0,0},{1,1,0},{1,1,1},{1,0,1}}, // +X
        {{0,0,0},{1,0,0},{1,0,1},{0,0,1}}, // -Y
        {{0,1,0},{0,1,1},{1,1,1},{1,1,0}}, // +Y
        {{0,0,0},{0,1,0},{1,1,0},{1,0,0}}, // -Z
        {{0,0,1},{1,0,1},{1,1,1},{0,1,1}}, // +Z
    };
    static const int dir_offset[6][3] = {
        {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
    };

    const int X = (int)grid.max_x;
    const int Y = (int)grid.max_y;
    const int Z = (int)grid.max_z;

    for (int vx = 0; vx < X; ++vx)
    for (int vy = 0; vy < Y; ++vy)
    for (int vz = 0; vz < Z; ++vz)
    {
        if (grid(vx, vy, vz) != SOLID_PART_ID) continue;

        for (int d = 0; d < 6; ++d) {
            int nx = vx + dir_offset[d][0];
            int ny = vy + dir_offset[d][1];
            int nz = vz + dir_offset[d][2];

            unsigned int neighbor_val;
            if (nx < 0 || nx >= X || ny < 0 || ny >= Y || nz < 0 || nz >= Z)
                neighbor_val = EXTERIOR_ID; // shouldn't happen given the padding, but be safe
            else
                neighbor_val = grid(nx, ny, nz);

            if (neighbor_val == SOLID_PART_ID) continue; // shared wall between two solids, no surface here

            SurfaceMesh* target = nullptr;
            std::map<std::tuple<int,int,int>, SurfaceMesh::Vertex_index>* vmap = nullptr;

            if (neighbor_val == EXTERIOR_ID) {
                target = &result.exterior;
                vmap   = &ext_vmap;
            } else if (neighbor_val >= FIRST_ROOM_ID) {
                unsigned int room_idx = neighbor_val - FIRST_ROOM_ID;
                if (room_idx >= result.interiors.size()) continue; // safety
                target = &result.interiors[room_idx];
                vmap   = &room_vmaps[room_idx];
            } else {
                continue; // a 0 here would mean flood fill missed a voxel
            }

            SurfaceMesh::Vertex_index vidx[4];
            for (int c = 0; c < 4; ++c) {
                int cx = vx + face_corners[d][c][0];
                int cy = vy + face_corners[d][c][1];
                int cz = vz + face_corners[d][c][2];
                vidx[c] = get_or_add_vertex(*target, *vmap, cx, cy, cz);
            }
            target->add_face(vidx[0], vidx[1], vidx[2], vidx[3]);
        }
    }

    return result;
}


// Used to convert a SurfaceMesh into the MultiSurface CityJSON boundaries.
// The offset here describes the first index of this mesh's vertices in the global list of all vertices
json mesh_to_boundaries(const SurfaceMesh& mesh, std::size_t offset)
{
    json boundaries = json::array();

    for (auto f : mesh.faces()) {
        json ring = json::array();

        // walk the halfedges around the face to collect vertex indices in order
        auto h = mesh.halfedge(f);
        auto h_start = h;
        do {
            ring.push_back(mesh.source(h).idx() + offset);
            h = mesh.next(h);
        } while (h != h_start);

        // CityJSON surface = [ ring ]  (outer ring only, no holes)
        boundaries.push_back({ ring });
    }

    return boundaries;
}

void write_cityjson(
    const SurfaceMesh& exterior,
    const std::vector<SurfaceMesh>& rooms,
    const std::vector<int>& room_storeys,
    int num_storeys,
    const std::string& filename)
{
    std::vector<Point3> all_pts;

    // get exterior vertices
    for (auto v : exterior.vertices())
        all_pts.push_back(exterior.point(v));

    // remember where each rooms vertices start in the global list
    std::vector<std::size_t> room_offset;  // so: room0 starts at v=195, room1 at v= 546...
    for (const SurfaceMesh& room : rooms) {
        room_offset.push_back(all_pts.size());
        for (auto v : room.vertices())
            all_pts.push_back(room.point(v));  // then add the vertices to all_pts
    }

    // CityJSON stores vertices as integers, so we extract scale and translate (bbox minimum)
    // The real coordinate would then be: integet * scale + translate
    double min_x = all_pts[0].x(), min_y = all_pts[0].y(), min_z = all_pts[0].z();
    for (const auto& p : all_pts) {
        min_x = std::min(min_x, p.x());
        min_y = std::min(min_y, p.y());
        min_z = std::min(min_z, p.z());
    }
    const double scale = 0.001;

    json j;
    j["type"] = "CityJSON";
    j["version"] = "2.0";
    j["transform"]["scale"] = { scale, scale, scale };
    j["transform"]["translate"] = { min_x, min_y, min_z };
    j["CityObjects"] = json::object();

    // safe the vertices as long long as the values can get very large
    j["vertices"] = json::array();
    for (const auto& p : all_pts) {
        long long ix = (long long)std::round((p.x() - min_x) / scale);
        long long iy = (long long)std::round((p.y() - min_y) / scale);
        long long iz = (long long)std::round((p.z() - min_z) / scale);
        j["vertices"].push_back({ ix, iy, iz });
    }

    json building;
    building["type"] = "Building";
    building["attributes"] = json::object();

    // BuildingStorey objects are children of Building; rooms are then
    // grouped as children of the storey they belong to (see
    // assign_room_storeys()). Collect each storey's room ids here, then
    // emit the BuildingStorey CityObjects after the rooms below.
    std::vector<json> storey_children(std::max(num_storeys, 0), json::array());

    if (num_storeys > 0) {
        json children = json::array();
        for (int s = 0; s < num_storeys; ++s)
            children.push_back("Storey_" + std::to_string(s));
        building["children"] = children;
    }

    building["geometry"] = json::array();
    if (exterior.num_faces() > 0) {
        json geom;
        geom["type"] = "MultiSurface";
        geom["lod"] = "1";
        geom["boundaries"] = mesh_to_boundaries(exterior, 0);
        building["geometry"].push_back(geom);
    }

    j["CityObjects"]["Building"] = building;

    // BuildingRoom objects
    for (std::size_t i = 0; i < rooms.size(); ++i) {
        std::string id = "Room_" + std::to_string(i);

        // a room's parent is its BuildingStorey if one was assigned, else the Building itself
        std::string parent_id = "Building";
        if (i < room_storeys.size() && room_storeys[i] >= 0) {
            parent_id = "Storey_" + std::to_string(room_storeys[i]);
            storey_children[room_storeys[i]].push_back(id);
        }

        json room_obj;
        room_obj["type"] = "BuildingRoom";
        room_obj["parents"] = { parent_id };
        room_obj["attributes"] = json::object();
        room_obj["geometry"] = json::array();

        if (rooms[i].num_faces() > 0) {
            json geom;
            geom["type"] = "MultiSurface";
            geom["lod"] = "1";
            geom["boundaries"] = mesh_to_boundaries(rooms[i], room_offset[i]);
            room_obj["geometry"].push_back(geom);
        }

        j["CityObjects"][id] = room_obj;
    }

    // BuildingStorey objects, now that we know which rooms belong to each storey
    for (int s = 0; s < num_storeys; ++s) {
        json storey_obj;
        storey_obj["type"] = "BuildingStorey";
        storey_obj["parents"] = { "Building" };
        storey_obj["attributes"] = json::object();
        storey_obj["attributes"]["storeyNumber"] = s;
        if (!storey_children[s].empty())
            storey_obj["children"] = storey_children[s];

        j["CityObjects"]["Storey_" + std::to_string(s)] = storey_obj;
    }

    // write to file
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: cannot write " << filename << std::endl;
        return;
    }
    out << j.dump(2) << std::endl;
    out.close();
    std::cout << "Written to: " << filename << std::endl;
}

// Remove flat/pseudo rooms that are only a few voxel tall
// these are usually floor slabs or the ceiling
// we assign those voxels the SOLID_PART_ID and also have to readjust the room numbering
unsigned int remove_flat_rooms(VoxelGrid& grid, unsigned int num_rooms)
{
    struct BBox {
        int min_z =  std::numeric_limits<int>::max();
        int max_z =  std::numeric_limits<int>::min();
    };

    std::vector<BBox> bbox(num_rooms);

    for (unsigned int x = 0; x < grid.max_x; ++x)
    for (unsigned int y = 0; y < grid.max_y; ++y)
    for (unsigned int z = 0; z < grid.max_z; ++z) {
        unsigned int v = grid(x, y, z);
        if (v < FIRST_ROOM_ID) continue;
        auto& b = bbox[v - FIRST_ROOM_ID];
        b.min_z = std::min(b.min_z, (int)z);
        b.max_z = std::max(b.max_z, (int)z);
    }

    std::vector<bool> is_flat(num_rooms, false);
    unsigned int n_flat = 0;
    for (unsigned int i = 0; i < num_rooms; ++i) {
        if (bbox[i].max_z - bbox[i].min_z + 1 == 1) {
            is_flat[i] = true;
            ++n_flat;
        }
    }

    // build remap old_id -> new_id (0-based), -1 for removed (flat) rooms
    std::vector<int> remap(num_rooms, -1);
    unsigned int new_id = 0;
    for (unsigned int i = 0; i < num_rooms; ++i)
        if (!is_flat[i]) remap[i] = (int)(new_id++);

    for (auto& v : grid.voxels) {
        if (v < FIRST_ROOM_ID) continue;
        unsigned int i = v - FIRST_ROOM_ID;
        if (is_flat[i]) v = SOLID_PART_ID;
        else            v = FIRST_ROOM_ID + (unsigned int)remap[i];
    }

    if (n_flat > 0)
        std::cout << "Removed " << n_flat << " one-voxel-tall pseudo-room(s); "
                  << new_id << " real room(s) remain.\n";

    return new_id;
}

std::vector<int> assign_room_storeys(const VoxelGrid& grid, unsigned int num_rooms)
{
    if (num_rooms == 0) return {};

    const int STOREY_GAP = 2; // voxel layers parameter

    // floor height (lowest z layer) of each room
    std::vector<int> z_min(num_rooms, std::numeric_limits<int>::max());
    for (unsigned int x = 0; x < grid.max_x; ++x)
        for (unsigned int y = 0; y < grid.max_y; ++y)
            for (unsigned int z = 0; z < grid.max_z; ++z) {
                unsigned int v = grid(x, y, z);
                if (v < FIRST_ROOM_ID) continue;
                unsigned int r = v - FIRST_ROOM_ID;
                z_min[r] = std::min(z_min[r], (int)z);
            }

    // sort rooms by floor height, then cut a new storey wherever the gap
    // to the previous floor height is more than STOREY_GAP voxel layers
    std::vector<unsigned int> order(num_rooms);
    for (unsigned int i = 0; i < num_rooms; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](unsigned int a, unsigned int b) { return z_min[a] < z_min[b]; });

    std::vector<int> room_storeys(num_rooms);
    int storey = 0;
    room_storeys[order[0]] = 0;
    for (unsigned int k = 1; k < num_rooms; ++k) {
        if (z_min[order[k]] - z_min[order[k - 1]] > STOREY_GAP) ++storey;
        room_storeys[order[k]] = storey;
    }

    return room_storeys;
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
            std::string rest;
            std::getline(ss, rest);

            size_t start = rest.find_first_not_of(" \t");
            rest = (start == std::string::npos) ? "" : rest.substr(start);
            size_t end = rest.find_last_not_of(" \t\r\n");
            if (end != std::string::npos) rest = rest.substr(0, end + 1);

            current_object_id = rest;
            objects[current_object_id] = Object();
            objects[current_object_id].id = current_object_id;
            objects[current_object_id].shells.emplace_back();
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

    for (auto& [id, obj] : objects) {
        size_t tris = 0;
        for (auto& sh : obj.shells) tris += sh.triangles.size();
        std::cout << id << ": shells=" << obj.shells.size()
                  << ", triangles=" << tris << "\n";
    }

    // -- SET UP VOXEL GRID --
    // bounding box extremes
    CGAL::Bbox_3 bbox = CGAL::bbox_3(vertices.begin(), vertices.end());
    double min_x = bbox.xmin();
    double min_y = bbox.ymin();
    double min_z = bbox.zmin();
    double max_x = bbox.xmax();
    double max_y = bbox.ymax();
    double max_z = bbox.zmax();

    // rows, columns, height
    double n = 0.25;
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

    // -- VOXELISATION --
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

                    //The 8 corners of this voxel
                    double x0 = origin_x + vx * n,       x1 = x0 + n;
                    double y0 = origin_y + vy * n,       y1 = y0 + n;
                    double z0 = origin_z + vz * n,       z1 = z0 + n;

                    CGAL::Bbox_3 bbox(x0, y0, z0, x1, y1, z1);

                    if (CGAL::do_intersect(bbox, tri)) {
                        grid(vx, vy, vz) = SOLID_PART_ID;  // 1
                    }
                    // Point_3 C(
                    //     origin_x + (vx + 0.5) * n,
                    //     origin_y + (vy + 0.5) * n,
                    //     origin_z + (vz + 0.5) * n
                    // );
                    // K::Segment_3 segX(C - K::Vector_3( n/2, 0, 0), C + K::Vector_3( n/2, 0, 0));
                    // K::Segment_3 segY(C - K::Vector_3(0,  n/2, 0), C + K::Vector_3(0,  n/2, 0));
                    // K::Segment_3 segZ(C - K::Vector_3(0, 0,  n/2), C + K::Vector_3(0, 0,  n/2));
                    //
                    // if (CGAL::do_intersect(tri, segX) ||
                    //     CGAL::do_intersect(tri, segY) ||
                    //     CGAL::do_intersect(tri, segZ)
                    //     )
                    // {
                    //     grid(vx, vy, vz) = SOLID_PART_ID;  // 1
                    // }
                }
            }
        }
    }

    // first fill exterior with ID 2
    flood_fill_exterior(grid, EXTERIOR_ID);
    // then fill rooms one by one starting with ID 3
    unsigned int num_rooms = flood_fill_interior(grid, FIRST_ROOM_ID);

    // before extracting geometry we want to remove
    num_rooms = remove_flat_rooms(grid, num_rooms);

    std::cout << "Rooms found: " << num_rooms << "\n";

    std::vector<int> room_storeys = assign_room_storeys(grid, num_rooms);
    int num_storeys = 0;
    for (int s : room_storeys) num_storeys = std::max(num_storeys, s + 1);

    // INITIAL APPROACH: Extract rooms from Mesh
    SurfaceMesh mesh = voxelsToMesh(grid, n);

    std::cout << mesh.num_vertices() << " vertices, "
              << mesh.num_faces()    << " faces\n";
    //ExtractedBoundaries boundaries = extract_boundaries(mesh);


    // FINAL APPROACH: Here extract surfaces from labeled voxels
    ExtractedBoundaries boundaries = extract_surfaces(grid, n, origin_x, origin_y, origin_z, num_rooms);


    std::cout << "Exterior mesh: " << boundaries.exterior.num_vertices() << " vertices\n";
    for (const SurfaceMesh& m : boundaries.interiors) {
        std::cout << "Interior mesh: " << m.num_vertices() << " vertices\n";
    }

    write_cityjson(boundaries.exterior, boundaries.interiors, room_storeys, num_storeys, "../../data/mybuilding.city.json");
}