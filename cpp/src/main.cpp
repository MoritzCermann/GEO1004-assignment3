/*
+------------------------------------------------------------------------------+
|                                                                              |
|                   Ruben Vons | Moritz Cermann | Artemi Kurski                |
|                                                                              |
|                                  2026-05-31                                  |
|                                                                              |
+------------------------------------------------------------------------------+
*/


#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <queue>
#include <ctime>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Dimension.h>
#include <vector>
#include "json.hpp" // nlohmann::json
#include "structs.h" // defines our custom data structures

using json = nlohmann::json;

// Kernel and geometric types
typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  K;
typedef K::Line_3                   Line;
typedef K::Plane_3                  Plane;
typedef K::Point_3                  Point_3;
typedef K::Point_2                  Point_2;
typedef K::Triangle_3               Triangle;
typedef K::Line_3                   Line_3;
typedef K::Line_2                   Line_2;
typedef K::Triangle_3               Triangle;
typedef K::Vector_3                 Vector_3;

// Constrained Delaunay triangulation type

typedef CGAL::Constrained_triangulation_2<K> CDT;



std::map<std::string, Object> objects;


// -------------------------------------------------
//
// -------------------------------------------------
std::vector<Point_3> vertices;





/*// int main(int argc, const char * argv[]) {
//   //-- will read the file passed as argument or twobuildings.city.json if nothing is passed
//   const char* filename = (argc > 1) ? argv[1] : "../../data/nextbk_2b.city.json";
//   std::cout << "Processing: " << filename << std::endl;
//   std::ifstream input(filename);
//   if (!input.is_open()) {
//     std::cerr << "Error: cannot open " << filename << std::endl;
//     return 1;
//   }
//   json j;
//   input >> j;
//   input.close();
//
//   //-- get total number of RoofSurface in the file
//   int noroofsurfaces = get_no_roof_surfaces(j);
//   std::cout << "Total RoofSurface: " << noroofsurfaces << std::endl;
//
//   // list_all_vertices(j);
//
//   // visit_roofsurfaces(j);
//
//   //-- print out the number of Buildings in the file
//   int nobuildings = 0;
//   for (auto& co : j["CityObjects"]) {
//     if (co["type"] == "Building") {
//       nobuildings += 1;
//     }
//   }
//   std::cout << "There are " << nobuildings << " Buildings in the file" << std::endl;
//
//   //-- print out the number of vertices in the file
//   std::cout << "Number of vertices " << j["vertices"].size() << std::endl;
//
//   std::srand(std::time(nullptr));
//
//   //-- add an attribute "volume"
//   for (auto& co : j["CityObjects"]) {
//     if (co["type"] == "Building") {
//       co["attributes"]["volume"] = rand();
//     }
//   }
//
//   json j_copy = j;
//   lod22_to_parent(j_copy);
//
//   // Triangulation here
//   triangulation(j_copy);
//
//   std::cout << "Triangulation complete" << std::endl;
//
//   volume(j_copy);
//
//   std::cout << "Volume calculated" << std::endl;
//
//   calc_roofsurface_area(j_copy);
//
//   std::cout << "Surface area calculated" << std::endl;
//
//   //-- write to disk the modified city model (insert "_out" before ".city.json")
//   std::string outfile = filename;
//   size_t pos = outfile.rfind(".city.json");
//   if (pos != std::string::npos) {
//     outfile.insert(pos, "_out");
//   } else {
//     outfile = "out.city.json";
//   }
//   std::ofstream o(outfile);
//   if (!o.is_open()) {
//     std::cerr << "Error: cannot write " << outfile << std::endl;
//     return 1;
//   }
//   std::cout << "Written to: " << outfile << std::endl;
//   o << j_copy.dump(2) << std::endl;
//   o.close();
//
//   return 0;
// }
//
// // ----------------
// // define functions
// // ----------------
//
// // helper function to get the point coords from the vertex id
// Point_3 get_pt(const json& j, int vertex_id) {
//   auto v = j["vertices"][vertex_id];
//   double x = v[0].get<double>() * j["transform"]["scale"][0].get<double>()
//            + j["transform"]["translate"][0].get<double>();
//   double y = v[1].get<double>() * j["transform"]["scale"][1].get<double>()
//            + j["transform"]["translate"][1].get<double>();
//   double z = v[2].get<double>() * j["transform"]["scale"][2].get<double>()
//            + j["transform"]["translate"][2].get<double>();
//   return Point_3(x, y, z);
// }
//
// void lod22_to_parent(json& j) {
//   // iterate through CityObjects
//   for (auto& co : j["CityObjects"].items()) {
//     // check if current object is child (has attribute = parent or type = BuildingPart)
//     if (co.value()["type"] == "BuildingPart") {
//       // get parent ref
//       std::string parent_name = co.value()["parents"][0];
//
//       // get LoD2.2 geometry and replace geomerty in parent
//       for (auto& geom : co.value()["geometry"]) {
//
//         if (geom["lod"] == "2.2") {
//           json lod22 = geom;
//           json& parent_geom = j["CityObjects"][parent_name]["geometry"];
//           parent_geom.clear();
//           parent_geom.push_back(lod22);
//         }
//       }
//     }
//   }
//
//   // collect ids of BuildingPart
//   std::vector<std::string> to_remove;
//   for (auto& [id, obj] : j["CityObjects"].items()) {
//     if (obj["type"] == "BuildingPart") {
//       // remove children reference
//       std::string parent_name = obj["parents"][0];
//       j["CityObjects"][parent_name].erase("children");
//       // add children id (name) to to_remove list
//       to_remove.push_back(id);
//     }
//   }
//   // erase them
//   for (auto& id : to_remove) {
//     j["CityObjects"].erase(id);
//   }
// }
//
// // Implementation of triangulation
// // extracts vertices
// void extract_vertices(json& j, json& surface,
//                       std::vector<int>& global_ids,
//                       std::unordered_map<int, Point_3>& global_vertices)
// {
//   // loops through rings
//   for (unsigned r = 0; r < surface.size(); r++) {
//     auto& ring = surface[r];
//
//     // loops through ring vertices
//     for (unsigned h = 0; h < ring.size(); h++) {
//
//       int vertex_id = ring[h].get<int>();
//       global_ids.push_back(vertex_id);
//
//       // add if new
//       if (global_vertices.find(vertex_id) == global_vertices.end()) {
//         auto v = j["vertices"][vertex_id];
//
//         double x = v[0].get<double>() * j["transform"]["scale"][0].get<double>()
//                  + j["transform"]["translate"][0].get<double>();
//         double y = v[1].get<double>() * j["transform"]["scale"][1].get<double>()
//                  + j["transform"]["translate"][1].get<double>();
//         double z = v[2].get<double>() * j["transform"]["scale"][2].get<double>()
//                  + j["transform"]["translate"][2].get<double>();
//
//         global_vertices[vertex_id] = Point_3(x, y, z);
//       }
//     }
//   }
// }
//
//
// void triangulation(json& j)
// {
//   std::unordered_map<int, Point_3> global_vertices; // global map
//   json j_out = j;   // output copy
//
//   // loop objs
//   for (auto& co : j["CityObjects"].items()) {
//
//     auto& co_out = j_out["CityObjects"][co.key()];
//
//     // loop geom
//     for (unsigned gi = 0; gi < co.value()["geometry"].size(); gi++) {
//
//       auto& g_in  = co.value()["geometry"][gi];
//       auto& g_out = co_out["geometry"][gi];
//
//       if (g_in["type"] != "Solid") continue;
//
//       // loop shells
//       for (unsigned si = 0; si < g_in["boundaries"].size(); si++) {
//
//         auto& shell_in  = g_in["boundaries"][si];
//         auto& shell_out = g_out["boundaries"][si];
//
//
//         json newShellSurfaces = json::array();   // new tris
//         json newShellSemantics = json::array();  // new sems
//
//         int nSurfaces = shell_in.size();         // fixed count
//
//
//         // loop surfaces
//         for (int k = 0; k < nSurfaces; k++) {
//
//           auto& surface = shell_in[k];   // [rings][verts]
//
//           // collect verts
//           std::vector<int> global_ids;
//           extract_vertices(j, surface, global_ids, global_vertices);
//
//           // ring offsets
//           std::vector<int> ring_offsets;
//           std::vector<int> ring_sizes;
//           int cursor = 0;
//
//           for (unsigned r = 0; r < surface.size(); r++) {
//             auto& ring = surface[r];
//             ring_offsets.push_back(cursor);
//             ring_sizes.push_back((int)ring.size());
//             cursor += (int)ring.size();
//           }
//
//           // edges
//           std::vector<std::pair<int,int>> edges;
//           for (unsigned r = 0; r < ring_offsets.size(); r++) {
//             int offset = ring_offsets[r];
//             int count  = ring_sizes[r];
//             for (int h = 0; h < count; h++) {
//               int curr = offset + h;
//               int next = offset + ((h + 1) % count);
//               edges.push_back({curr, next});
//             }
//           }
//
//           // 3d verts
//           std::vector<Point_3> verts3d;
//           verts3d.reserve(global_ids.size());
//           for (unsigned v = 0; v < global_ids.size(); v++) {
//             verts3d.push_back(global_vertices[global_ids[v]]);
//           }
//
//           // plane fit
//           Plane plane;
//           CGAL::linear_least_squares_fitting_3(
//             verts3d.begin(), verts3d.end(), plane, CGAL::Dimension_tag<0>()
//           );
//
//           // project 2d
//           std::vector<Point_2> verts2d;
//           verts2d.reserve(verts3d.size());
//           for (unsigned v = 0; v < verts3d.size(); v++) {
//             verts2d.push_back(plane.to_2d(verts3d[v]));
//           }
//
//           // CDT
//           CDT cdt;
//           std::vector<CDT::Vertex_handle> vh;
//           vh.reserve(verts2d.size());
//           for (unsigned v = 0; v < verts2d.size(); v++) {
//             vh.push_back(cdt.insert(verts2d[v]));
//           }
//
//           // constraints
//           for (int e = 0; e < (int)edges.size(); e++) {
//             int a = edges[e].first;
//             int b = edges[e].second;
//             // skip degenerate constraints where va = vb
//             if (vh[a] == vh[b]) continue;
//             cdt.insert_constraint(vh[a], vh[b]);
//           }
//
//           // odd-even
//           std::unordered_map<CDT::Face_handle, bool> visited;
//           std::unordered_map<CDT::Face_handle, bool> inside;
//
//           std::queue<CDT::Face_handle> q;
//           CDT::Face_handle start = cdt.infinite_face();
//           visited[start] = true;
//           inside[start]  = false;
//           q.push(start);
//
//           while (!q.empty()) {
//             CDT::Face_handle fh = q.front(); q.pop();
//
//             for (int ei = 0; ei < 3; ei++) {
//               CDT::Face_handle nh = fh->neighbor(ei);
//               if (nh == nullptr) continue;           // truly null: skip
//               if (visited[nh]) continue;             // already done
//               bool crossing = cdt.is_constrained(CDT::Edge(fh, ei));
//               inside[nh] = inside[fh];
//               if (crossing) inside[nh] = !inside[nh];
//               visited[nh] = true;
//               q.push(nh);                            // push regardless of infinite/finite
//             }
//           }
//
//           // map vh → vid
//           std::unordered_map<CDT::Vertex_handle, int> vh_to_vid;
//           for (unsigned v = 0; v < vh.size(); v++) {
//             vh_to_vid[vh[v]] = global_ids[v];
//           }
//
//           json newBoundaries = json::array(); // rings only
//
//           // plane normal
//           Vector_3 ref_n(0,0,0);
//
//           for (unsigned r = 0; r < surface.size(); r++) {
//             auto& ring = surface[r];
//             int n = ring.size();
//
//             // loops through edges
//             for (int i = 0; i < n; i++) {
//               int ia = ring[i].get<int>();
//               int ib = ring[(i+1) % n].get<int>();
//
//               const Point_3& A = global_vertices[ia];
//               const Point_3& B = global_vertices[ib];
//
//               // Newell normal method for polygons
//               ref_n = ref_n + Vector_3(
//                   (A.y() - B.y()) * (A.z() + B.z()),
//                   (A.z() - B.z()) * (A.x() + B.x()),
//                   (A.x() - B.x()) * (A.y() + B.y())
//               );
//             }
//           }
//           //normalized normal vector
//           if (ref_n.squared_length() > 0)
//             ref_n = ref_n / std::sqrt(ref_n.squared_length());
//
//
//
//           // tri loop
//           for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
//
//             if (!inside[f]) continue;
//
//             std::vector<int> tri;
//             tri.reserve(3);
//
//             Point_3 pts[3];
//
//             // collect triangle vertices
//             for (int vi = 0; vi < 3; vi++) {
//               CDT::Vertex_handle vfh = f->vertex(vi);
//               int vid;
//
//               auto it = vh_to_vid.find(vfh);
//               if (it != vh_to_vid.end()) {
//                 vid = it->second;
//               } else {
//                 Point_2 p2 = vfh->point();
//                 Point_3 q3 = plane.to_3d(p2);
//                 j_out["vertices"].push_back({q3.x(), q3.y(), q3.z()});
//                 vid = (int)j_out["vertices"].size() - 1;
//                 vh_to_vid[vfh] = vid;
//                 global_vertices[vid] = q3;
//               }
//
//               tri.push_back(vid);
//               pts[vi] = global_vertices[vid];
//             }
//
//             // compute triangle normal
//             Vector_3 triN = CGAL::cross_product(pts[1] - pts[0], pts[2] - pts[0]);
//
//             //  flip triangle if orientation different from parent face
//             if (triN * ref_n < 0) {
//               std::reverse(tri.begin(), tri.end());
//             }
//
//             newBoundaries.push_back(tri);
//           }
//
//           // parent sem
//           int parentSem = g_in["semantics"]["values"][si][k].get<int>();
//
//           // add tris to shell_out
//           for (auto& ring : newBoundaries) {
//             json surf = json::array();
//             surf.push_back(ring);
//             newShellSurfaces.push_back(surf);
//             newShellSemantics.push_back(parentSem);
//           }
//         }
//
//         // write shell
//         shell_out = newShellSurfaces;
//         g_out["semantics"]["values"][si] = newShellSemantics;
//       }
//     }
//   }
//
//   j = j_out; // commit
// }
//
//
// void volume(json& j)
// {
//   // initialise some util vars to be used for cout later
//   int count = 0;
//   int count_correct= 0;
//   float volume_diff;
//
//   json j_out = j;   // output copy
//
//   // iterate through objects [buildings]
//   for (auto& co : j["CityObjects"].items()) {
//     std::vector<int> volume_ids;
//     std::unordered_map<int, Point_3> volume_vertices;
//
//     auto& co_out = j_out["CityObjects"][co.key()];
//
//
//     // iterate through building's geometry
//     for (unsigned gi = 0; gi < co.value()["geometry"].size(); gi++) {
//       auto& g_in = co.value()["geometry"][gi];
//       // std::cerr << "geom[" << gi << "] type=" << g_in["type"]
//       //           << " lod=" << g_in.value("lod", "N/A")
//       //           << " boundaries.size()=" << g_in["boundaries"].size() << std::endl;
//
//       if (g_in["type"] != "Solid") continue;
//
//       // iterate through building's boundaries
//       for (unsigned si = 0; si < g_in["boundaries"].size(); si++) {
//         auto& shell_in  = g_in["boundaries"][si];
//         // std::cerr << "  shell[" << si << "].size()=" << shell_in.size() << std::endl;
//
//         if (shell_in.size() == 0) continue;
//
//         int nTris = shell_in.size();
//         // extract all the corners of the building
//         for (int k = 0; k < nTris; k++) {
//           auto&tri = shell_in[k];
//           extract_vertices(j, tri, volume_ids, volume_vertices);
//         }
//       }
//     }
//
//     // init list of all vertices of a solid
//     std::vector<Point_3> vertices_list;
//
//     // populate list from json
//     for (auto & [vid, vertex] : volume_vertices) {
//       vertices_list.push_back(vertex);
//     }
//
//     // calculate centroid
//     if (vertices_list.empty()) {
//     std::cerr << "SKIP " << co.key() << ": empty vertices\n";
//     continue;  // skip this building
//     }
//     Point_3 centroid = CGAL::centroid(vertices_list.begin(), vertices_list.end(), CGAL::Dimension_tag<0>());
//
//     // initialise volume
//     float building_volume = 0;
//
//     // again, iterate through building's geometry
//     for (unsigned gi = 0; gi < co.value()["geometry"].size(); gi++) {
//       auto& g_in  = co.value()["geometry"][gi];
//
//       // skip solids
//       if (g_in["type"] != "Solid") continue;
//
//       // iterate through building's triangles
//       for (unsigned si = 0; si < g_in["boundaries"].size(); si++) {
//         auto& shell_in  = g_in["boundaries"][si];
//
//         int nTris = shell_in.size();
//
//         for (int k = 0; k < nTris; k++) {
//           auto&tri = shell_in[k];
//           std::vector<Vector_3> vertices;
//           for (unsigned r = 0; r < tri.size(); r++) {
//             auto& ring = tri[r];
//             // iterate through ring vertices
//             for (unsigned h = 0; h < ring.size(); h++) {
//               int vertex_id = ring[h].get<int>();
//               vertices.push_back(volume_vertices[vertex_id] - centroid);
//             }
//           }
//         float t_volume = vertices[0] * CGAL::cross_product(vertices[1], vertices[2]);
//         building_volume += t_volume;
//         }
//       }
//     }
//     co_out["attributes"]["geo1004_volume"] = std::abs(building_volume / 6.0);
//     double json_volume = co.value()["attributes"]["b3_volume_lod22"].get<double>();
//
//
//     volume_diff = std::abs((std::abs(building_volume / 6.0) - json_volume) / json_volume);
//
//     if (volume_diff < 0.01) {
//       count_correct++;
//     }
//     count++;
//   }
//   j = j_out;
//   std::cout << "Volume matches 3DBAG for " << count_correct << " / " << count << " buildings" << std::endl;
// }
//
//
// void calc_roofsurface_area(json& j) {
//   int count = 0;
//   int count_close = 0;
//
//   for (auto& co : j["CityObjects"].items()) {
//     // in the case we dont have a Building
//     if (co.value()["type"] != "Building") { continue; }
//
//     double total_roof_area = 0.0;
//
//     for (auto& g : co.value()["geometry"]) {
//       // in the case our geometry is not Solid
//       if (g["type"] != "Solid") { continue; }
//
//       for (unsigned i = 0; i < g["boundaries"].size(); i++) {  // iterate through shells (just 1)
//         for (unsigned k = 0; k < g["boundaries"][i].size(); k++) {  // iterate through surfaces (each is 1 triangle)
//
//           // get the corresponding semantic index from "semantics" (this is now for each triangle instead of each plane)
//           int sem_index = g["semantics"]["values"][i][k].get<int>();
//           // std::cout << "k=" << k << " sem_index=" << sem_index
//           //          << " type=" << g["semantics"]["surfaces"][sem_index]["type"] << std::endl;
//
//           // if the current triangle has type RoofSurface
//           if (g["semantics"]["surfaces"][sem_index]["type"] == "RoofSurface")  // or hardcode: if(sem_index == 3)
//           {
//             // access triangles for each roofsurface
//             json& triangle_indices = g["boundaries"][i][k][0];
//
//             Point_3 a = get_pt(j, triangle_indices[0].get<int>());
//             Point_3 b = get_pt(j, triangle_indices[1].get<int>());
//             Point_3 c = get_pt(j, triangle_indices[2].get<int>());
//
//             // calculate the crossproduct
//             double ax = b.x()-a.x(), ay = b.y()-a.y(), az = b.z()-a.z();
//             double bx = c.x()-a.x(), by = c.y()-a.y(), bz = c.z()-a.z();
//
//             double cx = ay*bz - az*by;  // cross product x
//             double cy = az*bx - ax*bz;  // cross product y
//             double cz = ax*by - ay*bx;  // cross product z
//             double area = 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
//             // std::cout << "k=" << k << " sem=" << sem_index << " area=" << area << " running=" << total_roof_area << std::endl;
//
//             total_roof_area += area;
//           }
//         }
//       }
//     }
//     // now save area in obj
//     co.value()["attributes"]["geo1004_total_roof_area"] = total_roof_area;
//
//     // Check area difference
//
//     // just for safety
//     if (co.value()["attributes"].contains("b3_opp_dak_plat") &&
//         co.value()["attributes"].contains("b3_opp_dak_schuin") &&
//         !co.value()["attributes"]["b3_opp_dak_plat"].is_null() &&
//         !co.value()["attributes"]["b3_opp_dak_schuin"].is_null()) {
//
//       double bag_roof = co.value()["attributes"]["b3_opp_dak_plat"].get<double>()
//                       + co.value()["attributes"]["b3_opp_dak_schuin"].get<double>();
//
//       if (bag_roof > 0.0) {
//         double area_diff = std::abs((total_roof_area - bag_roof) / bag_roof);
//         if (area_diff < 0.10) {   // loose 10% tolerance
//           count_close++;
//         }
//       }
//       count++;
//     }
//   }
//
//   std::cout << "Roof area within 10% of 3DBAG (dak_plat + dak_schuin) for "
//             << count_close << " / " << count << " buildings" << std::endl;
// }
//
// // Visit every 'RoofSurface' in the CityJSON model and output its geometry (the arrays of indices)
// // Useful to learn to visit the geometry boundaries and at the same time check their semantics.
// void visit_roofsurfaces(const json& j) {
//   for (auto& co : j["CityObjects"].items()) {
//     for (auto& g : co.value()["geometry"]) {
//       if (g["type"] == "Solid") {
//         for (int i = 0; i < g["boundaries"].size(); i++) {
//           for (int k = 0; k < g["boundaries"][i].size(); k++) {
//             int sem_index = g["semantics"]["values"][i][k];
//             // if (g["semantics"]["surfaces"][sem_index]["type"].get<std::string>().compare("RoofSurface") == 0) {
//             //   std::cout << "RoofSurface: " << g["boundaries"][i][k] << std::endl;
//             // }
//           }
//         }
//       }
//     }
//   }
// }
//
//
// // Returns the number of 'RoofSurface' in the CityJSON model
// int get_no_roof_surfaces(const json& j) {
//   int total = 0;
//   for (auto& co : j["CityObjects"].items()) {
//     for (auto& g : co.value()["geometry"]) {
//       if (g["type"] == "Solid") {
//         for (auto& shell : g["semantics"]["values"]) {
//           for (auto& s : shell) {
//             if (g["semantics"]["surfaces"][s.get<int>()]["type"].get<std::string>().compare("RoofSurface") == 0) {
//               total += 1;
//             }
//           }
//         }
//       }
//     }
//   }
//   return total;
// }
//
//
// // CityJSON files have their vertices compressed: https://www.cityjson.org/specs/1.1.1/#transform-object
// // this function visits all the surfaces and print the (x,y,z) coordinates of each vertex encountered
// void list_all_vertices(const json& j) {
//   for (auto& co : j["CityObjects"].items()) {
//     std::cout << "= CityObject: " << co.key() << std::endl;
//     for (auto& g : co.value()["geometry"]) {
//       if (g["type"] == "Solid") {
//         for (auto& shell : g["boundaries"]) {
//           for (auto& surface : shell) {
//             for (auto& ring : surface) {
//               std::cout << "---" << std::endl;
//               for (auto& v : ring) {
//                 std::vector<int> vi = j["vertices"][v.get<int>()];
//                 double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
//                 double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
//                 double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
//                 std::cout << std::setprecision(2) << std::fixed << v << " (" << x << ", " << y << ", " << z << ")" << std::endl;
//               }
//             }
//           }
//         }
//       }
//     }
//   }
// }*/
