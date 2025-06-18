#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <algorithm>
#include <iomanip>

struct Vec2f {
    float x, y;
};

struct Vec3f {
    float x, y, z;
    float norm() const { return std::sqrt(x*x + y*y + z*z); }
};

Vec3f operator-(const Vec3f& a, const Vec3f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

float dot(const Vec3f& a, const Vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3f cross(const Vec3f& a, const Vec3f& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

struct Face {
    int v[3];
};

void load_obj(const std::string& filename, std::vector<Vec3f>& vertices, std::vector<Face>& faces) {
    std::ifstream in(filename);
    std::string line;
    while (std::getline(in, line)) {
        if (line.substr(0, 2) == "v ") {
            std::istringstream s(line.substr(2));
            Vec3f v;
            s >> v.x >> v.y >> v.z;
            vertices.push_back(v);
        } else if (line.substr(0, 2) == "f ") {
            std::istringstream s(line.substr(2));
            Face f;
            char slash;
            int v_idx, vt_idx, vn_idx;

            for (int i = 0; i < 3; ++i) {
                s >> v_idx;
                s >> slash >> slash >> vn_idx;
                f.v[i] = v_idx - 1;
            }
            faces.push_back(f);
        }
    }
}

void save_obj(const std::string& filename, const std::vector<Vec2f>& uvs, const std::vector<Face>& faces) {
    std::ofstream out(filename);
    out << std::fixed << std::setprecision(6);
    for (const auto& uv : uvs) {
        out << "v " << uv.x << " " << uv.y << " 0.0\n";
    }
    for (const auto& face : faces) {
        out << "f " << face.v[0] + 1 << " " << face.v[1] + 1 << " " << face.v[2] + 1 << "\n";
    }
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <input.obj> <output.obj> <iterations>\n";
        return 1;
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];
    int nb_iter = std::stoi(argv[3]);

    std::vector<Vec3f> vertices;
    std::vector<Face> faces;
    load_obj(input_filename, vertices, faces);

    std::map<std::pair<int, int>, std::vector<int>> edge_to_opposite_verts;
    for (const auto& face : faces) {
        for (int i = 0; i < 3; ++i) {
            int u = face.v[i];
            int v = face.v[(i + 1) % 3];
            int opposite = face.v[(i + 2) % 3];
            if (u > v) std::swap(u, v);
            edge_to_opposite_verts[{u, v}].push_back(opposite);
        }
    }

    std::vector<int> boundary_loop;
    int start_node = -1;
    
    std::map<int, std::vector<int>> boundary_adj;
    for (const auto& edge_pair : edge_to_opposite_verts) {
        if (edge_pair.second.size() == 1) {
            boundary_adj[edge_pair.first.first].push_back(edge_pair.first.second);
            boundary_adj[edge_pair.first.second].push_back(edge_pair.first.first);
        }
    }

    if (!boundary_adj.empty()) {
        start_node = boundary_adj.begin()->first;
        int current_node = start_node;
        int prev_node = -1;

        do {
            boundary_loop.push_back(current_node);
            int next_node = -1;
            for (int neighbor : boundary_adj[current_node]) {
                if (neighbor != prev_node) {
                    next_node = neighbor;
                    break;
                }
            }
            prev_node = current_node;
            current_node = next_node;
        } while (current_node != start_node && current_node != -1);
    }
    
    std::vector<Vec2f> uvs(vertices.size());
    std::vector<bool> is_boundary(vertices.size(), false);
    float boundary_len = 0.0f;
    for (size_t i = 0; i < boundary_loop.size(); ++i) {
        int u = boundary_loop[i];
        int v = boundary_loop[(i + 1) % boundary_loop.size()];
        boundary_len += (vertices[u] - vertices[v]).norm();
        is_boundary[u] = true;
    }

    float current_len = 0.0f;
    for (size_t i = 0; i < boundary_loop.size(); ++i) {
        int u_idx = boundary_loop[i];
        int v_idx = boundary_loop[(i + 1) % boundary_loop.size()];
        float angle = 2.0f * M_PI * current_len / boundary_len;
        uvs[u_idx] = {std::cos(angle), std::sin(angle)};
        current_len += (vertices[u_idx] - vertices[v_idx]).norm();
    }
    
    std::vector<std::vector<int>> adj(vertices.size());
    for(const auto& face : faces) {
        for(int i = 0; i < 3; ++i) {
            adj[face.v[i]].push_back(face.v[(i+1)%3]);
            adj[face.v[(i+1)%3]].push_back(face.v[i]);
        }
    }

    for(auto& v_adj : adj) {
        std::sort(v_adj.begin(), v_adj.end());
        v_adj.erase(std::unique(v_adj.begin(), v_adj.end()), v_adj.end());
    }

    for (int iter = 0; iter < nb_iter; ++iter) {
        std::vector<Vec2f> next_uvs = uvs;
        for (size_t i = 0; i < vertices.size(); ++i) {
            if (!is_boundary[i]) {
                Vec2f avg = {0.0f, 0.0f};
                for (int neighbor_idx : adj[i]) {
                    avg.x += uvs[neighbor_idx].x;
                    avg.y += uvs[neighbor_idx].y;
                }
                if (!adj[i].empty()) {
                    avg.x /= adj[i].size();
                    avg.y /= adj[i].size();
                    next_uvs[i] = avg;
                }
            }
        }
        uvs = next_uvs;
    }

    save_obj(output_filename, uvs, faces);

    return 0;
} 