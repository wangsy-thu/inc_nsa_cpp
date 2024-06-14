//
// Created by Wang Songyao on 2024/5/20.
//

#include "Graph.h"

unsigned int Graph::get_vertex_num() const {
    return this->vertex_num;
}

unsigned int Graph::get_edge_num() const {
    return this->edge_num;
}

void Graph::add_vertex(Vertex u) {
    if(has_vertex(u)) {
        std::cerr << "Vertex " << u << " already exist!";
        return;
    }
    this->adj_list.insert(std::make_pair(u, VertexList(0)));
    this->vertex_num ++;
}

void Graph::add_edge(Vertex u1, Vertex u2) {
    // 首先检查 u1 和 u2 是否存在
    if(!has_vertex(u1)) {
        this->add_vertex(u1);
    }
    if(!has_vertex(u2)) {
        this->add_vertex(u2);
    }

    // 再检查这条无向边是否存在
    if(!this->has_edge(u1, u2) && !this->has_edge(u2, u1)) {
        this->adj_list[u1].emplace_back(u2);
        this->adj_list[u2].emplace_back(u1);
        this->edge_num ++;
    }
}

bool Graph::has_edge(Vertex u1, Vertex u2){
    // Find u1 and u2.
    if(!has_vertex(u1) || !has_vertex(u2)){
        return false;
    }

    // u1 <-> u2
    return std::any_of(this->adj_list[u1].begin(), this->adj_list[u1].end(), [&](int item) {
        return item == u2;
    });
}

void Graph::print_graph() {
    std::cout << "Graph Size: " << this->vertex_num << std::endl;
    for(int u = 0; u < this->vertex_num; u ++) {
        std::cout << u << ": ";
        for(const unsigned int & v : this->adj_list[u]) {
            std::cout << v << " -> ";
        }
        std::cout << "END" << std::endl;
    }
}

bool Graph::has_vertex(Vertex u1) {
    return this->adj_list.find(u1) != this->adj_list.end();
}

Graph::Graph(const std::string &file_name) {
    this->adj_list = std::map<Vertex, std::vector<Vertex>>();
    std::ifstream edges_file(file_name);
    if(!edges_file.is_open()){
        std::cerr << "Edge list file open failed!" << std::endl;
        return;
    }

    std::string buffer;

    while(std::getline(edges_file, buffer)) {
        std::istringstream iss(buffer);
        int src, dest;
        char t, r;
        iss >> src >> dest;
        this->add_edge(src - 1, dest - 1);
    }
}

VertexList Graph::get_neighbor_list(Vertex u) {
    VertexList v_list;
    if(!this->has_vertex(u)) {
        std::cerr << "Vertex " << u << " not exist" << std::endl;
        return v_list;
    }
    return this->adj_list[u];
}

VertexSet Graph::get_neighbor_set(Vertex u) {
    VertexSet v_set;
    if(!this->has_vertex(u)) {
        std::cerr << "Vertex " << u << " not exist 1" << std::endl;
        return v_set;
    }
    for(auto v : this->adj_list[u]) {
        v_set.insert(v);
    }
    return v_set;
}

float Graph::similarity(Graph &g, Vertex u1, Vertex u2) {
    VertexSet v_set1 = g.get_neighbor_set(u1);
    VertexSet v_set2 = g.get_neighbor_set(u2);
    VertexSet v_set_intersection;
    std::set_intersection(v_set1.begin(), v_set1.end(), v_set2.begin(), v_set2.end(),
                          std::inserter(v_set_intersection, v_set_intersection.begin()));
    return (float)((float) v_set_intersection.size() / (float) (v_set1.size() + v_set2.size() - v_set_intersection.size()));
}

float Graph::group_similarity(Graph &g, const VertexGroup& group1, const VertexGroup& group2) {
    float total_similarity = 0.0;
    for(Vertex v1 : group1){
        for(Vertex v2 : group2) {
            if(g.has_edge(v1, v2)) {
                total_similarity += similarity(g, v1, v2);
            }
        }
    }
    return total_similarity / (float)(group2.size());
}

float Graph::ratio(Graph &g, const VertexGroup &group) {
    float in_out_ratio;
    unsigned int inner_edge_num = 0, outer_edge_num = 0;
    VertexSet temp_set(group.begin(), group.end());
    for (Vertex v : group) {
        VertexGroup v_neighbor_group = g.get_neighbor_list(v);
        for(auto v_neighbor : v_neighbor_group) {
            if(temp_set.count(v_neighbor) > 0) {
                inner_edge_num += 1;
            } else {
                outer_edge_num += 1;
            }
        }
    }

    if(outer_edge_num == 0) {
        in_out_ratio = 1000;
    } else {
        in_out_ratio = (float) ((float) inner_edge_num / (float) (outer_edge_num * 2)) *
                ((float) group.size() / (float) g.get_vertex_num());
    }
    return in_out_ratio;
}

DegreeMap Graph::get_degree_map() {
    DegreeMap deg_map;
    for(Vertex u = 0; u < this->vertex_num; u++) {
        deg_map.insert(std::make_pair(u, this->adj_list[u].size()));
    }
    return deg_map;
}

unsigned int Graph::get_degree(Vertex u) {
    if(!has_vertex(u)) {
        std::cerr << "Vertex not exist!" << std::endl;
        return 0;
    }
    return this->adj_list[u].size();
}

VertexList Graph::get_vertex_list() {
    VertexList result_list;
    for(const auto& p : this->adj_list) {
        result_list.emplace_back(p.first);
    }
    return result_list;
}

VertexSet Graph::get_vertex_set() {
    VertexSet result_set;
    for(const auto& p : this->adj_list) {
        result_set.insert(p.first);
    }
    return result_set;
}

EdgeList Graph::get_edge_list() {
    EdgeList result;
    for(const auto & p : this->adj_list) {
        for(Vertex v : p.second) {
            if(p.first <= v) {
                result.emplace_back(std::make_pair(p.first, v));
            }
        }
    }
    return result;
}

float Graph::compute_similarity(Graph &g, Vertex u1, Vertex u2) {
    VertexSet v_set1 = g.get_neighbor_set(u1);
    VertexSet v_set2 = g.get_neighbor_set(u2);
    v_set1.insert(u1);
    v_set2.insert(u2);
    VertexSet v_set_intersection;
    std::set_intersection(v_set1.begin(), v_set1.end(), v_set2.begin(), v_set2.end(),
                          std::inserter(v_set_intersection, v_set_intersection.begin()));
    auto s12 = (float) v_set_intersection.size();
    return (float)(s12 / (float) sqrt((float)(g.get_degree(u1) * g.get_degree(u2))));
}

Graph::Graph() = default;
