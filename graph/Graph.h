//
// Created by Wang Songyao on 2024/5/20.
//
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <filesystem>
#ifndef INC_NAS_CPP_GRAPH_H
#define INC_NAS_CPP_GRAPH_H

typedef unsigned int Vertex;
typedef std::pair<Vertex, Vertex> Edge;
typedef std::map<Vertex, unsigned int> DegreeMap;
typedef std::vector<Vertex> VertexList, NeighborList, VertexGroup, Community;
typedef std::map<Vertex, VertexList> AdjList;
typedef std::vector<VertexList> CommunitySplit;
typedef std::set<Vertex> VertexSet;
typedef std::vector<std::string> VertexName;
typedef std::vector<Edge> EdgeList;
typedef std::set<Edge> EdgeSet;

/**
 * 这是一个自定义的Graph类，这里默认是Undirected Graph.
 */
class Graph {
private:
    unsigned int vertex_num{};  // 节点数
    unsigned int edge_num{};  // 边数
    AdjList adj_list;  // 邻接表
    DegreeMap degree_map;  // 节点度
public:
    explicit Graph();
    explicit Graph(const std::string& file_name);
    [[nodiscard]] unsigned int get_vertex_num() const;
    [[nodiscard]] unsigned int get_edge_num() const;
    void add_edge(Vertex u1, Vertex u2);
    bool has_edge(Vertex u1, Vertex u2);
    void print_graph();
    VertexList get_neighbor_list(Vertex u);
    VertexSet get_neighbor_set(Vertex u);
    static float similarity(Graph & g, Vertex u1, Vertex u2);
    static float group_similarity(Graph & g, const VertexGroup& group1, const VertexGroup& group2);
    static float ratio(Graph& g, const VertexGroup& group);
    DegreeMap get_degree_map();
    unsigned int get_degree(Vertex u);
    bool has_vertex(Vertex u1);
    void add_vertex(Vertex u);
    VertexList get_vertex_list();
    VertexSet get_vertex_set();
    static float compute_similarity(Graph & g, Vertex u1, Vertex u2);
    EdgeList get_edge_list();
};

typedef Graph SnapshotGraph;


#endif //INC_NAS_CPP_GRAPH_H
