//
// Created by Wang Songyao on 2024/5/22.
//

#include "../graph/Graph.h"
#include <chrono>
#ifndef INC_NAS_CPP_DYNCOMDETECTOR_H
#define INC_NAS_CPP_DYNCOMDETECTOR_H

typedef std::vector<CommunitySplit> CommunitySplitSnapshots;
typedef std::vector<Graph> GraphSnapshots;
typedef unsigned int TimeStep;
typedef std::vector<unsigned int> TimeSteps;

class DynComDetector {
public:
    GraphSnapshots graph_snapshots;
    CommunitySplitSnapshots communities;
    TimeSteps time_steps;
    VertexList added_vertices;
    VertexList deleted_vertices;
    EdgeList added_edges;
    EdgeList deleted_edges;
    std::vector<VertexList> disappeared_vertices;
    std::vector<std::chrono::milliseconds> time_durations;
public:
    explicit DynComDetector(const std::string& snapshot_path, unsigned int time_step_num);
    DynComDetector();
    CommunitySplit initial_merge();
    void initialize(float merge_ratio);
    CommunitySplit merge_closure_ratio(CommunitySplit & groups, float merge_ratio);
    static float compute_modularity(Graph & g, const CommunitySplit & partition);
    CommunitySplit get_community_snapshot(TimeStep t);
    Graph get_graph_snapshot(TimeStep t);
    CommunitySplit compute_inherited_communities(TimeStep current_step);
    CommunitySplit compute_community_cores(TimeStep current_step, CommunitySplit & inherited_communities);
    static int edge_differ(Graph &g, Community &c, Vertex u);
    CommunitySplit detect_communities_snapshot(TimeStep current_step);
    CommunitySplit initial_merge_snapshot(TimeStep current_step, CommunitySplit & comm_cores);
    CommunitySplit final_merge_snapshot(TimeStep current_step, CommunitySplit & final_community);
    void perform_dynamic_detection();
    void evaluate_modularity();
    void write_result(const std::string & result_path);
};


#endif //INC_NAS_CPP_DYNCOMDETECTOR_H
