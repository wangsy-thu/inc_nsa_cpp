//
// Created by Wang Songyao on 2024/5/22.
//

#include "DynComDetector.h"

void DynComDetector::initialize(float merge_ratio) {
    auto start = std::chrono::high_resolution_clock::now();
    CommunitySplit initial_groups = initial_merge();
    CommunitySplit final_groups = merge_closure_ratio(initial_groups, merge_ratio);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // 将时间添加到vector容器中
    this->time_durations.emplace_back(duration);
    this->communities.emplace_back(final_groups);
    // compute disappeared vertices;
    for(TimeStep t = 1; t < this->time_steps.size(); t++) {
        this->disappeared_vertices.emplace_back(0);
        Graph &previous_graph = this->graph_snapshots[t - 1];
        Graph &current_graph = this->graph_snapshots[t];
        VertexSet previous_vertex_set = previous_graph.get_vertex_set();
        VertexSet current_vertex_set = current_graph.get_vertex_set();
        for(Vertex pre_v : previous_vertex_set) {
            if(current_vertex_set.count(pre_v) == 0) {
                this->disappeared_vertices[t - 1].emplace_back(pre_v);
            }
        }
    }
}

CommunitySplit DynComDetector::initial_merge() {
    Graph &g = this->graph_snapshots[0];
    VertexList sorted_vertex_list = g.get_vertex_list();
    CommunitySplit merge;

    std::sort(sorted_vertex_list.begin(), sorted_vertex_list.end(),
              [&](auto u1, auto u2) {
        return g.get_degree(u1) > g.get_degree(u2);
    });

    int for_count = 0;
    std::vector<bool> is_deleted(sorted_vertex_list.size(), false);
    for(unsigned int vid = 0; vid < sorted_vertex_list.size(); vid++) {
        if(is_deleted[vid]) {
            continue;
        }
        for_count ++;
        Vertex v = sorted_vertex_list[vid];
        float max_similarity = 0;
        Vertex max_similarity_vertex = 0;
        for(Vertex u : g.get_neighbor_list(v)) {
            float similarity = Graph::similarity(g, v, u);
            if(similarity >= max_similarity) {
                max_similarity = similarity;
                max_similarity_vertex = u;
            }
        }
        if(merge.empty()) {
            if(max_similarity_vertex == 68) {
                std::cout << "Find Error Place!" << std::endl;
            }
            VertexList t = {v, max_similarity_vertex};
            merge.emplace_back(t);
            auto max_sim_vid = std::find(sorted_vertex_list.begin(), sorted_vertex_list.end(), max_similarity_vertex)
                    - sorted_vertex_list.begin();
            is_deleted[max_sim_vid] = true;
        } else {
            int merge_index = -1;
            for (int i = 0; i < merge.size(); ++i) {
                if(std::find(merge[i].begin(), merge[i].end(), max_similarity_vertex) != merge[i].end()) {
                    merge_index = i;
                    break;
                }
            }
            if(merge_index == -1) {
                VertexList t = {v, max_similarity_vertex};
                merge.emplace_back(t);
                auto max_sim_vid = std::find(sorted_vertex_list.begin(), sorted_vertex_list.end(), max_similarity_vertex)
                                   - sorted_vertex_list.begin();
                is_deleted[max_sim_vid] = true;
            } else {
                merge[merge_index].emplace_back(v);
            }
        }

    }
    return merge;
}

CommunitySplit DynComDetector::merge_closure_ratio(CommunitySplit &groups, float merge_ratio) {
    Graph &g = this->graph_snapshots[0];
    std::vector<float> closure_ratio;
    closure_ratio.reserve(groups.size());
    for (const auto & group : groups) {
        closure_ratio.emplace_back(Graph::ratio(g, group));
    }
    auto min_ratio_idx = std::min_element(closure_ratio.begin(), closure_ratio.end());
    float min_ratio = *min_ratio_idx;

    while(min_ratio < merge_ratio) {
        float max_similarity = 0;
        int max_index = -1;
        long min_groups_index = min_ratio_idx - closure_ratio.begin();
        for (int j = 0; j < groups.size(); ++j) {
            if(j != min_groups_index) {
                float total_similarity = Graph::group_similarity(
                        g, groups[min_groups_index], groups[j]
                        );
                if(total_similarity > max_similarity) {
                    max_similarity = total_similarity;
                    max_index = j;
                }
            }
        }
        groups[max_index].insert(groups[max_index].end(), groups[min_groups_index].begin(),
                                 groups[min_groups_index].end());
        closure_ratio[max_index] = Graph::ratio(g, groups[max_index]);
        closure_ratio.erase(closure_ratio.begin() + min_groups_index);
        groups.erase(groups.begin() + min_groups_index);
        min_ratio_idx = std::min_element(closure_ratio.begin(), closure_ratio.end());
        min_ratio = *min_ratio_idx;
    }
    return groups;
}

float DynComDetector::compute_modularity(Graph &g, const CommunitySplit &partition) {
    unsigned int m = g.get_edge_num();
    float q = 0;
    for(const auto & group : partition) {
        unsigned int inner_edge_num = 0, outer_edge_num = 0;
        for(Vertex u : group) {
            for(Vertex u_neighbor : g.get_neighbor_list(u)) {
                if(std::find(group.begin(), group.end(), u_neighbor) != group.end()) {
                    inner_edge_num += 1;
                } else {
                    outer_edge_num += 1;
                }
            }
        }
        float eii = (float)inner_edge_num / (float)(2 * m);
        float ai = (float)(outer_edge_num + inner_edge_num) / (float)(2 * m);
        q += eii - ai * ai;
    }
    return q;
}

CommunitySplit DynComDetector::get_community_snapshot(TimeStep t) {
    if(t >= this->communities.size()) {
        std::cerr << "Snapshot not exist!" << std::endl;
        return {};
    }
    return this->communities[t];
}

DynComDetector::DynComDetector() {
    this->communities = CommunitySplitSnapshots();
    this->graph_snapshots = std::vector<Graph>(0);
    this->time_steps = TimeSteps();
}

DynComDetector::DynComDetector(const std::string& snapshot_path, unsigned int time_step_num) {
    this->communities = CommunitySplitSnapshots();
    this->graph_snapshots = std::vector<Graph>(0);
    this->time_steps = TimeSteps();
    for(TimeStep i = 1; i <= time_step_num; i++) {
        Graph g(snapshot_path + std::to_string(i));
        this->time_steps.emplace_back(i - 1);
        this->graph_snapshots.emplace_back(g);

    }
}

Graph DynComDetector::get_graph_snapshot(TimeStep t) {
    if (t >= this->graph_snapshots.size()){
        std::cerr << "timestep not exist" << std::endl;
        return Graph();
    }

    return this->graph_snapshots[t];
}

CommunitySplit DynComDetector::compute_inherited_communities(TimeStep current_step) {
    CommunitySplit previous_community_split = this->communities[current_step - 1];
    CommunitySplit inherited_community_split;
    for(const auto & group : previous_community_split) {
        Community inherited_group;
        for(Vertex v : group) {
            if(std::find(this->disappeared_vertices[current_step - 1].begin(),
                         this->disappeared_vertices[current_step - 1].end(), v) ==
                         this->disappeared_vertices[current_step - 1].end()) {
                inherited_group.emplace_back(v);
            }
        }
        if(!inherited_group.empty()) {
            inherited_community_split.emplace_back(inherited_group);
        }
    }
    return inherited_community_split;
}

CommunitySplit DynComDetector::compute_community_cores(TimeStep current_step, CommunitySplit& inherited_communities) {
    Graph &g = this->graph_snapshots[current_step];
    CommunitySplit comm_cores;
    std::map<Vertex, int> dic_inherit;
    for(int i = 0; i < inherited_communities.size(); i++) {
        for(Vertex u : inherited_communities[i]){
            dic_inherit[u] = i;
        }
    }
    for(auto & group : inherited_communities) {
        std::map<Vertex, int> dic_v_edge;
        for(Vertex v : group) {
            dic_v_edge[v] = edge_differ(g, group, v);
        }

        // Find max value of dic_v_edge;
        int max_differ = INT_MIN;
        unsigned int max_vertex = 0;
        for(const auto & p : dic_v_edge) {
            if(max_differ < p.second) {
                max_differ = p.second;
                max_vertex = p.first;
            }
        }
        while(max_differ >= 0) {
            group.erase(std::remove(group.begin(), group.end(), max_vertex), group.end());
            if (group.empty()) {
                break;
            }
            dic_v_edge.erase(max_vertex);

            int max_vertex_belong = dic_inherit[max_vertex];
            for(Vertex u : g.get_neighbor_list(max_vertex)) {
                if(dic_inherit.find(u) != dic_inherit.end() && dic_inherit[u] == max_vertex_belong) {
                    dic_v_edge[u] = edge_differ(g, inherited_communities[max_vertex_belong], u);
                } else {
                    continue;
                }
            }
            dic_inherit[max_vertex] = -1;

            max_differ = INT_MIN;
            for(const auto & p : dic_v_edge) {
                if(max_differ < p.second) {
                    max_differ = p.second;
                    max_vertex = p.first;
                }
            }
        }
    }

    for(auto & group : inherited_communities) {
        if(!group.empty()) {
            comm_cores.emplace_back(group);
        }
    }
    return comm_cores;
}

int DynComDetector::edge_differ(Graph &g, Community &c, Vertex u) {
    int in_edge_num = 0;
    int out_edge_num = 0;
    int result_edge_differ;

    VertexSet u_neighbor_set = g.get_neighbor_set(u);
    if(u_neighbor_set.empty()) {
        return -1;
    }

    for(Vertex v : u_neighbor_set) {
        if(std::find(c.begin(), c.end(), v) != c.end()) {
            in_edge_num++;
        } else {
            out_edge_num++;
        }
    }
    result_edge_differ = out_edge_num - in_edge_num;
    return result_edge_differ;
}

CommunitySplit DynComDetector::detect_communities_snapshot(TimeStep current_step) {
    CommunitySplit inherited_comm = compute_inherited_communities(current_step);
    CommunitySplit comm_cores = compute_community_cores(current_step, inherited_comm);
    CommunitySplit initial_comm = initial_merge_snapshot(current_step, comm_cores);
    CommunitySplit result_comm = final_merge_snapshot(current_step, initial_comm);
    this->communities.emplace_back(result_comm);
    return result_comm;
}

CommunitySplit DynComDetector::initial_merge_snapshot(TimeStep current_step, CommunitySplit & comm_cores) {
    Graph &g = this->graph_snapshots[current_step];
    VertexList sorted_remain_list = g.get_vertex_list();
    std::sort(sorted_remain_list.begin(), sorted_remain_list.end(), [&](Vertex u1, Vertex u2) {
        return g.get_degree(u1) > g.get_degree(u2);
    });

    for(const auto & group : comm_cores) {
        for(Vertex u : group) {
            if(std::find(sorted_remain_list.begin(), sorted_remain_list.end(), u) != sorted_remain_list.end()) {
                sorted_remain_list.erase(std::remove(sorted_remain_list.begin(),
                                                     sorted_remain_list.end(), u), sorted_remain_list.end());
            }
        }
    }

    for(Vertex v : sorted_remain_list) {
        float max_similarity = 0;
        int max_vertex = -1;

        for(Vertex u : g.get_neighbor_list(v)) {
            float similarity = Graph::compute_similarity(g, v, u);

            if(similarity > max_similarity) {
                max_similarity = similarity;
                max_vertex = (int) u;
            }
        }

        if(max_vertex == -1) {
            int max_degree = 0;
            for(Vertex u : g.get_neighbor_list(v)) {
                if(max_degree < g.get_degree(u)) {
                    max_degree  = (int) g.get_degree(u);
                }
            }
        } else {
            int merge_index = -1;
            for (int i = 0; i < comm_cores.size(); ++i) {
                if(std::find(comm_cores[i].begin(), comm_cores[i].end(), max_vertex) != comm_cores[i].end()) {
                    merge_index = i;
                    break;
                }
            }
            if(merge_index == -1) {
                VertexList t = {v, (Vertex) (max_vertex)};
                comm_cores.emplace_back(t);
                sorted_remain_list.erase(std::remove(sorted_remain_list.begin(),
                                                     sorted_remain_list.end(), max_vertex), sorted_remain_list.end());
            } else {
                comm_cores[merge_index].emplace_back(v);
            }
        }
    }
    return comm_cores;

}

CommunitySplit DynComDetector::final_merge_snapshot(TimeStep current_step, CommunitySplit &final_community) {
    Graph &g = this->graph_snapshots[current_step];
    std::map<unsigned int, Community> cs;
    for (int i = 0; i < final_community.size(); ++i) {
        cs[i] = final_community[i];
    }

    std::map<Vertex, unsigned int> vertex_to_community;
    for (int i = 0; i < final_community.size(); ++i) {
        for(Vertex v : final_community[i]) {
            vertex_to_community[v] = i;
        }
    }

    std::map<std::pair<unsigned int, unsigned int>, float> e;
    std::map<unsigned int, float> a;
    for (int i = 0; i < final_community.size(); ++i) {
        e[std::make_pair(i, i)] = 0;
        a[i] = 0;
    }
    unsigned int m = g.get_edge_num();
    for(Edge edge : g.get_edge_list()) {
        unsigned int comm_i = vertex_to_community[edge.first];
        unsigned int comm_j = vertex_to_community[edge.second];

        if(e.find(std::make_pair(comm_i, comm_j)) != e.end()) {
            e[std::make_pair(comm_i, comm_j)] = e[std::make_pair(comm_i, comm_j)] + 1;
        } else {
            e[std::make_pair(comm_i, comm_j)] = 1;
        }

        if(e.find(std::make_pair(comm_j, comm_i)) != e.end()) {
            e[std::make_pair(comm_j, comm_i)] = e[std::make_pair(comm_j, comm_i)] + 1;
        } else {
            e[std::make_pair(comm_j, comm_i)] = 1;
        }
    }

    for(auto & p : e) {
        p.second = p.second / (float)(2 * m);
        if(a.find(p.first.first) != a.end()) {
            a[p.first.first] = a[p.first.first] + e[p.first];
        } else {
            a[p.first.first] = 0;
        }
    }

    while(true) {
        float dq;
        float dq_max = -100;
        int si = -1, sj = -1;
        for(auto & p : e) {
            if(p.first.first >= p.first.second) {
                continue;
            }
            dq = e[p.first] - a[p.first.first] * a[p.first.second];
            if(dq > dq_max && dq > 0) {
                si = (int) p.first.first;
                sj = (int) p.first.second;
                dq_max = dq;
            }
        }
        if(si == -1 || sj == -1) {
            break;
        }

        // Union si and sj;
        for(Vertex u : cs[sj]) {
            if(std::find(cs[si].begin(), cs[si].end(), u) == cs[si].end()) {
                cs[si].emplace_back(u);
            }
        }
        cs.erase(sj);
        dq += dq_max * 2;
        if(dq >= dq_max) {
            dq_max = dq;
        }
        a[si] = a[si] + a[sj];
        a.erase(sj);
        e[std::make_pair(si, si)] += e[std::make_pair(sj, sj)];

        for(const auto & kv : cs) {
            auto e_jk_iter = e.find(std::make_pair(sj, kv.first));
            if(e_jk_iter != e.end()) {
                if(e.find(std::make_pair(si, kv.first)) != e.end()) {
                    e[std::make_pair(si, kv.first)] += e_jk_iter->second;
                } else {
                    e[std::make_pair(si, kv.first)] = e_jk_iter->second;
                }
            }

            auto e_kj_iter = e.find(std::make_pair(sj, kv.first));
            if(e_kj_iter != e.end()) {
                if(e.find(std::make_pair(kv.first, si)) != e.end()) {
                    e[std::make_pair(kv.first, si)] += e_kj_iter->second;
                } else {
                    e[std::make_pair(kv.first, si)] = e_kj_iter->second;
                }
            }

            if(e.find(std::make_pair(kv.first, sj)) != e.end()) {
                e.erase(std::make_pair(kv.first, sj));
            }
            if(e.find(std::make_pair(sj, kv.first)) != e.end()) {
                e.erase(std::make_pair(sj, kv.first));
            }
        }
    }
    // Process CS
    CommunitySplit result;
    for(const auto & kv : cs) {
        Community c;
        for(Vertex v : kv.second) {
            c.emplace_back(v);
        }
        result.emplace_back(c);
    }
    return result;
}

void DynComDetector::perform_dynamic_detection() {
    for (int i = 1; i < this->time_steps.size(); ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        this->detect_communities_snapshot(i);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // 将时间添加到vector容器中
        this->time_durations.emplace_back(duration);
    }
}

void DynComDetector::evaluate_modularity() {
    TimeStep t = 0;
    for(auto & cs : this->communities) {
        std::cout << "=====Time Step: " << t << "=====" << std::endl;
        std::cout << "# Communities: " << cs.size() << std::endl;
        std::cout << "Modularity: " << compute_modularity(
                this->graph_snapshots[t], cs) << std::endl;
        std::cout << "Time Consumption(Milliseconds): " << this->time_durations[t].count() << std::endl;
        ++t;
        std::cout << std::endl;
    }
}

void DynComDetector::write_result(const std::string & result_path) {
    unsigned int comm_id = 1;
    for(const auto & cs : this->communities) {
        std::ofstream comm_writer(result_path + "Community_" + std::to_string(comm_id++) + ".txt");
        if(!comm_writer.is_open()) {
            std::cerr << "Result file open failed." << std::endl;
            return;
        }

        // For each community, write a line.
        for(const Community & c : cs) {
            // For each vertex.
            for(Vertex v : c) {
                comm_writer << ++v << " ";
            }
            comm_writer << std::endl;
        }
    }
}
