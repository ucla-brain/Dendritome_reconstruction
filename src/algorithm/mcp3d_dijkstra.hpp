//
// Created by muyezhu on 9/19/19.
//

#ifndef MCP3D_MCP3D_DIJKSTRA_HPP
#define MCP3D_MCP3D_DIJKSTRA_HPP

#include <cstdint>
#include <limits>
#include <utility>
#include <unordered_set>
#include "common/mcp3d_macros.hpp"
#include "mcp3d_vertex.hpp"
#include "mcp3d_heap.hpp"
#include "mcp3d_graph.hpp"


namespace mcp3d
{
/// MDijkstra does not own the graph
/// it is assumed there is at most one edge between each ordered pair of vertices
/// currently source_id is assumed to be different from any target_ids
/// Graph type should implement:
/// NeighborVertexIds(vertex_id), Dsrc(vertex_id), SetDsrc(vertex_id, dsrc),
/// SetParentVertexId(vertex_id, parent_id)
template <typename Graph>
class MDijkstra
{
public:
    using GraphType = Graph;

    MDijkstra(): heap_(MMinHeap<typename Graph::VertexType>{}),
                 boundary_(std::unordered_map<int64_t, int64_t>{}),
                 target_ids_(std::unordered_set<int64_t>{}),
                 unreached_ids_(std::unordered_set<int64_t>{}),
                 source_id_(std::numeric_limits<int64_t>::lowest()),
                 track_boundary_(false), find_path_complete_(false)  {}

    /// if target_ids is not empty, terminates when all reachable targets are reached
    /// otherwise, shortest paths to all graph vertices are found. source_id is
    /// removed from target_ids if present
    bool FindPaths(Graph& graph, int64_t source_id,
                   const std::unordered_set<int64_t>& target_ids = std::unordered_set<int64_t>{},
                   bool track_boundary = false);

    void Clear();

    /// target_ids - unreached_ids. not same as the set of all vertices marked
    /// by the instance
    std::unordered_set<int64_t> ReachedTargetIds() const;

    int64_t source_id() const  { return source_id_; };

    const std::unordered_set<int64_t>& target_ids() const  { return target_ids_; }

    /// the id of vertices in target_ids set, to which the source vertex has no
    /// paths to
    const std::unordered_set<int64_t>& unreached_ids() const  { return unreached_ids_; }

    /// in the edge case where source vertex can not reach any vertices in graph,
    /// the boundary will contain only the source vertex
    std::unordered_map<int64_t, int64_t>& boundary() { return boundary_; };

    bool track_boundary() const  { return track_boundary_; }

    bool find_path_complete() const  { return find_path_complete_; }
private:
    void Init(Graph& graph, int64_t source_id,
              const std::unordered_set<int64_t>& target_ids = std::unordered_set<int64_t>{},
              bool track_frontier = false);

    void UpdateBoundary(int64_t vertex_id, const Graph &graph);

    MMinHeap<typename Graph::VertexType> heap_;
    // vertex_id: componenent_id
    std::unordered_map<int64_t, int64_t> boundary_;
    std::unordered_set<int64_t> target_ids_, unreached_ids_;
    int64_t source_id_;
    bool track_boundary_, find_path_complete_;

};
}

template <typename Graph>
bool mcp3d::MDijkstra<Graph>::FindPaths(Graph& graph, int64_t source_id,
                                        const std::unordered_set<int64_t>& target_ids,
                                        bool track_boundary)
{
    Init(graph, source_id, target_ids, track_boundary);
    typename Graph::EdgeType edge_cost, dsrc_candidate, dsrc_top;
    std::unique_ptr<int64_t []> non_frozen_neighbor_ids;
    int64_t n_non_frozen_neighbors, non_frozen_neighbor_id;
    INIT_TIMER(0)
    INIT_TIMER(1)
    INIT_TIMER(2)
    INIT_TIMER(3)
    INIT_TIMER(4)
    INIT_TIMER(5)
    INIT_TIMER(6)
    TIC(0)
    double time_heap_pop = 0.0, time_heap_update = 0.0, time_find_neighbor = 0.0,
            time_set_neighbor = 0.0, time_unreached_erase = 0.0, time_distance_calc = 0.0;

    while (!heap_.Empty())
    {
        TIC(1)
        int64_t top_id = heap_.Pop();
        TOC(1)
        time_heap_pop += ELAPSE(1);
        graph.SetIsFrozen(top_id, true);
        TIC(5)
        unreached_ids_.erase(top_id);
        TOC(5)
        time_unreached_erase += ELAPSE(5);
        if (track_boundary_)
            UpdateBoundary(top_id, graph);
        if (unreached_ids_.empty())
            break;
        // retrieve neighbors
        dsrc_top = graph.Dsrc(top_id);
        TIC(3)
        graph.NonFrozenNeighborVertexIds(top_id, non_frozen_neighbor_ids);
        TOC(3)
        time_find_neighbor += ELAPSE(3);
        auto non_frozen_neighbor_ids_ptr = non_frozen_neighbor_ids.get();
        n_non_frozen_neighbors = non_frozen_neighbor_ids_ptr[0];
        for (int64_t i = 1; i <= n_non_frozen_neighbors; ++i)
        {
            non_frozen_neighbor_id = non_frozen_neighbor_ids_ptr[i];
            // retrieve needed edge costs. some graph types compute edge costs
            // on the fly, so only do so if needed
            TIC(6)
            edge_cost = graph.EdgeCost(top_id, non_frozen_neighbor_id);
            dsrc_candidate = dsrc_top + edge_cost;
            TOC(6)
            time_distance_calc += ELAPSE(6);
            if (edge_cost < 0)
                MCP3D_RUNTIME_ERROR("dijkstra algorithm requires non negative weights")
            if (dsrc_candidate < graph.Dsrc(non_frozen_neighbor_id))
            {
                TIC(4)
                graph.SetDsrc(non_frozen_neighbor_id, dsrc_candidate);
                graph.SetParentVertexId(non_frozen_neighbor_id, top_id);
                TOC(4)
                time_set_neighbor += ELAPSE(4);
                TIC(2)
                heap_.Update(&(graph[non_frozen_neighbor_id]));
                TOC(2)
                time_heap_update += ELAPSE(2);
            }
        }
    }
    TOC(0)
    std::cout << "heap pop time = " << time_heap_pop << " seconds" << std::endl;
    std::cout << "heap update time = " << time_heap_update << " seconds" << std::endl;
    std::cout << "find neighbor time = " << time_find_neighbor << " seconds" << std::endl;
    std::cout << "set neighbor time = " << time_set_neighbor << " seconds" << std::endl;
    std::cout << "set unreached erase = " << time_unreached_erase << " seconds" << std::endl;
    std::cout << "set distance cal = " << time_distance_calc << " seconds" << std::endl;
    REPORT_TIME_TO_COMPLETION("dijkstra findpaths", 0)

    find_path_complete_ = true;
    return unreached_ids_.empty();
}


template <typename Graph>
std::unordered_set<int64_t> mcp3d::MDijkstra<Graph>::ReachedTargetIds() const
{
    if (!find_path_complete_)
        return std::unordered_set<int64_t> {};
    std::unordered_set<int64_t> reached_ids(target_ids_);
    for (const auto& unreachable_id: unreached_ids_)
        reached_ids.erase(unreachable_id);
    return reached_ids;
}

template <typename Graph>
void mcp3d::MDijkstra<Graph>::Init(Graph& graph, int64_t source_id,
                                   const std::unordered_set<int64_t>& target_ids,
                                   bool track_frontier)
{
    Clear();
    if (!graph.HasVertex(source_id))
        MCP3D_RUNTIME_ERROR("source vertex id " + std::to_string(source_id) + " is not in input graph")
    if (!target_ids.empty())
    {
        for (int64_t target_id: target_ids)
            if (!graph.HasVertex(target_id))
                MCP3D_RUNTIME_ERROR("target vertex id " + std::to_string(target_id) + " is not in input graph")
    }
    source_id_ = source_id;
    graph.SetDsrc(source_id_, (typename Graph::EdgeType)0);
    target_ids_ = target_ids;
    // if target_ids_ is empty, find paths to all vertices
    if (target_ids_.empty())
        for (const auto& item: graph.Vertices())
            target_ids_.insert(item.first);
    // the source vertex is not allowed to be a target
    target_ids_.erase(source_id_);
    unreached_ids_ = target_ids_;
    track_boundary_ = track_frontier;
    // initialize heap
    heap_.Push(&(graph[source_id_]));
}

template <typename Graph>
void mcp3d::MDijkstra<Graph>::Clear()
{
    source_id_ = std::numeric_limits<int64_t>::lowest();
    target_ids_.clear();
    unreached_ids_.clear();
    boundary_.clear();
    heap_.Clear();
    track_boundary_ = false;
    find_path_complete_ = false;
}

template <typename Graph>
void mcp3d::MDijkstra<Graph>::UpdateBoundary(int64_t vertex_id, const Graph &graph)
{
    MCP3D_ASSERT(graph.HasVertex(vertex_id))
    // vertex_id not in frontier
    MCP3D_ASSERT(boundary_.find(vertex_id) == boundary_.end())
    boundary_.insert(std::pair<int64_t, int64_t>(vertex_id, std::numeric_limits<int64_t>::lowest()));
    // remove parent from frontier
    boundary_.erase(graph.ParentVertexId(vertex_id));
}


#endif //MCP3D_MCP3D_DIJKSTRA_HPP
