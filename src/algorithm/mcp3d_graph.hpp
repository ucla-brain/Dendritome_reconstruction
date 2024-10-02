//
// Created by muyezhu on 9/26/19.
//

#ifndef MCP3D_MCP3D_GRAPH_HPP
#define MCP3D_MCP3D_GRAPH_HPP

#include <cstdint>
#include <limits>
#include <array>
#include <vector>
#include <unordered_map>
#include "common/mcp3d_macros.hpp"
#include "image/mcp3d_image_utils.hpp"
#include "image/mcp3d_image_maths.hpp"
#include "mcp3d_algorithm_macros.hpp"
#include "mcp3d_vertex.hpp"


namespace mcp3d
{
template <typename Vertex>
class MGraph
{
public:
    using VertexType = Vertex;

    using EdgeType = typename Vertex::EdgeType;

    MGraph(): vertices_(std::unordered_map<int64_t, Vertex>{}) {}

    virtual void Clear() = 0;

    virtual bool IsDirected() const = 0;

    bool StoreEdges() const { return Vertex{}.StoreEdges(); }

    size_t Size() const  { return vertices_.size(); }

    std::unordered_map<int64_t, Vertex>& Vertices()  { return vertices_; }

    const std::unordered_map<int64_t, Vertex>& Vertices() const  { return vertices_; }

    bool HasVertex(int64_t vertex_id) const
    { return vertices_.find(vertex_id) != vertices_.end(); }

    // return default constructed vertex instance if vertex_id not found
    Vertex& operator[](int64_t vertex_id)  { return vertices_[vertex_id]; }

    // throw error if vertex id not found
    Vertex& At(int64_t vertex_id);

    const Vertex& At(int64_t vertex_id) const;

    virtual std::vector<int64_t> NeighborVertexIds(int64_t vertex_id) const = 0;

    virtual EdgeType EdgeCost(int64_t vertex_id0, int64_t vertex_id1) = 0;

    virtual EdgeType EdgeCost(int64_t vertex_id0, int64_t vertex_id1) const = 0;

    virtual std::unordered_map<int64_t, EdgeType> Edges(int64_t vertex_id) = 0;

private:
    std::unordered_map<int64_t, Vertex> vertices_;
};


template <typename Vertex, typename VType>
class MGraphAllPath: public MGraph<Vertex>
{
public:
    MCP3D_MGRAPH_ALIASES

    using DataType = VType;

    MGraphAllPath(): MGraphAllPath(3) {}

    explicit MGraphAllPath(int cnn_type): MGraph<Vertex> {},
                                          data_(nullptr), data_zdim_(0), data_ydim_(0), data_xdim_(0),
                                          srcz_(-1), srcy_(-1), srcx_(-1),
                                          cnn_type_(cnn_type), n_total_voxels_(0), n_yx_voxels_(0),
                                          data_max_(std::numeric_limits<EdgeType>::lowest()) {}

    void Clear() override;

    EdgeType EdgeCost(int64_t vertex_id0, int64_t vertex_id1) override;

    EdgeType EdgeCost(int64_t vertex_id0, int64_t vertex_id1) const override;

    std::vector<int64_t> NeighborVertexIds(int64_t vertex_id) const override;

    void NonFrozenNeighborVertexIds(int64_t vertex_id, std::unique_ptr<int64_t []>& non_frozen_neighbor_ids) const;

    std::unordered_map<int64_t, EdgeType> Edges(int64_t vertex_id) override;

    // ?? ambigous. directionality does not exist for the input graph defined
    // by voxel connectivity, but is established by graph algorithm
    bool IsDirected() const override { return false; }

    void Init(void* input, int input_zdim, int input_ydim, int input_xdim, int srcz, int srcy, int srcx, int cnn_type = 3);

    void Init(void* input, int input_zdim, int input_ydim, int input_xdim, double srcz, double srcy, double srcx, int cnn_type = 3)
    { Init(input, input_zdim, input_ydim, input_ydim, (int)std::round(srcz), (int)std::round(srcy), (int)std::round(srcx), cnn_type); }

    void Init(void* input, int input_zdim, int input_ydim, int input_xdim, float srcz, float srcy, float srcx, int cnn_type = 3)
    { Init(input, input_zdim, input_ydim, input_ydim, (int)std::round(srcz), (int)std::round(srcy), (int)std::round(srcx), cnn_type); }

    int64_t SrcId() const;

    const VType* data() const  { return data_; }

    std::vector<int> data_dims() const { return { data_zdim_, data_ydim_, data_xdim_ }; }

    int cnn_type() const  { return cnn_type_; }

    int64_t n_total_voxels() const  { return n_total_voxels_; }

    bool IsFrozen(int64_t vertex_id) const  { return At(vertex_id).is_frozen(); }

    EdgeType Dsrc(int64_t vertex_id) const  { return At(vertex_id).dsrc(); }

    EdgeType NormalizedDataValue(int64_t vertex_id) const;

    // returns the vertex value, which is equal to
    // exp(lambda *  pow((1.0 - data_ptr[address] / data_max), 2.0))
    EdgeType VertexValue(int64_t vertex_id) const  { return At(vertex_id).value(); }

    bool HasParent(int64_t vertex_id) const  { return At(vertex_id).has_parent(); }

    // return numeric_limits<int64_t>::lowest() if no parent
    int64_t ParentVertexId(int64_t vertex_id) const;

    // return nullptr if vertex has no parent
    Vertex* ParentVertex(int64_t vertex_id);

    void SetIsFrozen(int64_t vertex_id, bool is_frozen)  { At(vertex_id).set_is_frozen(is_frozen); }

    void SetParentVertexId(int64_t vertex_id, int64_t parent_vertex_id);

    void SetDsrc(int64_t vertex_id, EdgeType dsrc) { At(vertex_id).set_dsrc(dsrc); }

private:
    // mapping of (z, y, x) into linear address. if (z, y, x) is out of bounds, return -1
    int64_t GetVid(int z, int y, int x) const;

    void GetZyx(int64_t vertex_id, std::array<int, 3>& zyx) const;

    /// return euclidean distance if the two vertices are connected (defined by cnn_type_)
    /// otherwise return infinity
    EdgeType DistanceEuclidean(int64_t vertex_id0, int64_t vertex_id1) const;

    VType* data()  { return data_; }

    EdgeType data_max() const { return data_max_; }

    void set_data(void* input);

    void set_data_dims(int data_zdim, int data_ydim, int data_xdim);

    void set_src(int srcz, int srcy, int srcx);

    void set_cnn_type(int cnn_type);

    // does not manage data
    VType* data_;
    int data_zdim_, data_ydim_, data_xdim_, srcz_, srcy_, srcx_, cnn_type_;
    int64_t n_total_voxels_, n_yx_voxels_;
    EdgeType data_max_;
};

}

template <typename Vertex>
Vertex& mcp3d::MGraph<Vertex>::At(int64_t vertex_id)
{
    if (!HasVertex(vertex_id))
        MCP3D_OUT_OF_RANGE("vertex id " + std::to_string(vertex_id) + " not found")
    return vertices_.at(vertex_id);
}

template <typename Vertex>
const Vertex& mcp3d::MGraph<Vertex>::At(int64_t vertex_id) const
{
    if (!HasVertex(vertex_id))
        MCP3D_OUT_OF_RANGE("vertex id " + std::to_string(vertex_id) + " not found")
    return vertices_.at(vertex_id);
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::Clear()
{
    data_ = nullptr;
    data_zdim_ = data_ydim_ = data_xdim_ = 0;
    n_total_voxels_ = 0;
    n_yx_voxels_ = 0;
    srcz_ = srcy_ = srcx_ = -1;
    cnn_type_ = 3;
    data_max_ = std::numeric_limits<EdgeType>::lowest();
    Vertices().clear();
}

template <typename Vertex, typename VType>
typename mcp3d::MGraphAllPath<Vertex, VType>::EdgeType
mcp3d::MGraphAllPath<Vertex, VType>::EdgeCost(int64_t vertex_id0, int64_t vertex_id1)
{
    EdgeType d_l2 = DistanceEuclidean(vertex_id0, vertex_id1);
    return d_l2 * (VertexValue(vertex_id0) + VertexValue(vertex_id1)) / 2;
}

template <typename Vertex, typename VType>
typename mcp3d::MGraphAllPath<Vertex, VType>::EdgeType
mcp3d::MGraphAllPath<Vertex, VType>::EdgeCost(int64_t vertex_id0, int64_t vertex_id1) const
{
    EdgeType d_l2 = DistanceEuclidean(vertex_id0, vertex_id1);
    EdgeType d_intensity = (VertexValue(vertex_id0) + VertexValue(vertex_id1)) / 2;
    return d_l2 * d_intensity;
}

/// 6 connected: (x, y, z) <-> (x+-1, y, z), (x, y+-1, z), (x, y, z+-1)
/// 18 connected: (x, y, z) <-> (x+-1, y+-1, z), (x+-1, y, z+-1), (x, y+-1, z+-1) + 6 connected
/// 26 connected: (x, y, z) <-> (x+-1, y+-1, z+-1) + 18 connected
template <typename Vertex, typename VType>
std::vector<int64_t> mcp3d::MGraphAllPath<Vertex, VType>::NeighborVertexIds(int64_t vertex_id) const
{
    MCP3D_ASSERT(HasVertex(vertex_id))
    std::vector<int64_t> neighbor_vertex_ids;
    std::unique_ptr<int64_t []> neighbor_candidates;
    mcp3d::NeighborAddresses(data_dims(), vertex_id, cnn_type_, neighbor_candidates);
    int64_t n_candidates = neighbor_candidates[0];
    for (int64_t i = 1; i <= n_candidates; ++i)
    {
        if (data()[neighbor_candidates[i]] > 0)
            neighbor_vertex_ids.push_back(neighbor_candidates[i]);
    }
    return neighbor_vertex_ids;
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::NonFrozenNeighborVertexIds(
        int64_t vertex_id, std::unique_ptr<int64_t[]> &non_frozen_neighbor_ids) const
{
    MCP3D_ASSERT(HasVertex(vertex_id))
    non_frozen_neighbor_ids = std::make_unique<int64_t []>(27);
    auto non_frozen_neighbor_ids_ptr = non_frozen_neighbor_ids.get();
    std::unique_ptr<int64_t []> neighbor_candidates;
    mcp3d::NeighborAddresses(data_dims(), vertex_id, cnn_type_, neighbor_candidates);
    auto neighbor_candidates_ptr = neighbor_candidates.get();
    int64_t n_candidates = neighbor_candidates_ptr[0], n_non_frozen_neighbors = 0, neighbor_candidate_id;
    for (int64_t i = 1; i <= n_candidates; ++i)
    {
        neighbor_candidate_id = neighbor_candidates_ptr[i];
        if (data()[neighbor_candidate_id] > 0 && !At(neighbor_candidate_id).is_frozen())
            non_frozen_neighbor_ids_ptr[++n_non_frozen_neighbors] = neighbor_candidate_id;
    }
    non_frozen_neighbor_ids_ptr[0] = n_non_frozen_neighbors;
}

template <typename Vertex, typename VType>
std::unordered_map<int64_t, typename mcp3d::MGraphAllPath<Vertex, VType>::EdgeType>
mcp3d::MGraphAllPath<Vertex, VType>::Edges(int64_t vertex_id)
{
    if (!HasVertex(vertex_id))
        MCP3D_OUT_OF_RANGE("vertex id " + std::to_string(vertex_id) + " not found")
    if (StoreEdges())
        return At(vertex_id).edges();
    std::unordered_map<int64_t, EdgeType> edges;
    for (auto neighbor_vertex_id: NeighborVertexIds(vertex_id))
        edges[neighbor_vertex_id] = EdgeCost(vertex_id, neighbor_vertex_id);
    return edges;
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::Init(void *input, int input_zdim, int input_ydim, int input_xdim,
                                               int srcz, int srcy, int srcx, int cnn_type)
{
    Clear();
    if (!input)
        MCP3D_RUNTIME_ERROR("no input data")
    set_data_dims(input_zdim, input_ydim, input_xdim);
    set_data(input);
    set_src(srcz, srcy, srcx);
    // source vertex must have positive value
    int64_t src_id = mcp3d::LinearAddress(data_dims(), srcz_, srcy_, srcx_);
    if (data()[src_id] <= 0)
        MCP3D_RUNTIME_ERROR("source voxel must be positive")
    set_cnn_type(cnn_type);
    int64_t n_total = n_total_voxels();
    for (int64_t i = 0; i < n_total; ++i)
    {
        if (data()[i] > 0)
        {
            VType data_value = data()[i];
            auto lambda = (EdgeType)10;
            EdgeType vertex_value = exp(lambda *  std::pow(((EdgeType)1.0 -
                                                            (EdgeType)data_value / (EdgeType)data_max()),
                                                           (EdgeType)2));
            Vertices().emplace(i, Vertex(i, vertex_value));
        }
    }
    operator[](SrcId()).set_is_src(true);
}

template <typename Vertex, typename VType>
int64_t mcp3d::MGraphAllPath<Vertex, VType>::GetVid(int z, int y, int x) const
{
    if (z < 0 || z >= data_zdim_)
        return -1;
    if (y < 0 || y >= data_ydim_)
        return -1;
    if (x < 0 || x >= data_xdim_)
        return -1;
    return mcp3d::LinearAddress(std::vector<int>({data_zdim_, data_ydim_, data_xdim_}), z, y, x);
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::GetZyx(int64_t vertex_id, std::array<int, 3>& zyx) const
{
    MCP3D_ASSERT(HasVertex(vertex_id))
    zyx[0] = (int)(vertex_id / n_yx_voxels_);
    zyx[1] = (int)(vertex_id % n_yx_voxels_ / data_xdim_);
    zyx[2] = (int)(vertex_id % data_xdim_);
}

template <typename Vertex, typename VType>
typename mcp3d::MGraphAllPath<Vertex, VType>::EdgeType
mcp3d::MGraphAllPath<Vertex, VType>::DistanceEuclidean(int64_t vertex_id0, int64_t vertex_id1) const
{
    MCP3D_ASSERT(HasVertex(vertex_id0) && HasVertex(vertex_id1))
    std::array<int, 3> zyx0, zyx1;
    GetZyx(vertex_id0, zyx0);
    GetZyx(vertex_id1, zyx1);
    int dz = zyx0[0] - zyx1[0], dy = zyx0[1] - zyx1[1], dx = zyx0[2] - zyx1[2];
    int dl1 = std::abs(dz) + std::abs(dy) + std::abs(dx);
    if (std::abs(dz) > 1 || std::abs(dy) > 1 || std::abs(dx) > 1)
        return std::numeric_limits<EdgeType>::infinity();
    return dl1 <= cnn_type_ ? std::sqrt((EdgeType)dl1) : std::numeric_limits<EdgeType>::infinity();
}

template <typename Vertex, typename VType>
int64_t mcp3d::MGraphAllPath<Vertex, VType>::SrcId() const
{
    if (srcx_ < 0)
        return std::numeric_limits<int64_t>::lowest();
    return GetVid(srcz_, srcy_, srcx_);
}

template <typename Vertex, typename VType>
typename mcp3d::MGraphAllPath<Vertex, VType>::EdgeType
mcp3d::MGraphAllPath<Vertex, VType>::NormalizedDataValue(int64_t vertex_id) const
{
    MCP3D_ASSERT(vertex_id >= 0 && vertex_id < n_total_voxels_)
    return data_[vertex_id] / data_max_;
}

template <typename Vertex, typename VType>
int64_t mcp3d::MGraphAllPath<Vertex, VType>::ParentVertexId(int64_t vertex_id) const
{
    if (!HasVertex(vertex_id))
        MCP3D_OUT_OF_RANGE("vertex id " + std::to_string(vertex_id) + " not found")
    if (!At(vertex_id).has_parent())
        return std::numeric_limits<int64_t>::lowest();
    std::array<int, 3> parent_vertex_delta = At(vertex_id).parent_delta(), vertex_zyx;
    GetZyx(vertex_id, vertex_zyx);
    return GetVid(vertex_zyx[0] + parent_vertex_delta[0],
                  vertex_zyx[1] + parent_vertex_delta[1],
                  vertex_zyx[2] + parent_vertex_delta[2]);
}

template <typename Vertex, typename VType>
Vertex* mcp3d::MGraphAllPath<Vertex, VType>::ParentVertex(int64_t vertex_id)
{
    int64_t parent_vertex_id = ParentVertexId(vertex_id);
    if (parent_vertex_id < 0)
        return nullptr;
    return &At(parent_vertex_id);
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::SetParentVertexId(int64_t vertex_id,
                                                            int64_t parent_vertex_id)
{
    MCP3D_ASSERT(HasVertex(vertex_id) && HasVertex(parent_vertex_id))
    std::array<int, 3> vertex_zyx, parent_vertex_zyx;
    GetZyx(vertex_id, vertex_zyx);
    GetZyx(parent_vertex_id, parent_vertex_zyx);
    operator[](vertex_id).set_parent_delta(parent_vertex_zyx[0] - vertex_zyx[0],
                                           parent_vertex_zyx[1] - vertex_zyx[1],
                                           parent_vertex_zyx[2] - vertex_zyx[2]);
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::set_data(void* input)
{
    MCP3D_ASSERT(input)
    data_ = (VType*)input;
    data_max_ = (EdgeType)mcp3d::VolumeMax<VType>(data_, std::vector<int>({data_zdim_, data_ydim_, data_xdim_}));
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::set_data_dims(int data_zdim,
                                                        int data_ydim,
                                                        int data_xdim)
{
    MCP3D_ASSERT(data_zdim > 0 && data_ydim > 0 && data_xdim > 0)
    data_zdim_ = data_zdim;
    data_ydim_ = data_ydim;
    data_xdim_ = data_xdim;
    n_yx_voxels_ = (int64_t)data_zdim_ * (int64_t)data_ydim_;
    n_total_voxels_ = n_yx_voxels_ * (int64_t)data_xdim_;
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::set_src(int srcz, int srcy, int srcx)
{
    MCP3D_ASSERT(srcz >= 0 && srcz < data_zdim_)
    MCP3D_ASSERT(srcy >= 0 && srcy < data_ydim_)
    MCP3D_ASSERT(srcx >= 0 && srcz < data_xdim_)
    srcz_ = srcz;
    srcy_ = srcy;
    srcx_ = srcx;
}

template <typename Vertex, typename VType>
void mcp3d::MGraphAllPath<Vertex, VType>::set_cnn_type(int cnn_type)
{
    MCP3D_ASSERT(cnn_type == 1 || cnn_type == 2 || cnn_type == 3)
    cnn_type_ = cnn_type;
}


#endif //MCP3D_MCP3D_GRAPH_HPP
