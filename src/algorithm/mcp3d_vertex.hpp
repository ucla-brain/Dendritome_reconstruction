//
// Created by muyezhu on 9/22/19.
//

#ifndef MCP3D_MCP3D_VERTEX_BASE_HPP
#define MCP3D_MCP3D_VERTEX_BASE_HPP

#include <cstdint>
#include <limits>
#include <type_traits>
#include <cstdlib>
#include <array>
#include <unordered_map>
#include "common/mcp3d_macros.hpp"

namespace mcp3d
{
/// most general vertex class. T is arithmetic type for edge cost
template <typename EdgeCostType = double>
class MVertexBase
{
public:
    using EdgeType = EdgeCostType;

    static_assert(std::is_arithmetic<EdgeCostType>(), "edge cost type must be floating point");
    static_assert(std::numeric_limits<EdgeCostType>::has_infinity, "edge cost type must support infinity");

    MVertexBase() = default;

    virtual bool StoreEdges() const = 0;

    virtual std::unordered_map<int64_t, EdgeCostType> edges()  { return std::unordered_map<int64_t, EdgeCostType>{}; }

    virtual int64_t id() const = 0;

    virtual void set_id(int64_t id) = 0;
};


template <typename EdgeCostType = double>
class MVertexDijkstra: public MVertexBase<EdgeCostType>
{
public:
    MVertexDijkstra(): MVertexDijkstra(-1) {}

    explicit MVertexDijkstra(int64_t vertex_id, EdgeCostType value = (EdgeCostType)0,
                             EdgeCostType dsrc = std::numeric_limits<EdgeCostType>::infinity(),
                             EdgeCostType radius = (EdgeCostType)0):
                             id_(vertex_id), value_(value), dsrc_(dsrc), radius_(radius), attributes_(0)  {}

    bool StoreEdges() const override { return false; }

    bool operator<(const MVertexDijkstra<EdgeCostType>& other) const { return dsrc_ < other.dsrc_; }

    bool operator<=(const MVertexDijkstra<EdgeCostType>& other) const { return dsrc_ <= other.dsrc_; }

    bool operator>(const MVertexDijkstra<EdgeCostType>& other) const { return dsrc_ > other.dsrc_; }

    bool operator>=(const MVertexDijkstra<EdgeCostType>& other) const { return dsrc_ >= other.dsrc_; }

    bool is_frozen() const  { return is_attribute(10); }

    bool has_parent() const  { return is_attribute(9); }

    int64_t id() const override {return id_; }

    EdgeCostType value() const  { return value_; };

    std::array<int, 3> parent_delta() const;

    EdgeCostType dsrc() const { return dsrc_; }

    EdgeCostType radius() const  { return radius_; }

    bool is_src() const { return is_attribute(15);}

    void set_id(int64_t id)  { id_ = id; }

    void set_value(EdgeCostType value)  { value_ = value; }

    void set_is_frozen(bool is_frozen)  { set_is_attribute(is_frozen, 10); }

    void set_parent_delta(int delta_z, int delta_y, int delta_x);

    void set_dsrc(EdgeCostType dsrc)  { dsrc_ = dsrc; }

    void set_radius(EdgeCostType radius)  { radius_ = radius; }

    void set_is_src(bool src)  { set_is_attribute(src, 15); }

private:
    bool is_attribute(int n_left_shift) const;

    void set_is_attribute(bool attribute, int n_left_shift);

    int64_t id_;
    // vertex value and distance to source node
    EdgeCostType value_, dsrc_, radius_;
    // 16 bit field of attributes
    // is_src_ | is_soma_ | is_leaf_ | is_severed | reserve | is_frozen | has_parent
    // lowest 9 bits: | z+1  |  z  |  z-1  |  y+1 |  y |  y-1 |  x+1 |  x |  x-1
    //                  256    128     64     32    16     8      4     2     1
    // (vertex zyx + parent_delta_ = parent vertex zyx)
    uint16_t attributes_;
};

}

template <typename EdgeCostType>
std::array<int, 3> mcp3d::MVertexDijkstra<EdgeCostType>::parent_delta() const
{
    if (!has_parent())
        return std::array<int, 3> {{-2, -2, -2}};
    int16_t x_minus_one = 1;
    // precedence: bitshift, logical comparison, bitwise logical operator
    int delta_x = ((x_minus_one & attributes_) > 0) ?
                  -1 : (((x_minus_one << 1 & attributes_) > 0) ? 0 : 1);
    int16_t y_minus_one = x_minus_one << 3;
    int delta_y = ((y_minus_one & attributes_) > 0) ?
                  -1 : (((y_minus_one << 1 & attributes_) > 0) ? 0 : 1);
    int16_t z_minus_one = x_minus_one << 6;
    int delta_z = ((z_minus_one & attributes_) > 0) ?
                  -1 : (((z_minus_one << 1 & attributes_) > 0) ? 0 : 1);
    return std::array<int, 3> {{ delta_z, delta_y, delta_x }};
}

template <typename EdgeCostType>
void mcp3d::MVertexDijkstra<EdgeCostType>::set_parent_delta(int delta_z, int delta_y, int delta_x)
{
    MCP3D_ASSERT(std::abs(delta_z) <= 1 && std::abs(delta_y) <= 1 && std::abs(delta_x) <= 1)
    // parent delta bits
    uint16_t bit_mask = 1;
    uint16_t x_bit_mask = bit_mask << (delta_x + 1),
             y_bit_mask = bit_mask << (delta_y + 4),
             z_bit_mask = bit_mask << (delta_z + 7);
    // preserve non parent delta bits, set old parent delta bits to 0s
    uint16_t high_bits = attributes_ - (attributes_ & (uint16_t)511);
    // set parent delta bits
    attributes_ = high_bits + (x_bit_mask + y_bit_mask + z_bit_mask);
    // set has_parent bit to 1
    set_is_attribute(true, 9);
}

template <typename EdgeCostType>
bool mcp3d::MVertexDijkstra<EdgeCostType>::is_attribute(int n_left_shift) const
{
    MCP3D_ASSERT(n_left_shift >= 9 && n_left_shift <= 15)
    return (attributes_ & 1 << n_left_shift) > 0;
}

template <typename EdgeCostType>
void mcp3d::MVertexDijkstra<EdgeCostType>::set_is_attribute(bool attribute,
                                                            int n_left_shift)
{
    if (is_attribute(n_left_shift) == attribute)
        return;
    if (attribute)
        attributes_ |= (1 << n_left_shift);
    else
        attributes_ -= (1 << n_left_shift);
}





#endif //MCP3D_MCP3D_VERTEX_BASE_HPP
