//
// Created by muyezhu on 9/27/19.
//
#include <limits>
#include <array>
#include <gtest/gtest.h>
#include "common/mcp3d_exceptions.hpp"
#include "algorithm/mcp3d_vertex.hpp"

using namespace std;

TEST(MVertexAllPath, Constructor)
{
    cout << "size of mcp3d::MVertexDijkstra<float> = " << sizeof(mcp3d::MVertexDijkstra<float>) << endl;
    cout << "size of mcp3d::MVertexDijkstra<double> = " << sizeof(mcp3d::MVertexDijkstra<double>) << endl;
    mcp3d::MVertexDijkstra<double> vertex0 {};
    EXPECT_EQ(-1, vertex0.id());
    EXPECT_EQ(numeric_limits<double>::infinity(), vertex0.dsrc());
    EXPECT_TRUE(numeric_limits<double>::max() < vertex0.dsrc());
    EXPECT_FALSE(vertex0.is_src());
    EXPECT_EQ(-2, vertex0.parent_delta()[0]);
    EXPECT_EQ(-2, vertex0.parent_delta()[1]);
    EXPECT_EQ(-2, vertex0.parent_delta()[2]);
    EXPECT_FALSE(vertex0.has_parent());

    mcp3d::MVertexDijkstra<double> vertex1(3, 3.0, 4.0);
    EXPECT_EQ(3, vertex1.id());
    EXPECT_EQ(4.0, vertex1.dsrc());
    EXPECT_FALSE(vertex1.is_src());
    EXPECT_EQ(-2, vertex1.parent_delta()[0]);
    EXPECT_EQ(-2, vertex1.parent_delta()[1]);
    EXPECT_EQ(-2, vertex1.parent_delta()[2]);
    EXPECT_FALSE(vertex1.has_parent());
}

TEST(MVertexAllPath, Comparison)
{
    mcp3d::MVertexDijkstra<double> vertex0 {};
    mcp3d::MVertexDijkstra<double> vertex1(3, 3.0, 4.0);
    EXPECT_TRUE(vertex0 > vertex1);
    vertex0.set_dsrc(3.9);
    EXPECT_TRUE(vertex0 < vertex1);
}

TEST(MVertexAllPath, ParentBits)
{
    mcp3d::MVertexDijkstra<double> vertex(3, 4.0);
    EXPECT_FALSE(vertex.has_parent());
    // can correctly set and reset parent
    int delta_z = -1, delta_y = 0, delta_x = 1;
    vertex.set_parent_delta(delta_z, delta_y, delta_x);
    EXPECT_TRUE(vertex.has_parent());
    EXPECT_EQ(delta_z, vertex.parent_delta()[0]);
    EXPECT_EQ(delta_y, vertex.parent_delta()[1]);
    EXPECT_EQ(delta_x, vertex.parent_delta()[2]);
    delta_z = 2;
    EXPECT_THROW(vertex.set_parent_delta(delta_z, delta_y, delta_x), mcp3d::MCPAssertionError);
    delta_z = 0;
    delta_y = -1;
    delta_x = -1;
    vertex.set_parent_delta(delta_z, delta_y, delta_x);
    EXPECT_TRUE(vertex.has_parent());
    EXPECT_EQ(delta_z, vertex.parent_delta()[0]);
    EXPECT_EQ(delta_y, vertex.parent_delta()[1]);
    EXPECT_EQ(delta_x, vertex.parent_delta()[2]);
}

TEST(MVertexAllPath, IsFrozen)
{
    mcp3d::MVertexDijkstra<float> vertex;
    EXPECT_FALSE(vertex.is_frozen());
    EXPECT_FALSE(vertex.has_parent());
    EXPECT_FALSE(vertex.is_src());
    vertex.set_is_frozen(false);
    EXPECT_FALSE(vertex.is_frozen());
    EXPECT_FALSE(vertex.has_parent());
    EXPECT_FALSE(vertex.is_src());
    vertex.set_is_frozen(true);
    EXPECT_TRUE(vertex.is_frozen());
    EXPECT_FALSE(vertex.has_parent());
    EXPECT_FALSE(vertex.is_src());
}

TEST(MVertexAllPath, SrcBit)
{
    mcp3d::MVertexDijkstra<double> vertex(3, 4.0);
    EXPECT_FALSE(vertex.has_parent());
    EXPECT_FALSE(vertex.is_src());
    vertex.set_is_src(false);
    EXPECT_FALSE(vertex.has_parent());
    EXPECT_FALSE(vertex.is_src());
    vertex.set_is_src(true);
    EXPECT_FALSE(vertex.has_parent());
    EXPECT_TRUE(vertex.is_src());
}
