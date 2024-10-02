//
// Created by muyezhu on 3/26/18.
//
#include <sstream>
#include "common/mcp3d_common.hpp"
#include <gtest/gtest.h>

using namespace std;

void ThrowTwo()
{
    // the throw_with_nested can only be used in catch blocks. actual point
    // of exception needs to only throw, otherwise the inner most exception
    // is not correctly rethrown in the print nested exception function.
    throw runtime_error("throw2");
}

void ThrowOne()
{
    try
    {
        ThrowTwo();
    }
    catch (...)
    {
        throw_with_nested(runtime_error("throw1"));
    }
}

void ThrowZero()
{
    try
    {
        ThrowOne();
    }
    catch (...)
    {
        throw_with_nested(runtime_error("throw0"));
    }
}

void MCPThrowOne()
{
    MCP3D_RUNTIME_ERROR("mcp throw one")
}

void MCPThrowZero()
{
    try
    {
        MCPThrowOne();
    }
    catch (...)
    {
        MCP3D_RETHROW(current_exception())
    }
}

TEST(Exceptions, PrintExceptionImpl)
{
    stringstream expected, outcome;
    expected << "throw0\nthrow1\nthrow2\n";
    try
    {
        ThrowZero();
    }
    catch (const exception& e)
    {
        mcp3d::PrintNested(e, outcome);
    }
    EXPECT_EQ(expected.str(), outcome.str());
    try
    {
        MCPThrowZero();
    }
    catch (const mcp3d::MCPRuntimeError& e)
    {
        mcp3d::PrintNestedException(::std::current_exception(), __FILE__, __LINE__, __func__);
    }
}

void ThrowingFunction()
{
    MCP3D_RUNTIME_ERROR("throwing error");
}

TEST(Exceptions, MultiThreadedExceptions)
{
    mcp3d::MultiThreadExceptions me;
    EXPECT_TRUE(me.e_ptr() == nullptr);
    me.RunAndCaptureException(ThrowingFunction);
    EXPECT_FALSE(me.e_ptr() == nullptr);
    EXPECT_TRUE(me.HasCapturedException());
}
