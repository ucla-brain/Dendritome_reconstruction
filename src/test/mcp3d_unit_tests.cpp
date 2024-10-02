//
// Created by muyezhu on 10/9/17.
//
#include <iostream>
#include <tiffio.h>
#include <gtest/gtest.h>


using namespace std;

int main(int argc, char **argv)
{
    ::TIFFSetWarningHandler(nullptr);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

