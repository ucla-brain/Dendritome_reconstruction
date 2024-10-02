//
// Created by muyezhu on 9/27/19.
//
#include <limits>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>
#include "algorithm/mcp3d_heap.hpp"

using namespace std;

class TestElement
{
public:
    TestElement(int64_t element_id, int value): id_(element_id), value_(value) {}

    int64_t id() const  { return id_; }

    void set_value(int value)  { value_ = value; }

    bool operator< (const TestElement& other) const { return value_ < other.value_; }
private:
    int64_t id_;
    int value_;
};


class TestElementBag
{
public:
    TestElementBag() = default;

    explicit TestElementBag(int n);

    size_t Size() const { return random_values_.size(); }

    unordered_map<int64_t, TestElement>& bag()  { return bag_; };

    const unordered_set<int>& random_values() const  { return random_values_; }

private:
    void GenRandomValues(int n);
    unordered_set<int> random_values_;
    unordered_map<int64_t, TestElement> bag_;
};

TestElementBag::TestElementBag(int n)
{
    EXPECT_TRUE(n > 0);
    GenRandomValues(n);
    for (auto value: random_values_)
        bag_.emplace(value, TestElement(value, value));
}

void TestElementBag::GenRandomValues(int n)
{
    default_random_engine generator;
    uniform_int_distribution<int> distribution(0, numeric_limits<int>::max());
    int value;
    for (int i = 0; i < n; ++i)
        while (true)
        {
            value = distribution(generator);
            if (random_values_.find(value) == random_values_.end())
            {
                random_values_.insert(value);
                break;
            }
        }
}

TEST(MMinHeap, Push)
{
    int n = 100000;
    TestElementBag bag(n);
    mcp3d::MMinHeap<TestElement> heap;
    EXPECT_TRUE(heap.Empty());
    int c = 0;
    for (const auto& item: bag.bag())
    {
        heap.Push(&(item.second));
        EXPECT_TRUE(heap.Verify());

    }
    EXPECT_EQ(bag.Size(), heap.size());
    EXPECT_FALSE(heap.Empty());
    int min_value = *min_element(bag.random_values().cbegin(), bag.random_values().cend());
    // root id should be equal to minimum value of the set
    // (test element's value field is equal to id field)
    EXPECT_EQ(min_value, heap.RootId());
    // Clear() should restore empty heap
    heap.Clear();
    EXPECT_TRUE(heap.Empty());
    for (const auto& item: bag.bag())
    {
        heap.Push(&(item.second));
        EXPECT_TRUE(heap.Verify());
    }
    EXPECT_EQ(bag.Size(), heap.size());
}

TEST(MMinHeap, Pop)
{
    int n = 100000;
    TestElementBag bag(n);
    mcp3d::MMinHeap<TestElement> heap;
    for (const auto& item: bag.bag())
        heap.Push(&(item.second));
    vector<int> value_vec(bag.random_values().cbegin(), bag.random_values().cend());
    sort(value_vec.begin(), value_vec.end());
    for (int value: value_vec)
        EXPECT_EQ(value, heap.Pop());
    EXPECT_EQ((size_t)0, heap.size());
    EXPECT_TRUE(heap.Empty());
}

TEST(MMinHeap, Update)
{
    int n = 100000;
    TestElementBag bag(n);
    mcp3d::MMinHeap<TestElement> heap;
    for (const auto& item: bag.bag())
    {
        heap.Update(&(item.second));
        EXPECT_TRUE(heap.Verify());
    }
    EXPECT_EQ(bag.Size(), heap.size());
    vector<int> value_vec(bag.random_values().cbegin(), bag.random_values().cend());
    sort(value_vec.begin(), value_vec.end());
    for (const auto& value: bag.random_values())
    {
        bag.bag().at(value).set_value(-value);
        heap.Update(&(bag.bag().at(value)));
    }
    for (auto it = value_vec.rbegin(); it != value_vec.rend(); ++it)
        EXPECT_EQ(*it, heap.Pop());
    EXPECT_TRUE(heap.Empty());
}

