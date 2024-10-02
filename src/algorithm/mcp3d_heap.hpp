//
// Created by muyezhu on 9/19/19.
//

#ifndef MCP3D_MCP3D_HEAP_HPP
#define MCP3D_MCP3D_HEAP_HPP

#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include "common/mcp3d_macros.hpp"

namespace mcp3d
{
/// Element type should
/// (1) implement < operator
/// (2) has id() function that returns element id
/// the heap does not own Element instances but operates using const pointers to
/// Element instances
/// min heap. support push, pop and element updates.
/// no heapify of non heap array is implemented
/// 1 indexing. element_array_[0] is a placeholder element.
/// heap invariant: heap[k] <= heap[2k] and heap[k] <= heap[2k+1],
/// where heap[2k] and heap[2k + 1] are children of heap[k] in binary tree
template <typename Element>
class MMinHeap
{
public:
    // constructor places nullptr at elements_[0].
    // this placeholder element does not appear in element_indices_
    // capacity_: number of non place holder elements the heap can hold.
    //            elements_ is allocated capacity_ + 1
    // size_: number of non place holder elements the heap currently holds
    MMinHeap(): element_indices_(std::unordered_map<int64_t, int64_t>{}),
                elements_(std::make_unique<const Element* []>(2)), capacity_(1), size_(0)
    { elements_[0] = nullptr; }

    // increment size_
    void Push(const Element* element);
    // remove root element from heap and return its element id
    // first swap root element with last element in elements_. remove root
    // element, and swim the temp root swapped from elements_[-1] down
    // decrement size_
    int64_t Pop();
    int64_t RootId() const;
    // reposition element in the heap if pre-existing element's value field
    // has been modified. if element was not in heap, push it
    void Update(const Element* element);
    bool Empty() const  { return size_ == 0; }
    int64_t capacity() const  { return capacity_; }
    int64_t size() const  { return size_; }
    void Clear();
    const std::unique_ptr<const Element* []>& elements() const { return elements_; }
    bool Verify() const;

private:
    // when size_ is equal to capacity_, this function will double
    // capacity_ and call relloc to elements_.get(). Resize() will only increase
    // capacity_ and not decrease it
    void Resize();
    // swap elements with vid0 and vid1 in element_array
    // additionally swap their mappings in element_indices_
    void Swap(const Element* element0, const Element* element1);
    int64_t HeapIndex(int64_t element_id);
    const Element* Root();
    // parent node in the heap
    const Element* Parent(const Element* element);
    const Element* LesserChild(const Element* element);
    const Element* Last();
    void SwimUp(const Element* element);
    void SwimDown(const Element* element);
    // vertex id: index in elements_
    std::unordered_map<int64_t, int64_t> element_indices_;
    std::unique_ptr<const Element* []> elements_;
    int64_t capacity_, size_;
};

}

template <typename Element>
void mcp3d::MMinHeap<Element>::Push(const Element* element)
{
    if (element_indices_.find(element->id()) != element_indices_.end())
        MCP3D_RUNTIME_ERROR("can not push existing element to heap. use Update instead.")
    Resize();
    elements_[++size_] = element;
    element_indices_[element->id()] = size_;
    if (size_ == 1)
        return;
    SwimUp(elements_[size_]);
}

template <typename Element>
int64_t mcp3d::MMinHeap<Element>::Pop()
{
    if (size_ == 0)
        return RootId();
    int64_t root_id = RootId();
    Swap(Root(), Last());
    // remove root node index mapping
    element_indices_.erase(root_id);
    // decrement size_
    --size_;
    if (size_ > 0)
        SwimDown(elements_[1]);
    return root_id;
}

template <typename Element>
void mcp3d::MMinHeap<Element>::Update(const Element *element)
{
    if (element_indices_.find(element->id()) == element_indices_.end())
        Push(element);
    SwimUp(element);
    SwimDown(element);
}

template <typename Element>
int64_t mcp3d::MMinHeap<Element>::RootId() const
{
    if (size_ == 0)
        return std::numeric_limits<int64_t>::lowest();
    return elements_[1]->id();
}

template <typename Element>
const Element* mcp3d::MMinHeap<Element>::Root()
{
    if (size_ == 0)
        return nullptr;
    return elements_[1];
}

template <typename Element>
const Element* mcp3d::MMinHeap<Element>::Last()
{
    if (size_ == 0)
        return nullptr;
    return elements_[size_];
}

template <typename Element>
int64_t mcp3d::MMinHeap<Element>::HeapIndex(int64_t element_id)
{
    if (element_indices_.find(element_id) == element_indices_.end())
        MCP3D_RUNTIME_ERROR("vertex id " + std::to_string(element_id) + " not in heap")
    return element_indices_[element_id];
}

template <typename Element>
const Element* mcp3d::MMinHeap<Element>::Parent(const Element* element)
{
    int64_t eid = element->id();
    if (element_indices_.find(eid) == element_indices_.end())
        MCP3D_RUNTIME_ERROR("element id " + std::to_string(eid) + " not in heap")
    int64_t index = element_indices_[eid];
    // placeholder node return 0
    int64_t parent_index = index / 2;
    if (parent_index == 0)
        return nullptr;
    return elements_[parent_index];
}

template <typename Element>
const Element* mcp3d::MMinHeap<Element>::LesserChild(const Element* element)
{
    int64_t eid = element->id();
    if (element_indices_.find(eid) == element_indices_.end())
        MCP3D_RUNTIME_ERROR("vertex id " + std::to_string(eid) + " not in heap")
    int64_t index = element_indices_[eid];
    int64_t child_index0 = 2 * index;
    int64_t child_index1 = child_index0 + 1;
    if (child_index0 > size_)
        // element has no child
        return nullptr;
    else if (child_index1 > size_)
        // element has one child
        return elements_[child_index0];
    else
        // smaller of the two children
        return *(elements_[child_index0]) < *(elements_[child_index1]) ?
               elements_[child_index0] : elements_[child_index1];
}

template <typename Element>
void mcp3d::MMinHeap<Element>::Resize()
{
    if (capacity_ == size_)
    {
        // size of current elements_ in bytes
        size_t n_bytes = (size_t)(capacity_ + 1) * sizeof(const Element*);
        capacity_ *= 2;
        auto new_ptr = new const Element*[capacity_ + 1];
        // copy content of current elements_.get.
        // realloc calls frees old block if new block is allocated, this
        // will cause double free with elements_.reset
        std::memcpy(new_ptr, elements_.get(), n_bytes);
        elements_.reset(new_ptr);
    }
}

template <typename Element>
void mcp3d::MMinHeap<Element>::Swap(const Element* element0, const Element* element1)
{
    int64_t eid0 = element0->id();
    int64_t eid1 = element1->id();
    // old heap indices
    int64_t index0 = HeapIndex(eid0);
    int64_t index1 = HeapIndex(eid1);
    // swap heap elements at index0 and index1
    std::swap<const Element*>(elements_[index0], elements_[index1]);
    // swap mapping of vids in element_indices
    std::swap(element_indices_[eid0], element_indices_[eid1]);
}

template <typename Element>
void mcp3d::MMinHeap<Element>::SwimUp(const Element* element)
{
    const Element* parent = Parent(element);
    while (parent && *element < *parent)
    {
        Swap(element, parent);
        parent = Parent(element);
    }
}

template <typename Element>
void mcp3d::MMinHeap<Element>::SwimDown(const Element* element)
{
    const Element* child = LesserChild(element);
    while (child && !(*element < *child))
    {
        Swap(element, child);
        child = LesserChild(element);
    }
}

template <typename Element>
bool mcp3d::MMinHeap<Element>::Verify() const
{
    for (int64_t i = 2; i <= size_; ++i)
        if (!(*(elements_[i / 2]) < *(elements_[i])))
            return false;
    return true;
}

template <typename Element>
void mcp3d::MMinHeap<Element>::Clear()
{
    elements_ = std::make_unique<const Element* []>(2);
    elements_[0] = nullptr;
    capacity_ = 1;
    size_ = 0;
    element_indices_.clear();
}

#endif //MCP3D_MCP3D_HEAP_HPP
