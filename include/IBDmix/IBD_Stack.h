#ifndef INCLUDE_IBDMIX_IBD_STACK_H_
#define INCLUDE_IBDMIX_IBD_STACK_H_

#include <vector>
#include <iostream>

struct IBD_Node {
    double cumulative_lod, lod;
    uint64_t position;
    unsigned char bitmask;
    IBD_Node *next;
};

struct IBD_Stack {
    IBD_Node *start = nullptr, *end = nullptr, *top = nullptr;
    IBD_Stack() = default;
    explicit IBD_Stack(IBD_Node *top) : top(top) {}

    void push(IBD_Node * new_node);
    IBD_Node * pop();
    int size();
    void write(std::ostream &strm) const;
    void reverse();
    bool empty();
};

std::ostream& operator<<(std::ostream &strm, const IBD_Stack &stack);

class IBD_Pool {
 private:
    int buffer_size;
    IBD_Stack pool;
    std::vector<IBD_Node*> alloc_ptrs;

    void allocate(int nodes);

 public:
    explicit IBD_Pool(int initial_buffer = 1024);
    ~IBD_Pool();

    IBD_Node* get_node(uint64_t position, double lod = 0,
            unsigned char bitmask = 0);
    int size();
    void reclaim_node(IBD_Node *node);
    void reclaim_after(IBD_Node *start);
    void reclaim_between(IBD_Node *start, IBD_Node *end);
    void reclaim_all(IBD_Node **top);
};

#endif  // INCLUDE_IBDMIX_IBD_STACK_H_
