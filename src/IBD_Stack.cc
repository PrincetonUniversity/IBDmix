#include "IBDmix/IBD_Stack.h"

void IBD_Stack::push(IBD_Node *new_node){
    new_node->next = top;
    top = new_node;
}

IBD_Node* IBD_Stack::pop(){
    IBD_Node *result = top;
    top = top-> next;
    return result;
}

int IBD_Stack::size(){
    int count = 0;
    for(IBD_Node* ptr = top; ptr != nullptr; ptr = ptr->next)
        count++;
    return count;
}

void IBD_Stack::write(std::ostream &strm) const{
    for(IBD_Node* ptr = top; ptr != nullptr; ptr = ptr->next)
        strm << ptr-> position << " ";
    strm << '\n';
}

std::ostream& operator<<(std::ostream &strm, const IBD_Stack &stack){
    stack.write(strm);
    return strm;
}

void IBD_Stack::reverse(){
    // reverse the order of the list
    // assumes nothing before top
    // performed in place!
    IBD_Node *result = nullptr, *temp;
    while(top != nullptr){
        // push top onto result
        temp = pop();
        temp->next = result;
        result = temp;
    }

    top = result;
}

bool IBD_Stack::empty(){
    return top == nullptr;
}

IBD_Pool::IBD_Pool(int initial_buffer) : buffer_size(initial_buffer) {
    allocate(buffer_size);
}

void IBD_Pool::allocate(int nodes){
    // assumes only called when pool is nullptr
    IBD_Node* allocation = (IBD_Node*) malloc(sizeof(IBD_Node) * nodes);
    pool.top = allocation;
    for (int i = 0; i < nodes; i++){
        allocation[i].next = &allocation[i+1];
    }
    allocation[nodes-1].next = nullptr;
    alloc_ptrs.push_back(allocation);
}

IBD_Pool::~IBD_Pool(){
    for (auto &ptr : alloc_ptrs)
        free(ptr);
    alloc_ptrs.clear();
}

IBD_Node* IBD_Pool::get_node(unsigned long position, double lod, unsigned char bitmask){
    if (pool.empty()){
        allocate(buffer_size);
        buffer_size <<= 1;  // double next request
    }

    IBD_Node* result = pool.pop();

    result->position = position;
    result->lod = lod;
    result->bitmask = bitmask;
    result->cumulative_lod = lod;
    result->next = nullptr;
    return result;
}

void IBD_Pool::reclaim_node(IBD_Node *node){
    pool.push(node);
}

void IBD_Pool::reclaim_after(IBD_Node *start){
    reclaim_between(start, nullptr);
}

void IBD_Pool::reclaim_between(IBD_Node *start,
        IBD_Node *end){
    if(start == nullptr)
        return;
    // move all nodes after start into pool
    IBD_Node* ptr = start;
    //find last node
    for (; ptr->next != end; ptr = ptr->next)
        ;
    ptr->next = pool.top;
    pool.top = start->next;
    start->next = end;
}

void IBD_Pool::reclaim_all(IBD_Node *&top){
    if(top == nullptr)
        return;
    reclaim_after(top);
    top->next = pool.top;
    pool.top = top;
    top = nullptr;
}

int IBD_Pool::size(void){
    return pool.size();
}
