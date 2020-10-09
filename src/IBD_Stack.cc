#include "IBDmix/IBD_Stack.h"

void IBD_Stack::push(IBD_Node *new_node) {
  new_node->next = top;
  top = new_node;

  if (isSingleton()) {
    start = end = top;
    top->cumulative_lod = top->lod;
  } else {
    top->cumulative_lod = top->next->cumulative_lod + top->lod;
  }
}

IBD_Node *IBD_Stack::pop() {
  IBD_Node *result = top;
  top = top->next;
  if (result == start) start = nullptr;
  if (result == end) end = nullptr;
  return result;
}

int IBD_Stack::size() const {
  int count = 0;
  for (IBD_Node *ptr = top; ptr != nullptr; ptr = ptr->next) count++;
  return count;
}

void IBD_Stack::write(std::ostream &strm) const {
  for (struct IBD_Node *ptr = top; ptr != nullptr; ptr = ptr->next) {
    strm << ptr->position << "\t" << ptr->lod << "\t" << ptr->cumulative_lod;
    if (ptr == top) strm << " <- top";
    if (ptr == start) strm << " <- start";
    if (ptr == end) strm << " <- end";
    strm << "\n";
  }
}

std::ostream &operator<<(std::ostream &strm, const IBD_Stack &stack) {
  stack.write(strm);
  return strm;
}

void IBD_Stack::reverse() {
  // reverse the order of the list
  // performed in place, end kept at same node
  IBD_Node *result = nullptr, *temp, *old_end = end;
  while (top != nullptr) {
    // push top onto result
    temp = pop();
    if (result == nullptr) start = temp;  // first node, set start
    temp->next = result;
    result = temp;
  }

  top = result;
  end = old_end;
}

IBD_Stack IBD_Stack::getUnprocessed() {
  // get nodes between top and end as a new, reversed stack
  reverse();
  // get nodes after end
  IBD_Stack unprocessed(end->next, start);
  if (end != nullptr) end->next = nullptr;
  start = end;
  return unprocessed;
}

void IBD_Stack::getAllFrom(IBD_Stack *other) {
  if (other->empty()) return;

  // the start is always at the bottom of the stack
  other->start->next = top;
  top = other->top;

  if (start == nullptr) start = other->start;

  other->start = other->top = other->end = nullptr;
}

void IBD_Stack::getSegmentFrom(IBD_Stack *other) {
  // move between start and end (exclusive) to this
  if (other->end == nullptr || other->end == other->start) return;

  IBD_Node *ptr = other->end;
  while (ptr->next != other->start) ptr = ptr->next;

  if (ptr == other->end) return;

  ptr->next = top;
  top = other->end->next;
  other->end->next = other->start;

  if (start == nullptr) start = ptr;
}

IBD_Pool::IBD_Pool(int initial_buffer) : buffer_size(initial_buffer) {
  allocate(buffer_size);
}

void IBD_Pool::allocate(int nodes) {
  IBD_Node *allocation =
      reinterpret_cast<IBD_Node *>(malloc(sizeof(IBD_Node) * nodes));
  // setup allocation as a linked list
  for (int i = 0; i < nodes; i++) {
    allocation[i].next = &allocation[i + 1];
  }
  allocation[nodes - 1].next = nullptr;
  alloc_ptrs.push_back(allocation);

  // interpret as a stack and move to pool
  IBD_Stack stack(allocation, &allocation[nodes - 1]);
  pool.getAllFrom(&stack);
}

IBD_Pool::~IBD_Pool() {
  for (auto &ptr : alloc_ptrs) free(ptr);
  alloc_ptrs.clear();
}

IBD_Node *IBD_Pool::get_node(uint64_t position, double lod,
                             unsigned char bitmask) {
  if (pool.empty()) {
    allocate(buffer_size);
    buffer_size <<= 1;  // double next request
  }

  IBD_Node *result = pool.pop();

  result->position = position;
  result->lod = lod;
  result->bitmask = bitmask;
  result->cumulative_lod = lod;
  result->next = nullptr;
  return result;
}

void IBD_Pool::reclaim_node(IBD_Node *node) { pool.push(node); }

void IBD_Pool::reclaim_segment(IBD_Stack *stack) { pool.getSegmentFrom(stack); }

void IBD_Pool::reclaim_stack(IBD_Stack *stack) { pool.getAllFrom(stack); }

int IBD_Pool::size(void) const { return pool.size(); }
