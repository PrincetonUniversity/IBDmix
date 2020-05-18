#pragma once

#include <iostream>
#include <vector>

struct IBD_Node {
  double cumulative_lod, lod;
  uint64_t position;
  unsigned char bitmask;
  IBD_Node *next;
};

class IBD_Stack {
 public:
  IBD_Stack() = default;
  explicit IBD_Stack(IBD_Node *top) : top(top) {}
  IBD_Stack(IBD_Node *top, IBD_Node *bottom) : top(top), start(bottom) {}

  void push(IBD_Node *new_node);
  IBD_Node *pop();
  const IBD_Node *getTop() const { return top; }

  void write(std::ostream &strm) const;
  void reverse();
  IBD_Stack getUnprocessed();

  bool empty() const { return top == nullptr; }
  bool isSingleton() const { return top->next == nullptr; }
  int size() const;

  void setEnd() { end = top; }
  bool topIsNewMax() const {
    return top->cumulative_lod >= end->cumulative_lod;
  }
  bool reachedMax() const { return top->cumulative_lod < 0; }
  bool isEnd(const IBD_Node *node) { return node == end; }

  uint64_t startPosition() const { return start->position; }
  uint64_t endPosition() const { return end->position; }
  double endLod() const { return end->cumulative_lod; }

  void getAllFrom(IBD_Stack *other);
  void getSegmentFrom(IBD_Stack *other);

 private:
  IBD_Node *end = nullptr, *start = nullptr, *top = nullptr;
};

std::ostream &operator<<(std::ostream &strm, const IBD_Stack &stack);

class IBD_Pool {
 public:
  explicit IBD_Pool(int initial_buffer = 1024);
  ~IBD_Pool();

  IBD_Node *get_node(uint64_t position, double lod = 0,
                     unsigned char bitmask = 0);
  int size() const;
  void reclaim_node(IBD_Node *node);
  void reclaim_segment(IBD_Stack *stack);
  void reclaim_stack(IBD_Stack *stack);

 private:
  int buffer_size;
  IBD_Stack pool;
  std::vector<IBD_Node *> alloc_ptrs;

  void allocate(int nodes);
};
