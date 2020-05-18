#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "IBDmix/IBD_Stack.h"

TEST(IBDpool, CanGetNode) {
  IBD_Node *top = nullptr;
  IBD_Pool pool(5);
  top = pool.get_node(12);
  ASSERT_EQ(top->position, 12);
  ASSERT_EQ(top->next, nullptr);
  ASSERT_EQ(top->lod, 0);
  ASSERT_EQ(top->cumulative_lod, 0);
}

TEST(IBDstack, CanPush) {
  IBD_Stack stack;
  IBD_Pool pool(5);

  stack.push(pool.get_node(3));
  stack.push(pool.get_node(2));
  stack.push(pool.get_node(1));

  ASSERT_EQ(stack.size(), 3);
  int count = 1;
  for (const IBD_Node *ptr = stack.getTop(); ptr != nullptr; ptr = ptr->next)
    ASSERT_EQ(ptr->position, count++);
}

TEST(IBDpool, CanReclaimAlloc) {
  int pool_len = 5;
  IBD_Stack stack;
  IBD_Pool pool(pool_len);
  for (int i = 0; i < pool_len; i++) stack.push(pool.get_node(i));

  ASSERT_EQ(pool_len, stack.size());
  ASSERT_EQ(pool.size(), 0);

  for (int i = 0; i < 6; i++) stack.push(pool.get_node(i));
  ASSERT_EQ(pool.size(), 9);  // should have alloc'd 5 then 10

  ASSERT_EQ(stack.size(), 11);

  // reclaim all
  pool.reclaim_stack(&stack);
  ASSERT_EQ(stack.size(), 0);
  ASSERT_EQ(pool.size(), 20);

  pool.reclaim_stack(&stack);
  ASSERT_EQ(stack.size(), 0);
  ASSERT_EQ(pool.size(), 20);

  for (int i = 0; i < 6; i++) stack.push(pool.get_node(i));
  ASSERT_EQ(stack.size(), 6);
  ASSERT_EQ(pool.size(), 14);

  pool.reclaim_stack(&stack);
  ASSERT_EQ(stack.size(), 0);
  ASSERT_EQ(pool.size(), 20);
}

TEST(IBDpool, CanReclaimBetween) {
  IBD_Stack stack;
  IBD_Pool pool(5);

  for (int i = 0; i < 5; i++) stack.push(pool.get_node(i));
  stack.setEnd();
  for (int i = 5; i < 15; i++) stack.push(pool.get_node(i));

  // remove a few
  pool.reclaim_segment(&stack);
  ASSERT_EQ(stack.size(), 12);

  // remove to end
  stack.setEnd();
  pool.reclaim_segment(&stack);

  // contains top and end
  ASSERT_EQ(stack.getTop()->position, 14);
  ASSERT_EQ(stack.size(), 2);
}

TEST(IBDstack, CanReverse) {
  IBD_Stack stack;
  IBD_Pool pool(5);
  std::ostringstream print_out;

  ASSERT_EQ(stack.size(), 0);

  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(), "");
  print_out.str("");
  print_out.clear();

  stack.reverse();
  ASSERT_EQ(stack.size(), 0);

  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(), "");
  print_out.str("");
  print_out.clear();

  pool.reclaim_stack(&stack);

  stack.push(pool.get_node(1));
  ASSERT_EQ(stack.size(), 1);
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(), "1\t0\t0 <- top <- start <- end\n");
  print_out.str("");
  print_out.clear();

  stack.reverse();
  ASSERT_EQ(stack.size(), 1);
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(), "1\t0\t0 <- top <- start <- end\n");
  print_out.str("");
  print_out.clear();

  pool.reclaim_stack(&stack);

  for (int i = 0; i < 10; i++) stack.push(pool.get_node(i));
  stack.setEnd();
  for (int i = 10; i < 15; i++) stack.push(pool.get_node(i));
  ASSERT_EQ(stack.getTop()->position, 14);
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "14\t0\t0 <- top\n"
               "13\t0\t0\n"
               "12\t0\t0\n"
               "11\t0\t0\n"
               "10\t0\t0\n"
               "9\t0\t0 <- end\n"
               "8\t0\t0\n"
               "7\t0\t0\n"
               "6\t0\t0\n"
               "5\t0\t0\n"
               "4\t0\t0\n"
               "3\t0\t0\n"
               "2\t0\t0\n"
               "1\t0\t0\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  stack.reverse();
  ASSERT_EQ(stack.getTop()->position, 0);
  ASSERT_EQ(stack.size(), 15);
  print_out << stack;
  // top and start are at the top and bottom, end stays in the same spot
  ASSERT_STREQ(print_out.str().c_str(),
               "0\t0\t0 <- top\n"
               "1\t0\t0\n"
               "2\t0\t0\n"
               "3\t0\t0\n"
               "4\t0\t0\n"
               "5\t0\t0\n"
               "6\t0\t0\n"
               "7\t0\t0\n"
               "8\t0\t0\n"
               "9\t0\t0 <- end\n"
               "10\t0\t0\n"
               "11\t0\t0\n"
               "12\t0\t0\n"
               "13\t0\t0\n"
               "14\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  IBD_Stack stack2 = stack.getUnprocessed();
  print_out << stack;

  // original stack is reversed
  ASSERT_STREQ(print_out.str().c_str(),
               "14\t0\t0 <- top\n"
               "13\t0\t0\n"
               "12\t0\t0\n"
               "11\t0\t0\n"
               "10\t0\t0\n"
               "9\t0\t0 <- start <- end\n");
  print_out.str("");
  print_out.clear();

  print_out << stack2;
  // new stack has top and bottom but no end
  ASSERT_STREQ(print_out.str().c_str(),
               "8\t0\t0 <- top\n"
               "7\t0\t0\n"
               "6\t0\t0\n"
               "5\t0\t0\n"
               "4\t0\t0\n"
               "3\t0\t0\n"
               "2\t0\t0\n"
               "1\t0\t0\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();
}

TEST(IBDstack, CanPushPop) {
  IBD_Stack stack;
  IBD_Pool pool(5);

  stack.push(pool.get_node(1));
  stack.push(pool.get_node(2));
  stack.push(pool.get_node(3));
  ASSERT_EQ(stack.size(), 3);

  IBD_Node *result = stack.pop();
  ASSERT_EQ(result->position, 3);
  pool.reclaim_node(result);
  result = stack.pop();
  ASSERT_EQ(result->position, 2);
  pool.reclaim_node(result);
  result = stack.pop();
  ASSERT_EQ(result->position, 1);
  pool.reclaim_node(result);

  ASSERT_EQ(stack.size(), 0);
}

TEST(IBDstack, CanPrint) {
  IBD_Stack stack;
  IBD_Pool pool(5);
  std::ostringstream print_out;

  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(), "");
  print_out.str("");
  print_out.clear();

  stack.push(pool.get_node(1));
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(), "1\t0\t0 <- top <- start <- end\n");
  print_out.str("");
  print_out.clear();

  stack.push(pool.get_node(2));
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "2\t0\t0 <- top\n"
               "1\t0\t0 <- start <- end\n");
  print_out.str("");
  print_out.clear();

  stack.push(pool.get_node(3));
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "3\t0\t0 <- top\n"
               "2\t0\t0\n"
               "1\t0\t0 <- start <- end\n");
  print_out.str("");
  print_out.clear();

  ASSERT_EQ(stack.size(), 3);

  IBD_Node *result = stack.pop();
  pool.reclaim_node(result);
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "2\t0\t0 <- top\n"
               "1\t0\t0 <- start <- end\n");
}

TEST(IBDstack, CanGetAll) {
  IBD_Stack stack, stack2;
  IBD_Pool pool(15);
  std::ostringstream print_out;

  for (int i = 0; i < 5; i++) stack.push(pool.get_node(i));

  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "4\t0\t0 <- top\n"
               "3\t0\t0\n"
               "2\t0\t0\n"
               "1\t0\t0\n"
               "0\t0\t0 <- start <- end\n");
  print_out.str("");
  print_out.clear();

  stack2.getAllFrom(&stack);

  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(), "");
  print_out.str("");
  print_out.clear();

  print_out << stack2;
  // note that start is set properly but end is left undefined
  ASSERT_STREQ(print_out.str().c_str(),
               "4\t0\t0 <- top\n"
               "3\t0\t0\n"
               "2\t0\t0\n"
               "1\t0\t0\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  // stack is empty so this should no op
  stack2.getAllFrom(&stack);

  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(), "");
  print_out.str("");
  print_out.clear();

  print_out << stack2;
  // note that start is set properly but end is left undefined
  ASSERT_STREQ(print_out.str().c_str(),
               "4\t0\t0 <- top\n"
               "3\t0\t0\n"
               "2\t0\t0\n"
               "1\t0\t0\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  stack2.setEnd();
  for (int i = 5; i < 10; i++) stack.push(pool.get_node(i));

  stack2.getAllFrom(&stack);

  print_out << stack2;
  ASSERT_STREQ(print_out.str().c_str(),
               "9\t0\t0 <- top\n"
               "8\t0\t0\n"
               "7\t0\t0\n"
               "6\t0\t0\n"
               "5\t0\t0\n"
               "4\t0\t0 <- end\n"
               "3\t0\t0\n"
               "2\t0\t0\n"
               "1\t0\t0\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  stack.push(pool.get_node(11));
  stack2.getAllFrom(&stack);

  print_out << stack2;
  ASSERT_STREQ(print_out.str().c_str(),
               "11\t0\t0 <- top\n"
               "9\t0\t0\n"
               "8\t0\t0\n"
               "7\t0\t0\n"
               "6\t0\t0\n"
               "5\t0\t0\n"
               "4\t0\t0 <- end\n"
               "3\t0\t0\n"
               "2\t0\t0\n"
               "1\t0\t0\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();
}

TEST(IBDstack, CanGetSegmentFrom) {
  IBD_Stack stack, stack2;
  IBD_Pool pool(15);
  std::ostringstream print_out;

  stack2.getSegmentFrom(&stack);
  ASSERT_EQ(stack.size(), 0);
  ASSERT_EQ(stack2.size(), 0);

  for (int i = 0; i < 5; i++) stack.push(pool.get_node(i));
  stack.setEnd();
  stack.push(pool.get_node(5));

  stack2.getSegmentFrom(&stack);
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "5\t0\t0 <- top\n"
               "4\t0\t0 <- end\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  print_out << stack2;
  ASSERT_STREQ(print_out.str().c_str(),
               "3\t0\t0 <- top\n"
               "2\t0\t0\n"
               "1\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  // repeat with end and start adjacent
  stack2.getSegmentFrom(&stack);
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "5\t0\t0 <- top\n"
               "4\t0\t0 <- end\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  print_out << stack2;
  ASSERT_STREQ(print_out.str().c_str(),
               "3\t0\t0 <- top\n"
               "2\t0\t0\n"
               "1\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  pool.reclaim_stack(&stack2);

  // repeat with end and start adjacent
  stack2.getSegmentFrom(&stack);
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "5\t0\t0 <- top\n"
               "4\t0\t0 <- end\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  print_out << stack2;
  ASSERT_STREQ(print_out.str().c_str(), "");
  print_out.str("");
  print_out.clear();

  // now one node between
  stack.setEnd();
  stack2.getSegmentFrom(&stack);
  print_out << stack;
  ASSERT_STREQ(print_out.str().c_str(),
               "5\t0\t0 <- top <- end\n"
               "0\t0\t0 <- start\n");
  print_out.str("");
  print_out.clear();

  print_out << stack2;
  ASSERT_STREQ(print_out.str().c_str(), "4\t0\t0 <- top <- start\n");
  print_out.str("");
  print_out.clear();
}
