#include <gtest/gtest.h>
#include <gmock/gmock.h>

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
    for (struct IBD_Node* ptr=stack.top; ptr != nullptr; ptr=ptr->next)
        ASSERT_EQ(ptr->position, count++);
}

TEST(IBDpool, CanReclaimAlloc) {
    int pool_len = 5;
    IBD_Stack stack;
    IBD_Pool pool(pool_len);
    for (int i = 0; i < pool_len; i++)
        stack.push(pool.get_node(i));

    ASSERT_EQ(pool_len, stack.size());
    ASSERT_EQ(pool.size(), 0);

    for (int i = 0; i < 6; i++)
        stack.push(pool.get_node(i));
    ASSERT_EQ(pool.size(), 9);  // should have alloc'd 5 then 10

    ASSERT_EQ(stack.size(), 11);

    IBD_Node *ptr = stack.top;
    for (int i = 0; i < 4; i++)
        ptr = ptr->next;

    pool.reclaim_after(ptr);
    ASSERT_EQ(stack.size(), 5);
    ASSERT_EQ(pool.size(), 15);

    // test no ops
    pool.reclaim_after(ptr);
    ASSERT_EQ(stack.size(), 5);
    ASSERT_EQ(pool.size(), 15);

    pool.reclaim_after(nullptr);
    ASSERT_EQ(stack.size(), 5);
    ASSERT_EQ(pool.size(), 15);

    // reclaim all
    pool.reclaim_all(&stack.top);
    ASSERT_EQ(stack.size(), 0);
    ASSERT_EQ(pool.size(), 20);
}

TEST(IBDpool, CanReclaimBetween) {
    IBD_Stack stack;
    IBD_Pool pool(5);

    for (int i = 0; i < 15; i++)
        stack.push(pool.get_node(i));

    // remove a few
    IBD_Node *other = stack.top;
    for (int i = 0; i < 4; i++)
        other = other->next;
    pool.reclaim_between(stack.top, other);
    ASSERT_EQ(stack.size(), 12);

    // remove to end
    for (other = stack.top; other->next != nullptr; other = other->next) {}
    pool.reclaim_between(stack.top, other);

    // contains top and end
    ASSERT_EQ(stack.top->position, 14);
    ASSERT_EQ(other->position, 0);
    ASSERT_EQ(stack.size(), 2);
}

TEST(IBDstack, CanReverse) {
    IBD_Stack stack;
    IBD_Pool pool(5);

    ASSERT_EQ(stack.size(), 0);
    stack.reverse();
    ASSERT_EQ(stack.size(), 0);

    pool.reclaim_all(&stack.top);

    stack.push(pool.get_node(1));
    ASSERT_EQ(stack.size(), 1);
    stack.reverse();
    ASSERT_EQ(stack.size(), 1);
    ASSERT_EQ(stack.top->position, 1);

    pool.reclaim_all(&stack.top);

    for (int i = 0; i < 15; i++)
        stack.push(pool.get_node(i));
    ASSERT_EQ(stack.top->position, 14);
    stack.reverse();
    ASSERT_EQ(stack.top->position, 0);
    ASSERT_EQ(stack.size(), 15);
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
    ASSERT_STREQ(print_out.str().c_str(), "\n");
    print_out.str("");
    print_out.clear();

    stack.push(pool.get_node(1));
    print_out << stack;
    ASSERT_STREQ(print_out.str().c_str(), "1 \n");
    print_out.str("");
    print_out.clear();

    stack.push(pool.get_node(2));
    print_out << stack;
    ASSERT_STREQ(print_out.str().c_str(), "2 1 \n");
    print_out.str("");
    print_out.clear();

    stack.push(pool.get_node(3));
    print_out << stack;
    ASSERT_STREQ(print_out.str().c_str(), "3 2 1 \n");
    print_out.str("");
    print_out.clear();

    ASSERT_EQ(stack.size(), 3);

    IBD_Node *result = stack.pop();
    pool.reclaim_node(result);
    print_out << stack;
    ASSERT_STREQ(print_out.str().c_str(), "2 1 \n");
}
