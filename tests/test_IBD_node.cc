#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "IBDmix/IBD_Stack.h"

TEST(IBDnode, CanGetNode){
    IBD_Node *top = NULL;
    top = get_node(12);
    top->next = NULL;
    ASSERT_EQ(top->position, 12);
    reclaim_all(top);
}

TEST(IBDnode, CanInsert){
    IBD_Node *top = NULL;
    push(top, get_node(3));
    push(top, get_node(2));
    push(top, get_node(1));

    ASSERT_EQ(stack_length(top), 3);
    int count = 1;
    for(struct IBD_Node* ptr=top; ptr != NULL; ptr=ptr->next)
        ASSERT_EQ(ptr->position, count++);
    reclaim_all(top);
}

TEST(IBDnode, CanReclaimAlloc){
    // based on execution, will have some free nodes left
    int pool_len = buff_size;
    IBD_Node *top = NULL;
    for(int i = 0; i < pool_len; i++)
        push(top, get_node(i));

    ASSERT_EQ(pool_len, stack_length(top));
    ASSERT_EQ(pool_length(), 0);
    for(int i = 0; i < 15; i++)
        push(top, get_node(i));
    ASSERT_EQ(pool_length(), pool_len - 15);

    ASSERT_EQ(15+pool_len, stack_length(top));
    IBD_Node *ptr = top;
    for(int i = 0; i < 14; i++)
        ptr = ptr->next;
    reclaim_after(ptr);
    ASSERT_EQ(stack_length(top), 15);
    ASSERT_EQ(pool_len*2-15, pool_length());

    //test no ops
    reclaim_after(ptr);
    ASSERT_EQ(stack_length(top), 15);
    ASSERT_EQ(pool_len*2-15, pool_length());

    reclaim_after(NULL);
    ASSERT_EQ(stack_length(top), 15);
    ASSERT_EQ(pool_len*2-15, pool_length());
    
    //reclaim all
    reclaim_all(top);
    ASSERT_EQ(stack_length(top), 0);
    ASSERT_EQ(pool_len*2, pool_length());
}

TEST(IBDnode, CanReclaimBetween){
    IBD_Node *top = NULL;
    for(int i = 0; i < 15; i++)
        push(top, get_node(i));
    IBD_Node *other = top;
    //remove a few
    for(int i = 0; i < 4; i++)
        other = other->next;
    reclaim_between(top, other);
    ASSERT_EQ(stack_length(top), 12);
    //remove to end
    for(other = top; other->next != NULL; other = other->next);
    reclaim_between(top, other);
    ASSERT_EQ(stack_length(top), 2);
    reclaim_all(top);
}

TEST(IBDnode, CanReverse){
    IBD_Node *top = NULL;
    IBD_Node *result = reverse(top);
    ASSERT_EQ(stack_length(result), 0);

    reclaim_all(result);

    top = NULL;
    push(top, get_node(1));
    result = reverse(top);
    ASSERT_EQ(stack_length(result), 1);

    reclaim_all(result);
    top = NULL;
    for(int i = 0; i < 15; i++)
        push(top, get_node(i));
    result = reverse(top);
    ASSERT_EQ(stack_length(result), 15);
    reclaim_all(result);
}

TEST(IBDnode, CanPushPop){
    IBD_Node *top = NULL;
    IBD_Node *result = NULL;
    int pool_len = pool_length();
    int pool = buff_size > pool_len ? buff_size : pool_len;

    push(top, get_node(1));
    push(top, get_node(2));
    push(top, get_node(3));
    ASSERT_EQ(stack_length(top), 3);
    ASSERT_EQ(pool-3, pool_length());

    result = pop(top);
    ASSERT_EQ(result->position, 3);
    reclaim_node(result);
    result = pop(top);
    ASSERT_EQ(result->position, 2);
    reclaim_node(result);
    result = pop(top);
    ASSERT_EQ(result->position, 1);
    reclaim_node(result);

    ASSERT_EQ(stack_length(top), 0);
    ASSERT_EQ(pool, pool_length());
}
