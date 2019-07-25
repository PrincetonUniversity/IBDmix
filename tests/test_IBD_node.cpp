#define BOOST_TEST_MODULE IBD_Node
#include <boost/test/included/unit_test.hpp>
#include "../IBDmix/IBD_Stack.hpp"

BOOST_AUTO_TEST_CASE(test_get_node){
    IBD_Node *top = NULL;
    top = get_node(12);
    top->next = NULL;
    BOOST_REQUIRE_EQUAL(12, top->position);
    reclaim_all(top);
}

BOOST_AUTO_TEST_CASE(test_insert){
    IBD_Node *top = NULL;
    push(top, get_node(3));
    push(top, get_node(2));
    push(top, get_node(1));

    BOOST_REQUIRE_EQUAL(3, stack_length(top));
    int count = 1;
    for(struct IBD_Node* ptr=top; ptr != NULL; ptr=ptr->next)
        BOOST_REQUIRE_EQUAL(ptr->position, count++);
    reclaim_all(top);
}

BOOST_AUTO_TEST_CASE(test_alloc_reclaim){
    // based on execution, will have some free nodes left
    int pool_len = pool_length();
    IBD_Node *top = NULL;
    for(int i = 0; i < pool_len; i++)
        push(top, get_node(i));

    BOOST_REQUIRE_EQUAL(pool_len, stack_length(top));
    BOOST_REQUIRE_EQUAL(0, pool_length());
    for(int i = 0; i < 15; i++)
        push(top, get_node(i));
    BOOST_REQUIRE_EQUAL(5, pool_length());

    BOOST_REQUIRE_EQUAL(15+pool_len, stack_length(top));
    IBD_Node *ptr = top;
    for(int i = 0; i < 14; i++)
        ptr = ptr->next;
    reclaim_after(ptr);
    BOOST_REQUIRE_EQUAL(15, stack_length(top));
    BOOST_REQUIRE_EQUAL(5+pool_len, pool_length());

    //test no ops
    reclaim_after(ptr);
    BOOST_REQUIRE_EQUAL(15, stack_length(top));
    BOOST_REQUIRE_EQUAL(5+pool_len, pool_length());

    reclaim_after(NULL);
    BOOST_REQUIRE_EQUAL(15, stack_length(top));
    BOOST_REQUIRE_EQUAL(5+pool_len, pool_length());
    
    //reclaim all
    reclaim_all(top);
    BOOST_REQUIRE_EQUAL(0, stack_length(top));
    BOOST_REQUIRE_EQUAL(20+pool_len, pool_length());
}

BOOST_AUTO_TEST_CASE(test_reclaim_between){
    IBD_Node *top = NULL;
    for(int i = 0; i < 15; i++)
        push(top, get_node(i));
    IBD_Node *other = top;
    //remove a few
    for(int i = 0; i < 4; i++)
        other = other->next;
    reclaim_between(top, other);
    BOOST_REQUIRE_EQUAL(12, stack_length(top));
    //remove to end
    for(other = top; other->next != NULL; other = other->next);
    reclaim_between(top, other);
    BOOST_REQUIRE_EQUAL(2, stack_length(top));
    reclaim_all(top);
}

BOOST_AUTO_TEST_CASE(test_reverse){
    IBD_Node *top = NULL;
    IBD_Node *result = reverse(top);
    BOOST_REQUIRE_EQUAL(0, stack_length(result));

    reclaim_all(result);

    top = NULL;
    push(top, get_node(1));
    result = reverse(top);
    BOOST_REQUIRE_EQUAL(1, stack_length(result));

    reclaim_all(result);
    top = NULL;
    for(int i = 0; i < 15; i++)
        push(top, get_node(i));
    result = reverse(top);
    BOOST_REQUIRE_EQUAL(15, stack_length(result));
    reclaim_all(result);
}

BOOST_AUTO_TEST_CASE(test_push_pop){
    IBD_Node *top = NULL;
    IBD_Node *result = NULL;
    int pool = pool_length();

    push(top, get_node(1));
    push(top, get_node(2));
    push(top, get_node(3));
    BOOST_REQUIRE_EQUAL(3, stack_length(top));
    BOOST_REQUIRE_EQUAL(pool-3, pool_length());

    result = pop(top);
    BOOST_REQUIRE_EQUAL(3, result->position);
    reclaim_node(result);
    result = pop(top);
    BOOST_REQUIRE_EQUAL(2, result->position);
    reclaim_node(result);
    result = pop(top);
    BOOST_REQUIRE_EQUAL(1, result->position);
    reclaim_node(result);

    BOOST_REQUIRE_EQUAL(0, stack_length(top));
    BOOST_REQUIRE_EQUAL(pool, pool_length());
}
