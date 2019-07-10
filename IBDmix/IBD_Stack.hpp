#pragma once

struct IBD_Node {
    double cumulative_lod, lod;
    unsigned long int position;
    struct IBD_Node *next;
};
void allocate(int buffer_size);
struct IBD_Node* get_node(unsigned long int position, double lod=0);
void push(struct IBD_Node *& top, struct IBD_Node * new_node);
struct IBD_Node * pop(struct IBD_Node *& top);
void reclaim_node(struct IBD_Node *node);
void reclaim_after(struct IBD_Node *start);
void reclaim_between(struct IBD_Node *start, struct IBD_Node *end);
void reclaim_all(struct IBD_Node *& top);
struct IBD_Node* reverse(IBD_Node* top);
int stack_length(struct IBD_Node * top);
void display_stack(struct IBD_Node * top);
int pool_length(void);
void display_pool(void);
