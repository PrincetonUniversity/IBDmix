#include <iostream>
#include <stdlib.h>

struct Node {
    int data;
    struct Node *next;
};

struct Node* head = NULL;
int buff_size = 10;
void display(struct Node*);

void allocate(int buffer_size){
    std::cout << "allocate\n";
    struct Node* temp = (struct Node*) malloc(sizeof(struct Node) * buffer_size);
    head = temp;
    for (int i = 0; i < buffer_size; i++){
        temp[i].data = i;
        temp[i].next = &temp[i+1];
    }
    temp[buffer_size-1].next = NULL;
    display(head);
}

struct Node* get(int data){
    std::cout << "get\n";
    if (head == NULL)
        allocate(buff_size);

    struct Node* result = head;
    head = head->next;
    result->data = data;
    return result;
}

void insert(int data, struct Node *&start){
    std::cout << "insert\n";
    struct Node* new_node = get(data);
    new_node->next = start;
    start = new_node;
}

void display(struct Node* start){
    std::cout << "display\n";
    struct Node* ptr = start;
    while (ptr != NULL) {
        std::cout << ptr-> data << " ";
        ptr = ptr->next;
    }
}

int main(){
    std::cout << "starting\n";
    struct Node* list = NULL;
    insert(1, list);
    insert(2, list);
    insert(3, list);

    display(list);
}
