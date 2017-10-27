#include <stdbool.h>

typedef struct node{
    int val;
    struct node* next;
} node;

node* create(int val);
bool append(node* head, int val);
bool contains(node* head, int val);
bool n_remove(node* head, int val);
void print_list(node* list);
void free_list(node* list);
