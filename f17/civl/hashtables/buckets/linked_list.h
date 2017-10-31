#include <stdbool.h>

typedef struct node{
    int val;
    struct node* next;
} node;

node* create_list(int val);
bool append_list(node* head, int val);
bool contains_list(node* head, int val);
bool discard_list(node** head, int val);
void print_list(node* list);
void free_list(node* list);

node* list;
