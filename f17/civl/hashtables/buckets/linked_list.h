#include <stdbool.h>

typedef struct Node{
    int val;
    struct Node* next;
} Node;

Node* create_list(int val);
bool append_list(Node* head, int val);
bool contains_list(Node* head, int val);
bool discard_list(Node** head, int val);
void print_list(Node* list);
void free_list(Node* list);

Node* list;
