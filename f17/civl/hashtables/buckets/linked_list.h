typedef struct node{
    int val;
    struct node* next;
} node;

node* create(int val);

int append(node* head, int val);
int contains(node* head, int val);
int n_remove(node* head, int val);
void print_list(node* list);
void free_list(node* list);

