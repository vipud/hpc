#include <stdio.h>
#include <stdlib.h>

int SIZE; // store size

// each node stores a value
typedef struct node{     
    int val;
    struct node* next;
} node;

// initialize
node* create(int val){
    node* new_node;
    new_node = (node *) malloc(sizeof(node));
    new_node->val = val;
    new_node->next = NULL;
    return new_node;
}

// add a new node
int append(node* head, int val){
    if(head == NULL){
        return 0;
    }
    node* cursor = head; // traverse
    while(cursor->next != NULL){
        cursor = cursor->next;
    }
    node* next_node = create(val); 
    cursor->next = next_node;
    return 1; 
}

// check if value is in the list
int contains(struct node* head, int val){
    if(head == NULL){
        //printf("no values in list!\n");
        return 0;
    }
    node* cursor = head; // traverse
    while(cursor != NULL){
        if(cursor->val == val){
            //printf("Found %d in list.\n", cursor->val);
            return 1;
        }
        cursor = cursor->next;
    }
    //printf("Could not find %d in list.\n", val);
    return 0;
}

// remove a node, based on value
int n_remove(struct node* head, int val){
    if(!contains(head,val)){
        printf("%d not found.\n", val);
        return 0;
    } 
    // must be in list
   node* cursor = head; 
   while(cursor != NULL){
       printf("cursor = %d\n", cursor->val);
       if(cursor->val == val){
           // remove
           node* temp = cursor->next;       
           // shift 1 to right
           if(temp->next != NULL){
               cursor->val = temp->next->val;
               cursor->next = temp->next;
               free(temp);
           }
           else{
               
           return 1;
       }
       cursor = cursor->next;
   }
   printf("Something's wrong!\n");
   return 0;
}

// traverse and print vals
void print_list(node* list){
    node* temp = list;
    while(temp != NULL){
        printf("%d ", temp->val);
        temp = temp->next;
    }
    printf("\n");
}

// free memory of whole list
void free_list(node* list){
    node* temp;
    node* cursor;
    if(list != NULL){
        cursor = list->next;
        list->next = NULL;
        while(cursor != NULL){
            temp = cursor->next;
            free(cursor);
            cursor = temp;
        }
    }
}

int main(){
    SIZE = 0;
    struct node* head;
    head = create(1);
    print_list(head);
    append(head, 2);
    print_list(head);
    append(head, 4);
    print_list(head);
    append(head, 16);
    print_list(head);
    contains(head, 4);
    contains(head, 5);
    contains(head, 1);
    contains(head, 16);
    n_remove(head, 15);
    n_remove(head, 16);
    free_list(head);
    return 0;
}
