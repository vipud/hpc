#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "linked_list.h"

// initialize
node* create(int val){
    node* new_node;
    new_node = (node *) malloc(sizeof(node));
    new_node->val = val;
    new_node->next = NULL;
    return new_node;
}

// add a new node
bool append(node* head, int val){
    if(head == NULL){
        return false;
    }
    node* cursor = head; // traverse
    while(cursor->next != NULL){
        cursor = cursor->next;
    }
    node* next_node = create(val); 
    cursor->next = next_node;
    return true; 
}

// check if value is in the list
bool contains(node* head, int val){
    if(head == NULL){
        //printf("no values in list!\n");
        return false;
    }
    node* cursor = head; // traverse
    while(cursor != NULL){
        if(cursor->val == val){
            printf("Found %d in list.\n", cursor->val);
            return true;
        }
        cursor = cursor->next;
    }
    return false;
}

// remove a node, based on value
bool n_remove(node* head, int val){
    if(head == NULL){
        printf("List does not exist!\n");
        return false;
    }
    else if(contains(head,val)){ // remove
        node* cursor = head;
        while(cursor->next != NULL){
            node* temp = cursor->next;
            if(cursor->val == val){
                if(temp == NULL){
                    return true;
                } else {
                    return true;
                }
            } // else go to next
            cursor = cursor->next;
        }
    } 
    else{
        printf("Could not find %d.!\n", val);
    }
    return false;
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

int main(int argc, char* argv[]){
    int SIZE = 0;
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
