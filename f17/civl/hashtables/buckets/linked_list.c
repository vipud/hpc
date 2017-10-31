#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "linked_list.h"

// initialize
node* create_list(int val){
    node* new_node;
    new_node = (node *) malloc(sizeof(node));
    new_node->val = val;
    new_node->next = NULL;
    return new_node;
}

// add a new node
bool append_list(node* head, int val){
    if(head == NULL){
        return false;
    }
    node* cursor = head; // traverse
    while(cursor->next != NULL){
        cursor = cursor->next;
    }
    node* next_node = create_list(val); 
    cursor->next = next_node;
    return true; 
}

// check if value is in the list
bool contains_list(node* head, int val){
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
    printf("%d not found in list.\n", val);
    return false;
}

// discard_list a node, based on value
bool discard_list(node** head, int val){
    node* cursor = *head;
    node* temp;

    // check front of linked list
    if(cursor == NULL){
        printf("list does not exist!\n");
        return false;
    } 
    else if(cursor->val == val){
        temp = *head;
        *head = (*head)->next;
        printf("deleting %d\n", temp->val);
        free(temp);
        return true;
    }
    node* cursor_next = (*head)->next;;
    while(cursor != NULL && cursor_next != NULL){ // run through list
        if(cursor_next->val == val){
            printf("deleting %d\n", cursor_next->val);
            temp = cursor_next;
            cursor->next = cursor_next->next;
            free(temp);
            return true;
        }
        cursor = cursor_next;
        cursor_next = cursor_next->next;
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
    head = create_list(1);
    print_list(head);
    append_list(head, 2);
    print_list(head);
    append_list(head, 4);
    print_list(head);
    append_list(head, 16);
    print_list(head);
    contains_list(head, 4);
    contains_list(head, 5);
    contains_list(head, 1);
    contains_list(head, 16);
    discard_list(&head, 1);
    print_list(head);
    discard_list(&head, 16);
    print_list(head);
    discard_list(&head, 15);
    discard_list(&head, 4);
    print_list(head);
    discard_list(&head, 2);
    print_list(head);
    discard_list(&head, 1);
    print_list(head);
    discard_list(&head, 1);
    free_list(head);
    return 0;
}

