#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

typedef struct List{
    int val;
    struct List* next;
} List;

List* create_list(int val){
    List* new_node;
    new_node = (List *) malloc(sizeof(List));
    new_node->val = val; // -1 = null value
    new_node->next = NULL;
    return new_node;
}

// add a new List
bool append_list(List* head, int val){
    if(head == NULL){
        head = create_list(-1);
        return true;
    }
    if(head->val == -1){
        head->val = val;
        return true;
    } 
    List* cursor = head; // traverse
    while(cursor->next != NULL){
        cursor = cursor->next;
    }
    List* next_node = create_list(val);
    cursor->next = next_node;
    return true;
}

// check if value is in the list
bool contains_list(List* head, int val){
    if(head == NULL){
        //printf("no values in list!\n");
        return false;
    }
    List* cursor = head; // traverse
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

// discard a List, based on value
bool remove_list(List** head, int val){
    List* cursor = *head;
    List* temp;

    // check front of linked list
    if(cursor == NULL){
        printf("list does not exist!\n");
        return false;
    }
    else if(cursor->val == val){
        if(cursor->next == NULL){ // if head, replace w/ -1
            cursor->val = -1;
            return true;
        } else {
            temp = *head;
            *head = (*head)->next;
            printf("deleting %d\n", temp->val);
            free(temp);
            return true;
        }
    }
    List* cursor_next = (*head)->next;;
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
void print_list(List* list){
    List* temp = list;
    while(temp != NULL){
        printf("%d ", temp->val);
        temp = temp->next;
    }
    printf("\n");
}

// free memory of whole list
void free_list(List* list){
    List* temp;
    List* cursor;
    if(list != NULL){
        cursor = list->next;
        list->next = NULL;
        while(cursor != NULL){
            temp = cursor->next;
            free(cursor);
            cursor = temp;
        }
    }
    free(list);
}
