#ifndef __MIN_HEAP_H__
#define __MIN_HEAP_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define LCHILD(x) 2 * x + 1
#define RCHILD(x) 2 * x + 2
#define PARENT(x) (x - 1) / 2

typedef struct node {
    double data;
    int index;
    int max_index;
} node ;

typedef struct minHeap {
    int size ;
    node *elem ;
} minHeap ;

minHeap initMinHeap(int size);
void swap(node *n1, node *n2);
void heapify(minHeap *hp, int i);
void buildMinHeap(minHeap *hp, double *arr, int*displacement, int*counts, int size);
void insertNode(minHeap *hp, double data, int index, int max_index);
void deleteNode(minHeap *hp);
int getMaxNode(minHeap *hp, int i);
void deleteMinHeap(minHeap *hp);
void inorderTraversal(minHeap *hp, int i);
void preorderTraversal(minHeap *hp, int i);
void postorderTraversal(minHeap *hp, int i) ;
void levelorderTraversal(minHeap *hp);
#endif