/*
	File:   minHeap.c
	Desc:   Program showing various operations on a binary min heap
	Modified by me, original implementation by Robin Thomas at 
	https://github.com/robin-thomas/min-heap/blob/master/minHeap.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>


#define LCHILD(x) 2 * x + 1
#define RCHILD(x) 2 * x + 2
#define PARENT(x) (x - 1) / 2

typedef struct node {
	float data;
	int idx;
} node;

typedef struct minHeap {
	int size;
	int maxSize;
	node *elem;
} minHeap;


/*
	Function to initialize the min heap with size = 0, maxsize = maxSize
*/
struct minHeap * initMinHeap(int maxSize) {
	struct minHeap *hp;
	hp = (struct minHeap*) malloc(sizeof(struct minHeap));
	hp->size = 0;
	hp->maxSize = maxSize;
	hp->elem = (node*)malloc(sizeof(node) * maxSize);
	return hp;
}


/*
	Function to swap data within two nodes of the min heap using pointers
*/
void swap(node *n1, node *n2) {
	node temp = *n1;
	*n1 = *n2;
	*n2 = temp;
}


/*
	Heapify function is used to make sure that the heap property is never violated
	In case of deletion of a node, or creating a min heap from an array, heap property
	may be violated. In such cases, heapify function can be called to make sure that
	heap property is never violated
*/
void heapify(minHeap *hp, int i) {
	int smallest = (LCHILD(i) < hp->size && hp->elem[LCHILD(i)].data < hp->elem[i].data) ? LCHILD(i) : i;
	if (RCHILD(i) < hp->size && hp->elem[RCHILD(i)].data < hp->elem[smallest].data) {
		smallest = RCHILD(i);
	}
	if (smallest != i) {
		swap(&(hp->elem[i]), &(hp->elem[smallest]));
		heapify(hp, smallest);
	}
}

/*
	Function to insert a node into the min heap
*/
void insertNode(minHeap *hp, float data, int idx) {
	assert(hp->size < hp->maxSize);

	node nd;
	nd.data = data;
	nd.idx = idx;

	int i = (hp->size)++;
	while (i && nd.data < hp->elem[PARENT(i)].data) {
		hp->elem[i] = hp->elem[PARENT(i)];
		i = PARENT(i);
	}
	hp->elem[i] = nd;
}


/*
	Function to delete a node from the min heap
	It shall remove the root node, and place the last node in its place
	and then call heapify function to make sure that the heap property
	is never violated
*/
void deleteNode(minHeap *hp) {
	if (hp->size) {
		// printf("Deleting node %d\n\n", hp->elem[0].data) ;
		hp->elem[0] = hp->elem[--(hp->size)];
		heapify(hp, 0);
	}
	else {
		printf("\nMin Heap is empty!\n");
		free(hp->elem);
	}
}


/*
	Function to get maximum node from a min heap
	The maximum node shall always be one of the leaf nodes. So we shall recursively
	move through both left and right child, until we find their maximum nodes, and
	compare which is larger. It shall be done recursively until we get the maximum
	node
*/
int getMaxNode(minHeap *hp, int i) {
	if (LCHILD(i) >= hp->size) {
		return hp->elem[i].data;
	}

	int l = getMaxNode(hp, LCHILD(i));
	int r = getMaxNode(hp, RCHILD(i));

	if (l >= r) {
		return l;
	}
	else {
		return r;
	}
}

float getMinNode(minHeap *hp) {
	return hp->elem[0].data;
}
/*
	Function to clear the memory allocated for the min heap
*/
void deleteMinHeap(minHeap *hp) {
	free(hp->elem);
}

/*
	Function to display all the nodes in the min heap by doing a preorder traversal
*/
void preorderTraversal(minHeap *hp, int i) {
	printf("%f ", hp->elem[i].data);
	printf("%d ", hp->elem[i].idx);
	if (LCHILD(i) < hp->size) {
		preorderTraversal(hp, LCHILD(i));
	}
	if (RCHILD(i) < hp->size) {
		preorderTraversal(hp, RCHILD(i));
	}
}

