#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "minHeap.c"
#include <assert.h>


#define CUDACHECK(cmd) do { \
    cudaError_t e = cmd;    \
    if( e != cudaSuccess ) { \
    printf("Failed: Cuda error %s:%d '%s'\n", \
        __FILE__,__LINE__,cudaGetErrorString(e)); \
    exit(EXIT_FAILURE);     \
  } \
} while(0)

void swap(float *lhs, float *rhs)
{
    if (lhs == rhs)
        return;

    float tmp = *lhs;
    *lhs = *rhs;
    *rhs = tmp;
}

__global__ void fastJacc(unsigned int *allMols, unsigned int *queries, float *sims, int row, int size)
{
    // take steps over each fingerprint in allmols
    // aka grid-stride loop
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    // the cardinality is stored as the first int, always.
    int cardX = queries[65 * row];
    for (int i = index * 65; i < size * 65; i += 65 * stride) {
        int totalSize = 0;
        // unroll the loop into 16 loops of 4. 15% speedup.
#pragma unroll 4
        for (int j = 1; j < 65; j++)
        {
            unsigned int x = queries[(65 * row) + j];
            unsigned int y = allMols[i + j];
            totalSize += __popc(x & y);
        }
        int cardY = allMols[i];
        // alternate eq for jaccard score: (int(x, y) / |x| + |y| + int(x, y)).
        int jaccDenom = cardX + cardY - totalSize;
        // theoretically can have div by 0 error, however, jaccard score isn't defined for int(x, y) = 0,
        // and the chances of it happening are basically 0.
        // so don't waste the compute power to guard against it. if you're getting errors here,
        // your data is probably corrupted.
        sims[(size * row) + (i / 65)] = float(totalSize) / jaccDenom;
    }
}

// function to get the topk from an array of similarities
void populateMinHeap(float *array, minHeap *m, int size) {
    for (int i = 0; i < size; i++) {
        if (m->size < m->maxSize) {
            insertNode(m, array[i], i);
        }

        else if (array[i] > getMinNode(m)) {
            // heapify instead?
            deleteNode(m);
            insertNode(m, array[i], i);
        }
    }
}

int main(int argc, char* argv[])
{
    if (argc != 6) {
        printf("%s\n", "Must input 4 arguements, in order, input database, input queries, output file, int block size, int top_k.");
        printf("%s\n", "Example run: fastSearch_CUDA mols.bin inq.bin search.txt 1024 30");
        exit(EXIT_FAILURE);
    }

    FILE *database_ptr;
    database_ptr = fopen(argv[1], "rb"); // r for read, b for binary
    if (database_ptr == NULL) {
        printf("Couldn't find input database. Is the path correct?");
        exit(EXIT_FAILURE);
    }
    FILE *queries_ptr;
    queries_ptr = fopen(argv[2], "rb"); // r for read, b for binary
    if (queries_ptr == NULL) {
        printf("Couldn't find input queries. Is the path correct?");
        exit(EXIT_FAILURE);
    }

    // start reading the database file
    int size;
    // read num mols
    fread(&size, sizeof(int), 1, database_ptr);
    unsigned int *mols;
    // Allocate Unified Memory -- accessible from CPU or GPU
    cudaMallocManaged(&mols, sizeof(int) * size * (64 + 1));
    // actually read the file to unified memory
    int freadReturnVal = fread(mols, sizeof(int) * (64 + 1), size, database_ptr);
    if (freadReturnVal != size) {
        printf("%s", "Bad read on input database. Do you have enough memory? Is the file corrupted? Has it been preprocessed with convert_to_binary?");
        exit(EXIT_FAILURE);
    }
    fclose(database_ptr);

    // exact same as above
    int num_queries;
    fread(&num_queries, sizeof(int), 1, queries_ptr); // read num queries
    unsigned int *queries;
    cudaMallocManaged(&queries, sizeof(int) * num_queries * (64 + 1));
    freadReturnVal = fread(queries, sizeof(int) * (64 + 1), num_queries, queries_ptr);
    if (freadReturnVal != num_queries) {
        printf("%s", "Bad read on input queries. Do you have enough memory? Is the file corrupted? Has it been preprocessed with convert_to_binary?");
        exit(EXIT_FAILURE);
    }
    fclose(queries_ptr);

    int numToRun = num_queries;

    printf("%s", "Number of queries: ");
    printf("%d\n", num_queries);
    int k = atoi(argv[5]);

    clock_t start, end;
    double cpu_time_used;
    start = clock();
    float *sims;
    // width * height -- CUDA doesn't like 2d arrays
    cudaMallocManaged(&sims, (size * sizeof(float)) * (numToRun));
    // this value should probably be about 1024
    int blockSize = atoi(argv[4]);
    int numBlocks = (size + blockSize - 1) / blockSize;
    for (int i = 0; i < numToRun; i++)
    {
        fastJacc << <numBlocks, blockSize >> > (mols, queries, sims, i, size);
    }
    CUDACHECK(cudaDeviceSynchronize());

    struct minHeap **heaps = (struct minHeap**) malloc(sizeof(struct minHeap *) * numToRun);
    for (int i = 0; i < numToRun; i++) {
        heaps[i] = initMinHeap(k);
        populateMinHeap(&sims[(i * size)], heaps[i], size);
        // uncomment these if you want to look at your data
        // preorderTraversal(heaps[i], 0);
        // printf("%s\n", "");
    }

    end = clock();
    cpu_time_used = (((double)(end - start)) / CLOCKS_PER_SEC);
    printf("%s", "Time to run queries: ");
    printf("%f\n", cpu_time_used);

    FILE *out = fopen(argv[3], "w");
    for (int i = 0; i < numToRun; i++) {
        minHeap *currHeap = heaps[i];
        for (int j = 0; j < k; j++) {
            fprintf(out, "%s", "(");
            fprintf(out, "%.4f", currHeap->elem[j].data);
            fprintf(out, "%s", ", ");
            fprintf(out, "%d", currHeap->elem[j].idx);
            fprintf(out, "%s", ")");
            if (j + 1 < size)
                fprintf(out, "%s", ", ");
        }
        fprintf(out, "%s", "\n");
    }

    fclose(out);

    // Free memory
    free(heaps);
    cudaFree(queries);
    cudaFree(mols);
    cudaFree(sims);
    return 0;
}