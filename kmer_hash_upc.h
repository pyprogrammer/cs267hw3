#ifndef KMER_HASH_UPC_H
#define KMER_HASH_UPC_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include <upc.h>
#include "contig_generation.h"

typedef shared hash_table_t* hash_dir_t;
typedef shared memory_heap_t* mem_dir_t;

/* Creates hashtable for UPC */
/* Whenever something is called on a hash_table_t or memory_heap_t, instead call it on
 * hash_dir_t[MYTHREAD] and mem_dir_t[MYTHREAD]
 */
hash_dir_t *upc_create_hash_table(int64_t nEntries, mem_dir_t *memory_heap)
{
  hash_dir_t result;

  int64_t n_buckets = nEntries * LOAD_FACTOR;

  result = upc_all_alloc(THREADS,sizeof(hash_table_t));
  result[MYTHREAD].size = (n_buckets + THREADS - 1)/THREADS;
  result[MYTHREAD].table = (bucket_t*) upc_alloc( result[MYTHREAD].size * sizeof(bucket_t) );
  memset((void*)result[MYTHREAD].table,0,result[MYTHREAD].size * sizeof(bucket_t));

  if (result[MYTHREAD].table == NULL) {
     fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(bucket_t));
     exit(1);
  }

  /*
  memory_heap->heap = (kmer_t *) upc_all_alloc( (nEntries + THREADS - 1)/THREADS * sizeof(kmer_t) );
  if (memory_heap[MYTHREAD]->heap == NULL) {
     fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
     exit(1);
  }

  memory_heap[MYTHREAD]->posInHeap=0;
  */

  return result;
}

#endif
