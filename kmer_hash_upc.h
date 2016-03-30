#ifndef KMER_HASH_UPC_H
#define KMER_HASH_UPC_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include <upc.h>
#include "contig_upc.h"
#include "packingDNAseq.h"

typedef shared hash_table_t* hash_dir_t;
typedef shared memory_heap_t* mem_dir_t;

/* Creates hashtable for UPC */
/* Whenever something is called on a hash_table_t or memory_heap_t, instead call it on
 * hash_dir_t[MYTHREAD] and mem_dir_t[MYTHREAD]
 */
shared hash_table_t *upc_create_hash_table(int64_t nEntries, shared memory_heap_t **memory_heap)
{
  shared hash_table_t *result;
  shared bucket_t *global_tables; 
  shared bucket_t *heaps;

  int64_t n_buckets = nEntries * LOAD_FACTOR;
  int64_t my_size = (n_buckets + THREADS - 1)/THREADS;
  int64_t heap_size = ((nEntries + THREADS - 1)/THREADS);

  result = upc_all_alloc(THREADS,sizeof(hash_table_t));
  global_tables = (shared bucket_t*) upc_all_alloc( THREADS, my_size * sizeof(bucket_t) );
  (result+MYTHREAD)->size = my_size;
  (result+MYTHREAD)->table = (shared bucket_t*) &global_tables[MYTHREAD*my_size];

  *memory_heap = (shared memory_heap_t*) upc_all_alloc(THREADS,sizeof(memory_heap_t));
  heaps = (shared bucket_t *) upc_all_alloc( THREADS, heap_size * sizeof(kmer_t) );
  memory_heap[0][MYTHREAD].heap = (shared kmer_t*) (((shared char*) heaps) + (MYTHREAD * heap_size) * sizeof(kmer_t));
  memory_heap[0][MYTHREAD].posInHeap = 0;

  for(int i=0;i<THREADS;i++)
  {
    fprintf(stderr,"I am %d: addr 0x%lx posInHeap: %d\n",MYTHREAD,memory_heap[0][i].heap,memory_heap[0][i].posInHeap);
  }

  return result;
}

/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
  if(size < 8)
  {
    unsigned long hashval;
    hashval = 5381;
    for(int i = 0; i < size; i++) {
      hashval = seq[i] + (hashval << 5) + hashval;
    }
    return hashval % hashtable_size;
  }

  uint64_t lol = *((uint64_t*) seq);
  
  return lol % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
  return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
shared kmer_t* lookup_kmer_upc(shared hash_table_t *hashtable, shared memory_heap_t *memory_heap, const unsigned char *kmer)
{
  char packedKmer[KMER_PACKED_LENGTH];
  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
  int64_t hashval = hashkmer(hashtable->size * THREADS, (char*) packedKmer);
  int64_t which = hashval % THREADS;

  hashval = hashval / THREADS;

  hashtable += which;

  shared kmer_t *result;
  
  result = hashtable->table[hashval].head;
  unsigned char cmp[KMER_PACKED_LENGTH+1];
  cmp[KMER_PACKED_LENGTH] = (unsigned char) 0;
  
  while(result!=NULL) {
    upc_memget(&cmp,result->kmer,KMER_PACKED_LENGTH);
    fprintf(stderr,"cmp %s\n",cmp);
    if( memcmp(packedKmer, cmp, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
      return result;
    }
    result = result->next;
  }
  return NULL;
}

kmer_t* lookup_kmer(hash_table_t *hashtable, const unsigned char *kmer)
{
  return NULL;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
shared kmer_t* add_kmer(shared hash_table_t *tables, shared memory_heap_t *heaps,
    const unsigned char *kmer, char left_ext, char right_ext)
{
  /* Pack a k-mer sequence appropriately */
  shared hash_table_t *hashtable;
  shared memory_heap_t *memory_heap;
  char packedKmer[KMER_PACKED_LENGTH];

  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);

  int64_t hashval = hashkmer(tables->size * THREADS, (char*) packedKmer);

  int64_t which = hashval % THREADS;
  hashval = hashval/THREADS;
  hashtable = tables + which;
  memory_heap = heaps + which;

  int64_t pos = memory_heap->posInHeap;

  shared kmer_t *next_empty_kmer = (shared kmer_t*) (((shared char*) memory_heap->heap) + pos*sizeof(kmer_t));
  
  /* Add the contents to the appropriate kmer struct in the heap */
  upc_memput(next_empty_kmer->kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));

  next_empty_kmer->l_ext = left_ext;
  next_empty_kmer->r_ext = right_ext;
  
  /* Fix the next pointer to point to the appropriate kmer struct */
  next_empty_kmer->next = hashtable->table[hashval].head;
  /* Fix the head pointer of the appropriate bucket to point to the current kmer */
  hashtable->table[hashval].head = next_empty_kmer;
  
  /* Increase the heap pointer */
  memory_heap->posInHeap++;
  
  return next_empty_kmer;
}

/* Adds a k-mer in the start list by using the memory heap (the k-mer was "just added" in the memory heap at position posInHeap - 1) 
void addKmerToStartList(memory_heap_t *memory_heap, start_kmer_t **startKmersList)
{
  start_kmer_t *new_entry;
  kmer_t *ptrToKmer;
  
  int64_t prevPosInHeap = memory_heap->posInHeap - 1;
  ptrToKmer = &(memory_heap->heap[prevPosInHeap]);
  new_entry = (start_kmer_t*) malloc(sizeof(start_kmer_t));
  new_entry->next = (*startKmersList);
  new_entry->kmerPtr = ptrToKmer;
  (*startKmersList) = new_entry;
}
*/


const unsigned char *sample="AGCTAGCTGTAGCTGTAGC";
const unsigned char *left = "F";
const unsigned char *right = "F";

int main()
{
  shared memory_heap_t *mem;
  shared hash_table_t *tab = upc_create_hash_table(200,&mem);
  
  if(MYTHREAD==0)
    fprintf(stderr,"tab[%d].size: %d\ntab[%d].size: %d\ntab[%d].size: %d\n",0,tab[0].size,1,tab[1].size,2,tab[2].size);

  upc_barrier;

  for(int i=0;i<THREADS;i++)
  {
    if(i==MYTHREAD)
    {
     fprintf(stderr,"tab[%d].table: 0x%lx\ntab[%d].table: 0x%lx\ntab[%d].table: 0x%lx\n\n",
         0,(long int)tab[0].table,1,(long int)tab[1].table,2,(long int)tab[2].table);
    }
    sleep(1);
  }

  shared kmer_t* added;
  if(MYTHREAD==0)
  {
    added = add_kmer(tab, mem, sample, *left, *right);
  }

  upc_barrier;

  shared kmer_t* srced = lookup_kmer_upc(tab, mem, sample);
  fprintf(stderr,"I am %d: 0x%lx 0x%lx\n",MYTHREAD,(long int)added,(long int)srced);

  return 0;
}
#endif
