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

upc_lock_t *global_lock;


/* Creates hashtable for UPC */
/* Whenever something is called on a hash_table_t or memory_heap_t, instead call it on
 * hash_dir_t[MYTHREAD] and mem_dir_t[MYTHREAD]
 */
shared hash_table_t *upc_create_hash_table(int64_t nEntries, shared memory_heap_t *memory_heap)
{
  shared hash_table_t *result, *mytable;
  int64_t n_buckets = (uint64_t)(nEntries * LOAD_FACTOR);
  if(THREADS == 1)
    n_buckets = nEntries;

  result = (shared hash_table_t*) upc_all_alloc(THREADS,sizeof(hash_table_t));
  mytable = result + MYTHREAD;
  mytable->table = (shared bucket_t*) upc_all_alloc(n_buckets,sizeof(bucket_t));
  for(int i=0;i<n_buckets;i++)
  {
    shared bucket_t *b = mytable->table + i;
    fprintf(stderr,"THREAD %d doing dummy op\n",MYTHREAD);
    b->lock = NULL;
    fprintf(stderr,"THREAD %d making a lock???\n",MYTHREAD);
    b->lock = upc_all_lock_alloc();
    fprintf(stderr,"THREAD %d it worked!!!\n",MYTHREAD);
  }
  mytable->table = mytable->table + MYTHREAD;
  mytable = result + MYTHREAD;
  mytable->size = (n_buckets + THREADS - 1)/THREADS;

  shared memory_heap_t *myheap = memory_heap + MYTHREAD;
  myheap->heap = (shared kmer_t*) upc_all_alloc( ((n_buckets + THREADS - 1)/THREADS)*THREADS,sizeof(kmer_t));
  myheap->heap = myheap->heap + MYTHREAD;
  myheap->posInHeap = 0;
  myheap->lock = upc_all_lock_alloc();

  return result;
}

/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
  unsigned long hashval;
  hashval = 5381;
  for(int i = 0; i < size; i++) {
    hashval = seq[i] + (hashval << 5) + hashval;
  }
  return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
  return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

const char *important = "ATCTCGGCTTTGATAAACC";
int iwhich = -1;
int ipos = -1;
/* Looks up a kmer in the hash table and returns a pointer to that entry */
shared kmer_t* lookup_kmer_upc(shared hash_table_t *hashtable, shared memory_heap_t *memory_heap, const unsigned char *kmer)
{
  init_LookupTable();
  char packedKmer[KMER_PACKED_LENGTH];
  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
  int64_t hashval = hashkmer(hashtable->size, (char*) packedKmer);
  int64_t which = hashval % THREADS;
  hashval = hashval / THREADS;

  hashtable = hashtable + which;
  memory_heap = memory_heap + which;

  shared bucket_t *bucket = hashtable->table + hashval * THREADS;

  shared kmer_t *result;
  
  result = bucket->head;
  unsigned char cmp[KMER_PACKED_LENGTH+1];
  cmp[KMER_PACKED_LENGTH] = (unsigned char) 0;
  unsigned char buf1[KMER_LENGTH + 1];
  unsigned char buf2[KMER_LENGTH + 1];
  buf1[KMER_LENGTH] = '\0';
  buf2[KMER_LENGTH] = '\0';
  
  while(result!=NULL) {
    if( memcmp(important,kmer,KMER_LENGTH) == 0)
    {
      fprintf(stderr,"THREAD %d looking for important string at which %d pos %d next %d\n",
          MYTHREAD,which,result->pos,result->next_pos);
    }
    // fprintf(stderr,"THREAD %d while loop starting\n",MYTHREAD);
    upc_memget(&cmp,result->kmer,KMER_PACKED_LENGTH);
    // fprintf(stderr,"THREAD %d packedKmer %s || cmp %s\n",MYTHREAD,packedKmer,cmp);
    if( memcmp(packedKmer, cmp, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
      // fprintf(stdout,"+%8d\n",result->pos);
      return result;
    }
    // fprintf(stderr,"THREAD %d memory getting\n",MYTHREAD);
    unpackSequence(cmp,(unsigned char*) buf1,KMER_LENGTH);
    unpackSequence(packedKmer,(unsigned char*) buf2,KMER_LENGTH);
    // fprintf(stderr,"THREAD %d comparing %s and %s\n",MYTHREAD,buf1,buf2);
    // memcpy(buf2,kmer,KMER_LENGTH);
    // fprintf(stderr,"THREAD %d comparing %s and %s\n",MYTHREAD,cmp,buf2);

    int n_pos = result->next_pos;
    // fprintf(stderr,"THREAD %d no success, looking at pos %d next\n",MYTHREAD,n_pos);
    if(n_pos == -1) break;
    result = memory_heap->heap + n_pos * THREADS;
  }
  fprintf(stderr, "THREAD %d DID NOT FIND\n",MYTHREAD);
  exit(-1);
  return NULL;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
shared kmer_t* add_kmer(shared hash_table_t *hashtable, shared memory_heap_t *memory_heap,
    const unsigned char *kmer, char left_ext, char right_ext)
{
  /* Pack a k-mer sequence appropriately */
  char packedKmer[KMER_PACKED_LENGTH];

  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);

  int64_t hashval = hashkmer(hashtable->size, (char*) packedKmer);
  int64_t which = hashval % THREADS;
  hashval = hashval / THREADS;

  hashtable = hashtable + which;
  memory_heap = memory_heap + which;


  upc_lock(memory_heap->lock);
    int64_t pos = memory_heap->posInHeap;
    /* Increase the heap pointer */
    memory_heap->posInHeap++;
  upc_unlock(memory_heap->lock);


  // shared kmer_t *next_empty_kmer = (shared kmer_t*) (((shared char*) memory_heap->heap) + pos*sizeof(kmer_t));
  shared kmer_t *next_empty_kmer = memory_heap->heap + pos * THREADS;
  shared bucket_t *bucket = hashtable->table + hashval * THREADS;


  upc_lock(bucket->lock);
    if(bucket->head == NULL)
    {
      next_empty_kmer->next_pos = -1;
      next_empty_kmer->next_which = -1;
    }
    else
    {
      next_empty_kmer->next_pos = bucket->head->pos;
      next_empty_kmer->next_which = bucket->head->which;
    }
    next_empty_kmer->pos = pos;
    next_empty_kmer->which = which;
    /* Fix the head pointer of the appropriate bucket to point to the current kmer */
    bucket->head = next_empty_kmer;
  upc_unlock(bucket->lock);


  /* Add the contents to the appropriate kmer struct in the heap */
  upc_memput(next_empty_kmer->kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
  next_empty_kmer->l_ext = left_ext;
  next_empty_kmer->r_ext = right_ext;
  
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

#endif
