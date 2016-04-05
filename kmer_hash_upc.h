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
shared hash_table_t *upc_create_hash_table(uint64_t nEntries, shared memory_heap_t **memory_heap)
{
  shared hash_table_t *result;
  shared bucket_t *global_tables; 
  shared bucket_t *heaps;

  int64_t n_buckets = nEntries * LOAD_FACTOR;
  int64_t my_size = (n_buckets + THREADS - 1)/THREADS;
  int64_t heap_size = ((nEntries * LOAD_FACTOR + THREADS - 1)/THREADS);

  fprintf(stderr,"nEntries %ld my_size %ld heap_size %ld\n",nEntries,my_size,heap_size);

  if(n_buckets < 0 || my_size < 0 || heap_size <0)
  {
    fprintf(stderr,"sizes are zero! check yourself before you wreck yourself.\n");
    exit(1);
  }

  result = upc_all_alloc(THREADS,sizeof(hash_table_t));
  global_tables = (shared bucket_t*) upc_all_alloc( THREADS, my_size * sizeof(bucket_t) );
  // upc_memset(global_tables+MYTHREAD*my_size,0,my_size*sizeof(bucket_t));
  (result+MYTHREAD)->size = my_size;
  (result+MYTHREAD)->table = (shared bucket_t*) &global_tables[MYTHREAD*my_size];
  (result+MYTHREAD)->write_lock = upc_global_lock_alloc();
  for(int i=0;i<my_size;i++)
  {
    (result+MYTHREAD)->table[i].head = NULL;
  }


  *memory_heap = (shared memory_heap_t*) upc_all_alloc(THREADS,sizeof(memory_heap_t));
  heaps = (shared bucket_t *) upc_all_alloc( THREADS, heap_size * sizeof(kmer_t) );

  memory_heap[0][MYTHREAD].heap = (shared kmer_t*) (((shared char*) heaps) + (MYTHREAD * heap_size) * sizeof(kmer_t));
  memory_heap[0][MYTHREAD].posInHeap = 0;
  memory_heap[0][MYTHREAD].write_lock = upc_global_lock_alloc();

  /*
  fprintf(stderr,"size of memory_heap[0][i] %d size of ptr %d\n",sizeof(memory_heap[0][MYTHREAD]),sizeof(void*));
  fprintf(stderr,"total heap size %d\n",THREADS * heap_size * sizeof(kmer_t));
  for(int i=0;i<THREADS;i++)
  {
    fprintf(stderr,"heap for thread %d at 0x%ld\n",i,memory_heap[0][i].heap);
    fprintf(stderr,"DIFF thread %d at 0x%ld\n",i,(shared char*)memory_heap[0][i].heap - (shared char*)memory_heap[0][0].heap);
  }
  */

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
    // fprintf(stdout,"%7d bucket value \n",hashval % hashtable_size);
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
  char buf[KMER_LENGTH+1];
  buf[KMER_LENGTH] = '\0';
  memcpy(buf,kmer,KMER_LENGTH);

  fprintf(stderr,"THREAD %d are we here or what? %s\n",MYTHREAD,buf);
  char packedKmer[KMER_PACKED_LENGTH];
  packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);

  memcpy(buf,kmer,KMER_LENGTH);
  fprintf(stderr,"THREAD %d packing? %s\n",MYTHREAD,buf);

  int64_t hashval = hashkmer(hashtable->size * THREADS, (char*) packedKmer);
  int64_t which = hashval % THREADS;


  hashval = hashval / THREADS;
  hashtable += which;

  shared kmer_t *result;
  kmer_t res;

  fprintf(stderr,"THREAD %d hash? %d %d %s\n",MYTHREAD,hashval,which,buf);
  result = hashtable->table[hashval].head;
  fprintf(stderr,"THREAD %d we are here now\n",MYTHREAD);
  unsigned char cmp[KMER_PACKED_LENGTH+1];
  cmp[KMER_PACKED_LENGTH] = (unsigned char) 0;

  if(result==NULL)
    fprintf(stderr,"THREAD %d wtf why is result null\n",MYTHREAD);
  else
  {
    fprintf(stderr,"THREAD %d result is 0x%ld\n",MYTHREAD,result);
  }
  
  int isnull;
  while(result!=NULL) {
    fprintf(stderr,"THREAD %d while loop is starting\n",MYTHREAD);
    /*
    upc_memget(&cmp,result->kmer,KMER_PACKED_LENGTH);
    */
    upc_memget(&res,result,sizeof(kmer_t));
    memcpy(&cmp,res.kmer,KMER_PACKED_LENGTH);
    fprintf(stderr,"THREAD %d memcpy done\n",MYTHREAD);
    fprintf(stderr,"THREAD %d packedKmer %s || cmp %s || str %s\n",MYTHREAD,packedKmer,cmp,buf);
    isnull = 0;
    if(res.next==NULL)
    {
      isnull = 1;
      fprintf(stderr,"THREAD %d result is null apparently!\n",MYTHREAD);
    }
    // fprintf(stderr,"THREAD %d result %lx result->next %lx next is null? %d isnull %d\n",MYTHREAD,
    //     result,result->next,result->next == NULL);
    fprintf(stderr,"THREAD %d result->next is null? %d\n",MYTHREAD,result->next==NULL);
    fprintf(stderr,"THREAD %d isnull? %d\n",MYTHREAD,isnull);
    if( memcmp(packedKmer, cmp, KMER_PACKED_LENGTH) == 0 ) {
      memcpy(buf,kmer,KMER_LENGTH);
      fprintf(stderr,"THREAD %d ok! we return result now! %s %s\n",MYTHREAD,cmp,buf);
      return result;
    }
    else
    {
      fprintf(stderr,"THREAD %d whoops! that's not it! %c%c%c%c%c %c%c%c%c%c %d %d %d %d %d\n",MYTHREAD,
          packedKmer[0],packedKmer[1],packedKmer[2],packedKmer[3],packedKmer[4],
          cmp[0],cmp[1],cmp[2],cmp[3],cmp[4],
          packedKmer[0] == cmp[0],packedKmer[1] == cmp[1],packedKmer[2] == cmp[2],packedKmer[3] == cmp[3],packedKmer[4] == cmp[4]);
      fprintf(stderr,"THREAD %d whoops! that's not it! where %s {{{{%d}}}} %s 0x%ld\n",
          MYTHREAD,cmp,memcmp(packedKmer, cmp, KMER_PACKED_LENGTH),packedKmer,res.next);
      fprintf(stderr,"THREAD %d whoops! that's not it! 0x%ld\n",MYTHREAD,res.next);
    }

    if((isnull && !(result->next == NULL)) || (!isnull && result->next==NULL))
    {
      fprintf(stderr,"THREAD %d wtf isnull and result->next == NULL doesn't match LOL\n",MYTHREAD);
    }

    if(isnull)
    {
      fprintf(stderr,"THREAD %d next is NULL! (wtf?) kmer %s\n",MYTHREAD,cmp);
      break;
    }

    fprintf(stderr,"THREAD %d what the fuck is result: 0x%ld\n",MYTHREAD,&result->next);
    result = (shared kmer_t*)res.next;
    fprintf(stderr,"THREAD %d what the fuck is result: 0x%ld\n",MYTHREAD,result);
    fprintf(stderr,"THREAD %d what the fuck is result: %s %s 0x%ld\n",MYTHREAD,cmp,packedKmer,res.next);
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

  upc_lock(hashtable->write_lock);
  while(1)
  {
    if(!upc_lock_attempt(memory_heap->write_lock))
    {
      upc_unlock(hashtable->write_lock);
      upc_lock(hashtable->write_lock);
    }
    else break;
  }

  int64_t pos = memory_heap->posInHeap;
  shared kmer_t *next_empty_kmer = (shared kmer_t*) (((shared char*) memory_heap->heap) + pos*sizeof(kmer_t));
  
  /* Add the contents to the appropriate kmer struct in the heap */
  upc_memput(next_empty_kmer->kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
  next_empty_kmer->l_ext = left_ext;
  next_empty_kmer->r_ext = right_ext;
  
  /* Fix the next pointer to point to the appropriate kmer struct */
  next_empty_kmer->next = (shared void*)hashtable->table[hashval].head;
  /* Fix the head pointer of the appropriate bucket to point to the current kmer */
  hashtable->table[hashval].head = next_empty_kmer;


  /*
  char buf[KMER_LENGTH+1];
  memcpy(buf,kmer,KMER_LENGTH);
  buf[KMER_LENGTH] = '\0';
  // fprintf(stderr,"THREAD %d adding to table (which) %d at posInHeap %d: loc 0x%lx %s THREAD %d which %d pos %2d\n",
  //     MYTHREAD,which,memory_heap->posInHeap,(long int) next_empty_kmer,buf,MYTHREAD,which,pos);
  fprintf(stdout,"THREAD %d added kmer %s next at 0x%ld\n",MYTHREAD,buf,next_empty_kmer->next);
  */
  
  /* Increase the heap pointer */
  memory_heap->posInHeap++;
  fprintf(stderr,"posInHeap %10d on %d\n",memory_heap->posInHeap,which);

  upc_unlock(memory_heap->write_lock);
  upc_unlock(hashtable->write_lock);
  
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
