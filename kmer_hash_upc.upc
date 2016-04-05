#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include <upc.h>

#include "kmer_hash_upc.h"

#define TESTSIZE 4514197
#define tsize 800

const unsigned char * sample = "GCCGAAAGGGCCCTAGCGCTAAGGTCGTGATGGCGAGGTCGCGTCCACACAGCTTGAGGGAACTCGTGAAACACCGTACGTAGCACCTCGGAAGTAAGCTCGTCCACTCTAGTTCAAACGACAGGCACACGTATTACCGGTGGGACCAATGACGGTGTCACACGCATACTCCATCTATATCCAGCGCGCATTTGCGGATTAAATCAACTCTTCACTCCTGAATTTAAATCGCCGCCTTATGATTTAGTGTAGTCAACTTGACATGTTGCATGTTACTGATGGCGCTCGACCAAGTTGCGGTCTTTTAATTAGGCCTAGTCGGTGGATTAATTTTTGTAATATGACACCCGTTCTAGCCATGATAAACCTTTCTTAGGGCGTTGTTTTTAGGTGGGCATCCCTTTTAACGAGGGCTACAACTAAGATTACGTGGCCCTGTGGACCTAGTTCAAGCTCCTGAAATCGACTGCTTGGGTTGCCTGTAAACGAGCAGCCTAAACTAAAATTCTAATTTAGTCACTGCAGCACTCATACCGATGCGCTTTAAATGTGGCCCGTAAGGCTAGCCTGGAACGGACTACATCAACGCAGGGGCATCTCGGCTTTGATAAACCATATTGCGAAGAGCACCAAACCCGTTTCCAATGCTGTAGAAAGGCGCCCTTCAGGATTATCAGAAGGCGAGACTTACCCTTAATAGTATTCCCACTGTCGACTTCCTTGTTTGTGACGCTCGACGCGCACTAGTGCAAGAACTCCTTGCGACATCTTAGCACTCGAGATTTATTACCGGTGCGCTATACCGGTGCTAGAGCCCTAGATACGGTACGTAATATTGTCGTATCATAGGATAATGTGCCAGTTGAATTTCCACGATGTTTGGAATGTCAGCACAACAGGTTTGGAGGCCATGTATAAATGTTAGGTGAAGCCCAGTGCTACGTCCAGCTACCATTTGACCGGTAAACCAGTTCCCTCCCCGGCTTTGGTCGTTCGAGCG";
const unsigned char *left = "F";
const unsigned char *right = "F";

int main()
{
  fprintf(stderr,"lock check start\n");

  shared int *t = (shared int*) upc_all_alloc(1,sizeof(int));
  int i;
  upc_lock_t *shitty_lock;
  shitty_lock = upc_all_lock_alloc();
  upc_lock(shitty_lock);
  upc_memget(&i,t,sizeof(int));
  fprintf(stderr,"THREAD %d i %d\n",MYTHREAD,i);
  t[0] = i+1;
  upc_unlock(shitty_lock);

  upc_barrier;

  fprintf(stderr,"making a hash\n");

  shared memory_heap_t *mem = (shared memory_heap_t*) upc_all_alloc(1,sizeof(memory_heap_t));
  shared hash_table_t *tab = upc_create_hash_table(tsize*THREADS,mem);
  
  fprintf(stderr,"THREAD %d table size %d\n",MYTHREAD,tab->size);

  fprintf(stderr,"\n\n\n\n\n");

  upc_barrier;

  /*
  for(int i=0;i<THREADS;i++)
  {
    if(i==MYTHREAD)
    {
     fprintf(stderr,"tab[%d].table: 0x%lx\ntab[%d].table: 0x%lx\ntab[%d].table: 0x%lx\n\n",
         0,(long int)tab[0].table,1,(long int)tab[1].table,2,(long int)tab[2].table);
    }
  }
  */

  fprintf(stderr,"adding shit\n");

  shared kmer_t* added;
  unsigned char *s; 
  s = (unsigned char*) sample + MYTHREAD*0;
  unsigned char buf[KMER_LENGTH+1];
  buf[KMER_LENGTH] = (unsigned char) 0;

  upc_barrier;

  for(int i=0; i < tsize ; i++)
  {
    if(i%THREADS == MYTHREAD)
    {
      added = add_kmer(tab, mem, s + (i), *left, *right);
      memcpy(buf,s + (i),KMER_LENGTH);
      fprintf(stderr,"(%d) thread %d added %s at 0x%lx\n",i,MYTHREAD,buf,(long int)added);
    }
  }

  fprintf(stderr,"MYTHREAD %d all kmers added!\n",MYTHREAD);
  upc_barrier;

  fprintf(stderr,"finding shit\n");

  shared kmer_t* srced;

  for(int i=tsize-1;i>=0;i--)
  {
    memcpy(buf,(unsigned char*) sample + i,KMER_LENGTH);
    fprintf(stderr,"Thread %d looking at %s (len %d)     \n",MYTHREAD,buf,KMER_LENGTH);
    srced = lookup_kmer_upc(tab, mem, sample+i);
    // fprintf(stderr,"I am %d: 0x%lx\n",MYTHREAD,(long int)srced);
  }

  upc_barrier;

  fprintf(stderr,"we're done!\n");

  return 0;
}
