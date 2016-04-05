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
  shared memory_heap_t *mem = (shared memory_heap_t*) upc_all_alloc(1,sizeof(memory_heap_t));
  shared hash_table_t *tab = upc_create_hash_table(tsize,mem);
  
  if(MYTHREAD==0)
    fprintf(stderr,"tab[%d].size: %d\ntab[%d].size: %d\ntab[%d].size: %d\n",0,tab[0].size,1,tab[1].size,2,tab[2].size);

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

  shared kmer_t* added;
  unsigned char *s; 
  s = (unsigned char*) sample + MYTHREAD*0;
  unsigned char buf[KMER_LENGTH+1];
  buf[KMER_LENGTH] = (unsigned char) 0;

  upc_barrier;

  for(int i=0; i < tsize ; i++)
  {
    added = add_kmer(tab, mem, s + (i), *left, *right);
    memcpy(buf,s + (i),KMER_LENGTH);
    fprintf(stderr,"(%d) thread %d added %s at 0x%lx\n",i,MYTHREAD,buf,(long int)added);
  }

  fprintf(stderr,"MYTHREAD %d all kmers added!\n",MYTHREAD);
  upc_barrier;

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
