#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include <upc.h>

#include "kmer_hash_upc.h"

/*
#define TEST_SIZE 4514197
*/
#define TEST_SIZE 1000

const unsigned char *sample="AGCTAGCTGTAGCTGTAGCACGATGCTAGCATGCTAGCAGCTGTCATGCTGAGTCGTAGTATTATATGCTCGCGCGCGATCGTACGTACTGCTATCGTAGCTAGCTGAGCGCTAGCTAGTCGCTACGTAGTAGCTCGTATCGTAGCTGTCGATGCTGCCGCCGTATATCGATGCATATATTAGCGCGTCTCAGCGCGGCTCTCGCTGTAGTACGATGCTCGTGCGCGATACGT";
const unsigned char *left = "F";
const unsigned char *right = "F";
const unsigned char *whoops="AAAAAAAAAAAAAAAAAAA\0";

int main()
{
  char x = 0xfd;
  hex_dump((void*)&x,1);

  for(int i=0;i<20;i++)
  {
    hex_dump((void*)sample,i);
  }

  shared memory_heap_t *mem;
  shared hash_table_t *tab = upc_create_hash_table(TEST_SIZE,&mem);
  
  if(MYTHREAD==0)
    fprintf(stderr,"tab[%d].size: %d\ntab[%d].size: %d\ntab[%d].size: %d\n",0,tab[0].size,1,tab[1].size,2,tab[2].size);

  fprintf(stderr,"\n\n\n\n\n");

  upc_barrier;
  sleep(1);

  /*
  for(int i=0;i<THREADS;i++)
  {
    if(i==MYTHREAD)
    {
     fprintf(stderr,"tab[%d].table: 0x%lx\ntab[%d].table: 0x%lx\ntab[%d].table: 0x%lx\n\n",
         0,(long int)tab[0].table,1,(long int)tab[1].table,2,(long int)tab[2].table);
    }
    sleep(1);
  }
  */

  shared kmer_t* added;
  unsigned char *s; 
  s = (unsigned char*) sample + MYTHREAD*0;
  unsigned char buf[KMER_LENGTH+1];
  buf[KMER_LENGTH] = (unsigned char) 0;

  kmer_t kmerbuf;

  upc_barrier;
  sleep(1);

  for(int i=0; i < TEST_SIZE ; i++)
  {
    added = add_kmer(tab, mem, s + (i%200), *left, *right);
    upc_memget(&kmerbuf,added,sizeof(kmer_t));
    // hex_dump((void*)&kmerbuf,sizeof(kmer_t));
    // memcpy(buf,s + (i%200),KMER_LENGTH);
    // fprintf(stderr,"(%d) thread %d added %s at 0x%lx\n",i,MYTHREAD,buf,(long int)added);
  }

  fprintf(stderr,"MYTHREAD %d all kmers added!\n",MYTHREAD);
  upc_barrier;
  sleep(1);

  shared kmer_t* srced;

  if(MYTHREAD==1)
  /*
  for(int i=0;i<TEST_SIZE;i++)
  {
    memcpy(buf,(unsigned char*) sample + (i%200),KMER_LENGTH);
    fprintf(stderr,"Thread %d looking at %s (len %d) i %d\n",MYTHREAD,buf,KMER_LENGTH,i);
    srced = lookup_kmer_upc(tab, mem, sample+(i%200));
    // fprintf(stderr,"I am %d: 0x%lx\n",MYTHREAD,(long int)srced);
  }
  */

  for(int i=TEST_SIZE/2-100;i<TEST_SIZE/2+100;i++)
  {
    memcpy(buf,(unsigned char*) sample + (i%200),KMER_LENGTH);
    // fprintf(stderr,"Thread %d looking at %s (len %d)     ",MYTHREAD,buf,KMER_LENGTH);
    srced = lookup_kmer_upc(tab, mem, sample+(i%200));
    // fprintf(stderr,"I am %d: 0x%lx\n",MYTHREAD,(long int)srced);
  }

  // srced = lookup_kmer_upc(tab, mem, whoops);


  upc_barrier;
  
  fprintf(stderr,"THREAD %d is done.\n",MYTHREAD);
  fprintf(stdout,"THREAD %d is done.\n",MYTHREAD);

  return 0;
}
