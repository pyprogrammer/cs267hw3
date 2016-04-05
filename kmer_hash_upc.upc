#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include <upc.h>

#include "kmer_hash_upc.h"

const unsigned char *sample="AGCTAGCTGTAGCTGTAGCACGATGCTAGCATGCTAGCAGCTGTCATGCTGAGTCGTAGTATTATATGCTCGCGCGCGATCGTACGTACTGCTATCGTAGCTAGCTGAGCGCTAGCTAGTCGCTACGTAGTAGCTCGTATCGTAGCTGTCGATGCTGCCGCCGTATATCGATGCATATATTAGCGCGTCTCAGCGCGGCTCTCGCTGTAGTACGATGCTCGTGCGCGATACGT";
const unsigned char left = 'F';
const unsigned char right = 'F';

typedef shared int* sint_t;
		
		
int main()
{
//  shared memory_heap_t *mem;
//  shared hash_table_t *tab = upc_create_hash_table(4514197,&mem);
//  
//  if(MYTHREAD==0)
//    fprintf(stderr,"tab[%d].size: %d\ntab[%d].size: %d\ntab[%d].size: %d\n",0,tab[0].size,1,tab[1].size,2,tab[2].size);
//
//  for(int i=0;i<90;i++) fprintf(stderr,"thread %d i %d\n",MYTHREAD,i);
//
//  fprintf(stderr,"\n\n\n\n\n");
//
//  upc_barrier;
//  sleep(1);

  
//  shared int* sp;
  sint_t arr = upc_all_alloc(THREADS, THREADS*sizeof(int));
  shared sint_t* directory = upc_all_alloc(THREADS, THREADS*sizeof(sint_t));
  directory[MYTHREAD] = &arr[MYTHREAD];
  *(directory[MYTHREAD]) = MYTHREAD;
  
  
  {
	  for (int i = 0; i < THREADS; i++)
	  {
		  printf("PRINTING FROM: %d\tTHREAD: %d\tN:%d\n", MYTHREAD, i, arr[i]);
	  }
  }
  
  
  

  return 0;
}
