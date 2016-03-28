#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"

#include "util.h"
#include "kmer_hash_upc.h"

int main(int argc, char *argv[]){

  /** Declarations **/
  double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
  char *input_name;
  int64_t nKmers, total_chars_to_read;
  
  
  /** Read input **/
  upc_barrier;
  inputTime -= gettime();
  ///////////////////////////////////////////
  // Your code for input file reading here //
  ///////////////////////////////////////////

  input_name = argv[1];
  nKmers = getNumKmersInUFX(input_name);
  hash_table_t *hashtable;
  memory_heap_t memory_heap;

  hash_dir_t * hash_dir;
  mem_dir_t * mem_dir;

  hashtable = create_hash_table(nKmers, &memory_heap);

  total_chars_to_read = nKmers *LINE_SIZE;


  upc_barrier;
  inputTime += gettime();
  






  /** Graph construction **/
  constrTime -= gettime();
  ///////////////////////////////////////////
  // Your code for graph construction here //
  ///////////////////////////////////////////
  upc_barrier;
  constrTime += gettime();
  










  /** Graph traversal **/
  traversalTime -= gettime();
  ////////////////////////////////////////////////////////////
  // Your code for graph traversal and output printing here //
  // Save your output to "pgen.out"                         //
  ////////////////////////////////////////////////////////////
  upc_barrier;
  traversalTime += gettime();











  
  /** Print timing and output info **/
  /***** DO NOT CHANGE THIS PART ****/
  if(MYTHREAD==0){
  	printf("%s: Input set: %s\n", argv[0], argv[1]);
  	printf("Number of UPC threads: %d\n", THREADS);
  	printf("Input reading time: %f seconds\n", inputTime);
  	printf("Graph construction time: %f seconds\n", constrTime);
  	printf("Graph traversal time: %f seconds\n", traversalTime);
  }
  return 0;

}
