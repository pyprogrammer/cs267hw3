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
  int64_t ptr = LINE_SIZE*MYTHREAD;
  int64_t cur_chars_read = 0;
  char* input_UFX_name = argv[1];
  unsigned char* working_buffer;
  FILE *inputFile;
  init_LookupTable();
  
  
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
  working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char));
  inputFile = fopen(input_UFX_name, "r");
  cur_chars_read = fread(working_buffer, sizeof(unsigned char),total_chars_to_read , inputFile);
  fclose(inputFile);
  
  shared entrylist_t* entrylist = (shared entrylist_t*) upc_alloc(sizeof(entrylist_t));
  init_list(entrylist);
  
  shared entrylist_t* startlist = (shared entrylist_t*) upc_alloc(sizeof(entrylist_t));
  init_list(startlist);
  
  while (ptr < cur_chars_read) {
    /* working_buffer[ptr] is the start of the current k-mer                */
    /* so current left extension is at working_buffer[ptr+KMER_LENGTH+1]    */
    /* and current right extension is at working_buffer[ptr+KMER_LENGTH+2]  */
	char left_ext = (char) working_buffer[ptr+KMER_LENGTH+1];
	char right_ext = (char) working_buffer[ptr+KMER_LENGTH+2];
	
	/* Add k-mer to hash table */
	shared kmer_t* location = add_kmer_upc(hashtable, &memory_heap, &working_buffer[ptr], left_ext, right_ext); // puts in unpacked
	
	/* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
	if (left_ext == 'F') {
		append_list(startlist, location, left_ext, right_ext);
	}
    append_list(entrylist, location, left_ext, right_ext);
    
    

    /* Move to the next k-mer in the input working_buffer */
    ptr += LINE_SIZE*THREADS;
  }


  upc_barrier;
  inputTime += gettime();
  
  shared kmer_t* curr;
  char left_ext;
  char right_ext;
  kmer_t newkmer;
  while(entrylist->end != NULL)
  {
	  pop_list(entrylist, &curr, &left_ext, &right_ext);
	  shift_into_kmer(curr, &newkmer, right_ext);
	  shared kmer_t* next = find_kmer_upc(hashtable, &memory_heap, &newkmer);
	  curr->next = next;
  }




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
