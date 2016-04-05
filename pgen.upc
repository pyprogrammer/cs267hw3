#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>

#include "packingDNAseq.h"
/*
#include "kmer_hash.h"
*/

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
  char cur_contig[MAXIMUM_CONTIG_SIZE];
  char unpackedKmer[KMER_LENGTH+1];
  unpackedKmer[KMER_LENGTH] = '\0';
  shared kmer_t* cur_kmer_ptr;
  int64_t posInContig;
  int64_t contigID = 0;
  init_LookupTable();
  
  
  /** Read input **/
  upc_barrier;
  inputTime -= gettime();
  ///////////////////////////////////////////
  // Your code for input file reading here //
  ///////////////////////////////////////////

  input_name = argv[1];
  nKmers = getNumKmersInUFX(input_name);
  shared hash_table_t *hashtable;
  shared memory_heap_t *memory_heap = (shared memory_heap_t*) upc_all_alloc(1,sizeof(memory_heap_t));

  hashtable = upc_create_hash_table(nKmers, memory_heap);

  total_chars_to_read = nKmers *LINE_SIZE;
  working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char));
  inputFile = fopen(input_UFX_name, "r");
  cur_chars_read = fread(working_buffer, sizeof(unsigned char),total_chars_to_read , inputFile);
  fclose(inputFile);
  
  entrylist_t* entrylist = (entrylist_t*) malloc(sizeof(entrylist_t));
  init_list(entrylist);
  
  entrylist_t* startlist = (entrylist_t*) malloc(sizeof(entrylist_t));
  init_list(startlist);
  
  while (ptr < cur_chars_read) {
    /* working_buffer[ptr] is the start of the current k-mer                */
    /* so current left extension is at working_buffer[ptr+KMER_LENGTH+1]    */
    /* and current right extension is at working_buffer[ptr+KMER_LENGTH+2]  */
	char left_ext = (char) working_buffer[ptr+KMER_LENGTH+1];
	char right_ext = (char) working_buffer[ptr+KMER_LENGTH+2];
	
	/* Add k-mer to hash table */
	shared kmer_t* location = add_kmer(hashtable, memory_heap, working_buffer + ptr, left_ext, right_ext); // puts in unpacked
	
	/* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
	if (left_ext == 'F') {
		append_list(startlist, location, left_ext, right_ext);
	}
    append_list(entrylist, location, left_ext, right_ext);
    
    

    /* Move to the next k-mer in the input working_buffer */
    ptr += LINE_SIZE*THREADS;
  }
  fprintf(stderr,"done reading input\n");

  upc_barrier;
  inputTime += gettime();
  





  /** Graph construction **/
  constrTime -= gettime();
  ///////////////////////////////////////////
  // Your code for graph construction here //
  ///////////////////////////////////////////
  
  shared kmer_t* curr;
  char left_ext;
  char right_ext;
  char newkmer[KMER_LENGTH+2];
  newkmer[KMER_LENGTH+1] = '\0';
  while(entrylist->end != NULL)
  {
	  pop_list(entrylist, &curr, &left_ext, &right_ext);
	  if (right_ext == 'F') continue; // we done here
	  shift_into_kmer(curr, newkmer, right_ext);
	  shared kmer_t* next = lookup_kmer_upc(hashtable, memory_heap, newkmer+1);
      if(next==NULL)
      {
        fprintf(stderr,"THREAD %d wtf there's a null %s\n",MYTHREAD,newkmer+1);
      }
	  curr->next_kmer_pos = next->pos;
  }

  fprintf(stderr,"done linking\n");
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
  char outputFilename[32];
  sprintf(outputFilename, "pgen.%d.out", MYTHREAD);
  FILE* output = fopen(outputFilename, "w");
  for(int i=0;startlist->end != NULL;i++)
  {
	  pop_list(startlist, &cur_kmer_ptr, &left_ext, &right_ext);
	  unpackSequence((unsigned char*) cur_kmer_ptr->kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
	  /* Initialize current contig with the seed content */
	  memcpy(cur_contig ,unpackedKmer, KMER_LENGTH * sizeof(char));
	  posInContig = KMER_LENGTH;
	  right_ext = cur_kmer_ptr->r_ext;

	  /* Keep adding bases while not finding a terminal node */
	  while (right_ext != 'F') {
		cur_contig[posInContig] = right_ext;
		posInContig++;
		int next_pos = cur_kmer_ptr->next_kmer_pos;
		cur_kmer_ptr = memory_heap->heap + next_pos;
		//cur_kmer_ptr = cur_kmer_ptr->next;
		right_ext = cur_kmer_ptr->r_ext;
	  }

	  /* Print the contig since we have found the corresponding terminal node */
	  cur_contig[posInContig] = '\0';
	  fprintf(output,"%s\n", cur_contig);
	  if (i == 0) fprintf(stderr, "FIRST\n");
	  contigID++;
	  //totBases += strlen(cur_contig);

      if(!i%10000)
        fprintf(stderr,"iteration %d\n",i);

  }
  fprintf(stderr,"done travelling salesman\n");
  fclose(output);









  
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
