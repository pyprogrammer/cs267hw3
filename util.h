#ifndef UTILH
#define UTILH

#include <stdint.h>
#include <packingDNAseq.h>

#include "kmer_hash_upc.h"

#define HASH_DEPTH 2 // 2 chars deep for the hash.

typedef struct entrylistel {
	struct entrylistel* next;
	struct entrylistel* prev;
	shared kmer_t* entry;
	char left_ext;
	char right_ext;
} entrylistel_t;

typedef struct {
	entrylistel_t* start;
	entrylistel_t* end;
	int size;
} entrylist_t;

typedef struct {
	char* buf;
	int size;
	int bufsize;
} stringbuffer_t;

void init_stringbuffer(stringbuffer_t* sb, int size)
{
	sb->bufsize = size+1;
	sb->size = 0;
	sb->buf = (char*) malloc(sizeof(char) * sb->bufsize);
	memset(sb->buf, 0, sizeof(char) * sb->bufsize);
}

void expand_buffer(stringbuffer_t* sb)
{
	sb->buf = realloc(sb->buf, sb->bufsize*2-1);
	sb->bufsize = sb->bufsize*2-1;
}

void insert_buffer(stringbuffer_t* sb, char c)
{
	if (sb->size == sb->bufsize-1)
	{
		expand_buffer(sb);
	}
	sb->buf[sb->size++] = c;
	sb->buf[sb->size] = '\0';
}

void init_list(entrylist_t* elist)
{
	elist->start = NULL;
	elist->end = NULL;
	elist->size = 0;
}

void append_list(entrylist_t* elist, shared kmer_t* entry, char left_ext, char right_ext)
{
	elist->size++;
	entrylistel_t* el = (entrylistel_t*) malloc(sizeof(entrylistel_t));
	el->prev = NULL;
	el->next = NULL;
	el->entry = entry;
	el->left_ext = left_ext;
	el->right_ext = right_ext;
	if (elist->end != NULL)
	{
		el->prev = elist->end;
		elist->end->next = el;
	}
	elist->end = el;
	if (elist->start == NULL) elist->start = el;
}

void pop_list(entrylist_t* elist, shared kmer_t** ret, char* left_ext, char* right_ext)
{
	entrylistel_t* last = elist->end;
	elist->end = elist->end->prev;
	*ret = last->entry;
	*left_ext = last->left_ext;
	*right_ext = last->right_ext;
	free(last);
	elist->size--;
}

void shift_into_kmer(shared kmer_t* current, char* dest, char append)
{
	unsigned char cpy[KMER_LENGTH];
	upc_memget(cpy, current->kmer, KMER_LENGTH);
	unpackSequence(cpy, dest, KMER_LENGTH);
	dest[KMER_LENGTH] = append;
//	packSequence(buf+1, dest->kmer, KMER_LENGTH);
}

void print_kmer(shared kmer_t* kmer)
{
	unsigned char cpy[KMER_LENGTH+1];
	unsigned char packed[KMER_PACKED_LENGTH];
	cpy[KMER_LENGTH] = 0;
	upc_memget(packed, kmer->kmer, KMER_LENGTH);
	unpackSequence(packed, cpy, KMER_LENGTH);
	fprintf(stderr, "KMER: %s\n", cpy);
}
#endif
