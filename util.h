#include <stdint.h>
#include <packingDNAseq.h>

#ifndef UTILH
#define UTILH

#define HASH_DEPTH 2 // 2 chars deep for the hash.

typedef struct entrylistel {
	shared struct entrylistel* next;
	shared struct entrylistel* prev;
	shared kmer_t* entry;
	char left_ext;
	char right_ext;
} entrylistel_t;

typedef struct {
	shared entrylistel_t* start;
	shared entrylistel_t* end;
} entrylist_t;

void init_list(shared entrylist_t* elist)
{
	elist->start = NULL;
	elist->end = NULL;
}

void append_list(shared entrylist_t* elist, shared kmer_t* entry, char left_ext, char right_ext)
{
	shared entrylistel_t* el = (shared entrylistel_t*) upc_alloc(sizeof(entrylistel_t));
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

void pop_list(shared entrylist_t* elist, kmer_t** ret, char* left_ext, char* right_ext)
{
	shared entrylistel_t* last = elist->end;
	elist->end = elist->end->prev;
	*ret = last->entry;
	*left_ext = last->left_ext;
	*right_ext = last->right_ext;
	upc_free(last);
}

void shift_into_kmer(shared kmer_t* current, kmer_t* dest, char append)
{
	unsigned char buf[KMER_LENGTH+1];
	unpackSequence(current->kmer, buf, KMER_LENGTH);
	buf[KMER_LENGTH] = append;
	packSequence(buf+1, dest->kmer, KMER_LENGTH);
}
#endif
