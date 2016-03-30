#include <stdint.h>

#ifndef UTILH
#define UTILH

#define HASH_DEPTH 2 // 2 chars deep for the hash.

typedef struct entrylistel {
	shared struct entrylistel* next;
	shared struct entrylistel* prev;
	shared kmer_t* entry;
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

void append_list(shared entrylist_t* elist, shared kmer_t* entry)
{
	shared entrylistel_t* el = (shared entrylistel_t*) upc_alloc(sizeof(entrylistel_t));
	el->prev = NULL;
	el->next = NULL;
	el->entry = entry;
	if (elist->end != NULL)
	{
		el->prev = elist->end;
		elist->end->next = el;
	}
	elist->end = el;
	if (elist->start == NULL) elist->start = el;
}

shared kmer_t* pop_list(shared entrylist_t* elist)
{
	shared entrylistel_t* last = elist->end;
	elist->end = elist->end->prev;
	shared kmer_t* retval = last->entry;
	upc_free(last);
	return retval;
}
#endif
