#include <stdint.h>

#ifndef UTILH
#define UTILH


// KMER_LENGTH
// since we pack 1 base -> 2 bits, we have a total of KMER_LENGTH * 2 / 8 = KMER_LENGTH/4 bytes, rounded up.
typedef struct {
	uint8_t* data;
	int length;
} kmer_t;

typedef struct {
	kmer_t data;
	uint8_t ext;
} kentry_t;
#endif
