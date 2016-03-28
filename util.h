#include <stdint.h>

#ifndef UTILH
#define UTILH

#define HASH_DEPTH 2 // 2 chars deep for the hash.

// KMER_LENGTH
// since we pack 1 base -> 2 bits, we have a total of KMER_LENGTH * 2 / 8 = KMER_LENGTH/4 bytes, rounded up.
typedef struct {
	uint8_t data[KMER_PACKED_LENGTH]; // each uint8 has the first base as the highest 2 bits.
} pkmer_t;

typedef struct pkstruct {
	pkmer_t data;
	uint8_t ext; // we use high 4 for back ext, low 4 for forward ext.
	struct pkstruct* next;
	// 0-3 = ACGT, 4 = F
} pkentry_t;


uint8_t baseat(pkmer_t* kmer, int index)
{
	int access_offset = index / 4;
	int shift = 6 - (index % 4) * 2; // so remainder 0 gives a shift of 6 and remainder 3 gives a shift of 0.
	uint8_t section = kmer->data[access_offset];
	return (section >> shift) & 3;
}


char lookup_table[256][5] = {"AAAA", "AAAC", "AAAG", "AAAT", "AACA", "AACC", "AACG", "AACT", "AAGA", "AAGC", "AAGG", "AAGT", "AATA", "AATC", "AATG", "AATT", "ACAA", "ACAC", "ACAG", "ACAT", "ACCA", "ACCC", "ACCG", "ACCT", "ACGA", "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "ACTT", "AGAA", "AGAC", "AGAG", "AGAT", "AGCA", "AGCC", "AGCG", "AGCT", "AGGA", "AGGC", "AGGG", "AGGT", "AGTA", "AGTC", "AGTG", "AGTT", "ATAA", "ATAC", "ATAG", "ATAT", "ATCA", "ATCC", "ATCG", "ATCT", "ATGA", "ATGC", "ATGG", "ATGT", "ATTA", "ATTC", "ATTG", "ATTT", "CAAA", "CAAC", "CAAG", "CAAT", "CACA", "CACC", "CACG", "CACT", "CAGA", "CAGC", "CAGG", "CAGT", "CATA", "CATC", "CATG", "CATT", "CCAA", "CCAC", "CCAG", "CCAT", "CCCA", "CCCC", "CCCG", "CCCT", "CCGA", "CCGC", "CCGG", "CCGT", "CCTA", "CCTC", "CCTG", "CCTT", "CGAA", "CGAC", "CGAG", "CGAT", "CGCA", "CGCC", "CGCG", "CGCT", "CGGA", "CGGC", "CGGG", "CGGT", "CGTA", "CGTC", "CGTG", "CGTT", "CTAA", "CTAC", "CTAG", "CTAT", "CTCA", "CTCC", "CTCG", "CTCT", "CTGA", "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT", "GAAA", "GAAC", "GAAG", "GAAT", "GACA", "GACC", "GACG", "GACT", "GAGA", "GAGC", "GAGG", "GAGT", "GATA", "GATC", "GATG", "GATT", "GCAA", "GCAC", "GCAG", "GCAT", "GCCA", "GCCC", "GCCG", "GCCT", "GCGA", "GCGC", "GCGG", "GCGT", "GCTA", "GCTC", "GCTG", "GCTT", "GGAA", "GGAC", "GGAG", "GGAT", "GGCA", "GGCC", "GGCG", "GGCT", "GGGA", "GGGC", "GGGG", "GGGT", "GGTA", "GGTC", "GGTG", "GGTT", "GTAA", "GTAC", "GTAG", "GTAT", "GTCA", "GTCC", "GTCG", "GTCT", "GTGA", "GTGC", "GTGG", "GTGT", "GTTA", "GTTC", "GTTG", "GTTT", "TAAA", "TAAC", "TAAG", "TAAT", "TACA", "TACC", "TACG", "TACT", "TAGA", "TAGC", "TAGG", "TAGT", "TATA", "TATC", "TATG", "TATT", "TCAA", "TCAC", "TCAG", "TCAT", "TCCA", "TCCC", "TCCG", "TCCT", "TCGA", "TCGC", "TCGG", "TCGT", "TCTA", "TCTC", "TCTG", "TCTT", "TGAA", "TGAC", "TGAG", "TGAT", "TGCA", "TGCC", "TGCG", "TGCT", "TGGA", "TGGC", "TGGG", "TGGT", "TGTA", "TGTC", "TGTG", "TGTT", "TTAA", "TTAC", "TTAG", "TTAT", "TTCA", "TTCC", "TTCG", "TTCT", "TTGA", "TTGC", "TTGG", "TTGT", "TTTA", "TTTC", "TTTG", "TTTT"};

void print_to_file(pkentry_t* start, FILE* fp)
{
	static const char s[5] = "ACGT";
	while (1)
	{
		int i;
		for (i = 0; i < KMER_LENGTH / 4; i++)
		{
			fputs(lookup_table[start->data.data[i]], fp);
		}
		for (i *= 4; i < KMER_LENGTH; i++)
		{
			fputc(baseat(&start->data, i), fp);
		}
		if (start->next == NULL) return;

		start = start->next;
	}
}


// hashtable shenanigans

uint64_t hash(pkentry_t* pke)
{ // hashing order goes first char, back ext, rest of str.
	uint64_t h = pke->ext >> 4;
	if (KMER_PACKED_LENGTH < 8) {
		uint64_t t = 0;
		for (int i = 0; i < KMER_PACKED_LENGTH; i++)
		{
			t |= pke->data[i] << (i*8);
		}
		return (t << 3) | h;
	} else {
		uint64_t* d = (uint64_t*) &pke->data;
		return (*d << 3) | h;
	}
}

#endif
