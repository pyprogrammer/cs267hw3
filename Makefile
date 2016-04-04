CC = CC
UPCC = upcc

KMER_LENGTH 		= 19
KMER_PACKED_LENGTH 	= $(shell echo $$((($(KMER_LENGTH)+3)/4)))

# Add -std=gnu99 to CFLAGS if use gnu compiler
CFLAGS 	=  
DEFINE 	= -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH)
HEADERS	= contig_generation.h kmer_hash.h packingDNAseq.h
UPCHEADERS	= contig_upc.h kmer_hash_upc.h packingDNAseq.h
LIBS	=
UPCFLAGS = 

TARGETS	= serial pgen kh

all: 	$(TARGETS)

serial: serial.c $(HEADERS)
		$(CC) $(CFLAGS) -o $@ $< -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH) $(LIBS)

pgen:	pgen.upc $(UPCHEADERS)
		$(UPCC) $(UPCFLAGS) -Wc,"$(CFLAGS)" -o $@ $< $(DEFINE) $(LIBS)

kh: 	kmer_hash_upc.upc $(UPCHEADERS)
		$(UPCC) $(UPCFLAGS) -Wc,"$(CFLAGS)" -o $@ $< $(DEFINE) $(LIBS)

kh-test:
	make clean && make kh && sbatch job-kh

scratch:
	cp ~/code/hw3/* .
	
pgen-test:
	scancel -u ysg && make clean && make && sbatch job-p

s-pgen:
	make scratch && make pgen-test

s-kh:
	make scratch && make kh-test



clean :
	rm -f *.o
	rm -rf $(TARGETS)
	rm -rf *.err *.out
	scancel -u ysg
