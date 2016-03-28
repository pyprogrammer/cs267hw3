#ifndef UPC_HASH_H
#define UPC_HASH_H

typedef shared hash_table_t* hash_dir_t;
typedef shared memory_heap_t* mem_dir_t;

/* Creates hashtable for UPC */
/* Whenever something is called on a hash_table_t or memory_heap_t, instead call it on
 * hash_dir_t[MYTHREAD] and mem_dir_t[MYTHREAD]
 */
hash_dir_t *upc_create_hash_table(int64_t nEntries, mem_dir_t *memory_heap);

#endif
