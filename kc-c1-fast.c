
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fcntl.h>
#include <limits.h>
#include <string.h>

typedef struct __ReadSeqList {
	char* sequence;
	unsigned length;
	struct __ReadSeqList* next;
} ReadSeqList;

typedef struct HashTable {
	unsigned int bits;
	unsigned int count;
	unsigned int read_count;
	//unsigned int *collition;
	unsigned long long int *keys;
    unsigned int *values;
} HashTable;



const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

HashTable* HashTable_init(unsigned int bits, unsigned int read_count){
	unsigned int capacity = 1U << bits;
    HashTable *ht;
    ht = (HashTable*)calloc(1, sizeof(HashTable));
    ht->keys = (unsigned long long int*)calloc(capacity, sizeof(unsigned long long int));
    ht->values = (unsigned int*)calloc(capacity, sizeof(unsigned int));

	ht->bits = bits;
	ht->count = 0;

    return ht;
}


void HashTable_destory(HashTable *ht) {
	if (!ht) return;
	free((void *)ht->keys);
	free((void *)ht->values);
	free(ht);
}


// funcion para calcular un hash de 64 bits
unsigned int hash_uint64(unsigned long long int key) {

	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return (unsigned int)key;
}

static inline unsigned int h2b(unsigned int hash, unsigned int bits) {
    return hash * 2654435769U >> (32 - bits);
}

void hash_insert(HashTable *ht, unsigned long long int kmer) {

	unsigned int iKey, last;

    iKey = last = h2b(hash_uint64(kmer), ht->bits);

    while (ht->values[iKey] > 0 && ht->keys[iKey] != kmer) {
        iKey = (iKey + 1U) & ((1U << ht->bits) - 1);
        if (iKey == last) break;
    }
    // Comprobar si se ha encontrado un slot vacío
    if (ht->values[iKey] == 0) { // no se ha encontrado la llave

        ht->keys[iKey] = kmer;
        ht->values[iKey] = 1;
        ++ht->count;

    } else {
        ht->values[iKey]++;
    } 
}

// insert k-mers in $seq to hash table $ht
void count_seq_kmers(HashTable *ht, int k, int len, char *seq)
{

    int i, l;
    unsigned long long int x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;

    for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        int c = seq_nt4_table[(unsigned char)seq[i]];
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (unsigned long long int)(3 - c) << shift;  // reverse strand
            if (++l >= k) { // we find a k-mer

                unsigned long long int kmer = x[0] < x[1]? x[0] : x[1];
                hash_insert(ht, kmer); // only add one strand!

            }
        } else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
    }
}

void kernel_print_hist(const HashTable *ht)
{
    uint32_t k;
    uint64_t cnt[256];
    int i;
    for (i = 0; i < 256; ++i) cnt[i] = 0;
    for (k = 0; k < (1U<<(ht)->bits); ++k)
        if (ht->values[k] != 0)
            ++cnt[ht->values[k] < 256 ? ht->values[k] : 255];
    for (i = 1; i < 256; ++i)
        printf("%d\t%ld\n", i, (long)cnt[i]);	
}

static int count_file(const char *fn, int k, unsigned int p)
{
	HashTable *ht;
    unsigned int read_count = 0;


    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(fn, "r");
    if (fp == NULL) exit(EXIT_FAILURE);

   
	ReadSeqList *current, *head;
	head = current = NULL;

    while ((read = getline(&line, &len, fp)) != -1) {

        read_count++;

		ReadSeqList *node = (ReadSeqList*)malloc(sizeof(ReadSeqList));
        node->sequence = (char*)malloc(strlen(line));
        strcpy(node->sequence, line);
        node->length = read;
        node->next =NULL;

        if(head == NULL){
            current = head = node;
        } else {
            current = current->next = node;
        }


    }

    fclose(fp);
    if (line) free(line);


    // inicializar hashtable
	ht = HashTable_init(p, read_count);



    unsigned int i;
    char* reads = malloc(head->length * read_count * sizeof(char));

	for(i=0, current = head; current; current=current->next){
        memcpy(reads + (i * head->length), current->sequence, head->length);
		i++;
    }


    printf("total reads: %d\n", read_count);
    printf("COUNT: %d\n\n", ht->count);
    
    clock_t begin = clock();

    current = head;
    for(i = 0; i<read_count ; i++){
        count_seq_kmers(ht, k, head->length, reads + (i * head->length));
    }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("CPU time: %f\n", time_spent);

	kernel_print_hist(ht);


	// limpieza
	for(current = head; current; current=current->next){
        free(current->sequence);
        free(current);
    }


    free(reads);
	HashTable_destory(ht);

    return 0;
}


int main(int argc, char *argv[])
{
	int k = 31;
    unsigned int p = 27;

    k = (int)strtol(argv[1], NULL, 10);
    p = (unsigned int)strtol(argv[2], NULL, 10);
	count_file(argv[3], k, p);

	return 0;
}
