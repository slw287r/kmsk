#pragma once
#define _GNU_SOURCE
#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <limits.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <zlib.h>

#include "cgranges.h"
#include "thpool.h"
#include "ketopt.h" // command-line argument parser
#include "kseq.h" // FASTA/Q parser

extern const char *__progname;

KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
KHASHL_MAP_INIT(, km_t, km, uint64_t, uint8_t, kh_hash_uint64, kh_eq_generic)
#define VERSION "0.1.3"
#define KMER 31
#define THREADS 8
#define WRAP 60
#define NH 8

#define min(x,y) ((x)>(y)?(y):(x))
#define max(x,y) ((x)>(y)?(x):(y))
#define basename(a) (strrchr(a, '/') ? strrchr(a, '/') + 1 : a)
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);
// symbols
#define ARW "\e[32m\xe2\x9e\x9c\e[0m"
#define BAR "\e[2m\xe2\x94\x80\e[0m"
#define BUL "\e[2m\xE2\x80\xA2\e[0m"
#define TAB "\e[2m\xe2\x87\xa5\e[0m"
#define ERR "\e[1;31m\xE2\x9C\x96\e[0;0m"
#define INF "\e[1;34m\xE2\x84\xb9\e[0;0m"
#define SUC "\e[1;31m\xE2\x9C\x94\e[0;0m"

typedef struct
{
	int d;
	int h;
	int m;
	int s;
} timer;

typedef struct
{
	km_t **h;
	int k, v, r;
	unsigned w;
	char i[PATH_MAX], o[PATH_MAX];
} op_t;

const uint8_t kmertochar[5] = { 'A', 'C', 'G', 'T', 'N' };
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

const int N_EXT = 6;
const char EXT[][16] = {
	".fa.gz",
	".fna.gz",
	".fasta.gz",
	".fasta",
	".fna",
	".fa"
};

static void usage();
static char *now(void);
static timer *runtime (int diff_t);
static void horiz(const int n);
static int is_dir(const char *path);
static int mkdir_p(const char *path);
static void error(const char *format, ...);
static void info(const char *format, ...);
static void success(const int n);
static int ends_with(const char *str, const char *suffix);
static void rec_ks(km_t **h, const int k, kseq_t *ks);
static void rec_fa(const char *fn, const int k, const int v, km_t **h);
static void dump_ks1(kseq_t *ks, const unsigned wr, FILE *fp);
static int mkdir_p(const char *path);
static void msk_fna(void *_op);

