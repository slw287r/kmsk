#include "kmsk.h"

int main(int argc, char *argv[])
{
	if (argc == 1)
		usage();
	time_t start_t, end_t;
	time(&start_t);
	int c = 0, k = KMER, suc_offset = 0, n = 0, v = 0, r = 0;
	unsigned i, j = 0, t = THREADS, w = WRAP;
	char *sbj = 0, *pfx = 0, *out = 0, *z = NULL;
	ketopt_t opt = KETOPT_INIT;
	while ((c = ketopt(&opt, argc, argv, 1, ":s:p:o:k:t:w:rhVv", 0)) >= 0)
	{
		if (c == 'h')
			usage();
		else if (c == 'v')
			return(puts(VERSION));
		else if (c == 's')
			sbj = opt.arg;
		else if (c == 'p')
			pfx = opt.arg;
		else if (c == 'k')
			k = atoi(opt.arg);
		else if (c == 'r')
			r = 1;
		else if (c == 't')
			t = atoi(opt.arg);
		else if (c == 'w')
			w = atoi(opt.arg);
		else if (c == 'V')
			v = 1;
		else if (c == '?')
			error("unknown option -%c\n", opt.opt ? opt.opt : ':');
		else if (c == ':')
			error("missing option argument: -%c\n", opt.opt ? opt.opt : ':');
	}
	if (!sbj || access(sbj, R_OK) == -1)
		error("Error accessing file [%s]\n", sbj);
	if (argc == opt.ind)
		error("Error at least one query fasta is required\n");
	if (k < 1 || k > KMER)
		error("Invalid kmer size: %d, must be in (1, %d]\n", k, KMER);
	suc_offset = asprintf(&z, "Loading kmers from [%s]", sbj);
	if (v)
		info("%s\n", z);
	km_t **h = calloc(NH, sizeof(km_t *));
	for (i = 0; i < NH; ++i)
		h[i] = km_init();
	rec_fa(sbj, k, v, h);
	/*
	puts("count of kmers");
	for (i = 0; i < NH; ++i)
		printf("%d\t%d\n", i, kh_size(h[i]));
	 * 0       714399117
	 * 1       290877234
	 * 2       260368586
	 * 3       139206549
	 * 4       562175559
	 * 5       319088177
	 * 6       139396056
	 * 7       86878792
	*/
	if (v)
		success(suc_offset);
	n = argc - opt.ind;
	suc_offset = asprintf(&z, "Masking sequences in [%d] files", n);
	if (v)
		info("%s\n", z);
	if (out && access(out, R_OK))
		mkdir_p(out);
	char cwd[PATH_MAX];
	getcwd(cwd, sizeof(cwd));
	// init threads
	threadpool thpool = thpool_init(min(n, t));
	// add msk work
	op_t *op = calloc(n, sizeof(op_t));
	for (i = opt.ind; i < argc; ++i)
	{
		char *qry = argv[i], p[NAME_MAX], e[NAME_MAX];
		if (access(qry, R_OK))
			error("Error accessing query fasta [%s]\n", qry);
		op[j].h = h;
		op[j].k = k;
		op[j].v = v;
		op[j].w = w;
		strcpy(op[j].i, qry);
		strcpy(p, basename(qry));
		if (ends_with(p, ".fa.gz"))
		{
			*(p + strlen(p) - strlen(".fa.gz")) = '\0';
			strcpy(e, "fa");
		}
		else if (ends_with(p, ".fna.gz"))
		{
			*(p + strlen(p) - strlen(".fna.gz")) = '\0';
			strcpy(e, "fna");
		}
		else if (ends_with(p, ".fasta.gz"))
		{
			*(p + strlen(p) - strlen(".fasta.gz")) = '\0';
			strcpy(e, "fasta");
		}
		else if (ends_with(p, ".fasta"))
		{
			*(p + strlen(p) - strlen(".fasta")) = '\0';
			strcpy(e, "fasta");
		}
		else if (ends_with(p, ".fna"))
		{
			*(p + strlen(p) - strlen(".fna")) = '\0';
			strcpy(e, "fna");
		}
		else if (ends_with(p, ".fa"))
		{
			*(p + strlen(p) - strlen(".fa")) = '\0';
			strcpy(e, "fa");
		}
		sprintf(op[j].o, "%s/%s.%s.%s", out ? out : cwd, p, pfx ? pfx : __progname, r ? "bed" : e);
		op[j].r = r;
		thpool_add_work(thpool, msk_fna, (void *)(uintptr_t)(op + j++));
	}
	thpool_wait(thpool);
	thpool_destroy(thpool);
	free(op);
	if (v)
	{
		success(suc_offset);
		time(&end_t);
		timer *run_time = runtime((int)difftime(end_t, start_t));
		if (run_time->d)
			info("runtime: %02dd:%02dh:%02dm:%02ds\n", run_time->d, run_time->h,
					run_time->m, run_time->s);
		else
			info("runtime: %02dh:%02dm:%02ds\n", run_time->h, run_time->m, run_time->s);
		free(run_time);
	}
	for (i = 0; i < NH; ++i)
		km_destroy(h[i]);
	free(h);
	free(z);
}

static void rec_ks(km_t **h, const int k, kseq_t *ks)
{
	int i, l, absent;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < ks->seq.l; ++i)
	{
		int c = seq_nt4_table[(uint8_t)ks->seq.s[i]];
		if (c < 4) // not an "N" base
		{
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) // we find a k-mer
			{
				uint64_t y = x[0] < x[1] ? x[0] : x[1];
				km_put(h[(uint8_t)(y & ~(~0U << 3))], y, &absent);
			}
		}
		else // if there is an "N", restart
			l = 0, x[0] = x[1] = 0;
	}
}

static void rec_fa(const char *fn, int k, const int v, km_t **h)
{
	gzFile fp;
	kseq_t *ks;
	if ((fp = gzopen(fn, "r")) == 0)
		error("Error reading file [%s]\n", fn);
	ks = kseq_init(fp);
	if (v)
	{
		while (kseq_read(ks) >= 0)
		{
			info("Processing %s %s ", ARW, ks->name.s);
			rec_ks(h, k, ks);
		}
	}
	else
		while (kseq_read(ks) >= 0)
			rec_ks(h, k, ks);
	kseq_destroy(ks);
	gzclose(fp);
}

static void dump_ks1(kseq_t *ks, const unsigned wr, FILE *fp)
{
	fputc(ks->is_fastq ? '@' : '>', fp);
	fputs(ks->name.s, fp);
	if (ks->comment.l)
	{
		fputc(' ', fp);
		fputs(ks->comment.s, fp);
	}
	fputc('\n', fp);
	if (!ks->is_fastq) // wrap line for fa
	{
		int i = 0;
		char *s = ks->seq.s;
		while (i < ks->seq.l)
		{
			fputc(*s++, fp);
			if (++i % wr == 0)
				fputc('\n', fp);
		}
		if (ks->seq.l % wr)
			fputc('\n', fp);
	}
	else
	{
		fputs(ks->seq.s, fp);
		fputs("\n+\n", fp);
		fputs(ks->qual.s, fp);
		fputc('\n', fp);
	}
}

static void msk_fna(void *_op)
{
	op_t *op = _op;
	km_t **h = op->h;
	const int k = op->k;
	const int v = op->v;
	const char *in = op->i;
	const char *out = op->o;
	const int r = op->r;
	const unsigned w = op->w;

	gzFile fp;
	size_t nsq = 0;
	if (!(fp = gzopen(in, "r")))
		error("Error reading %s\n", in);
	FILE *fo = out ? fopen(out, "w"): stdout;
	kseq_t *ks = kseq_init(fp);
	while (kseq_read(ks) >= 0)
	{
		++nsq;
		if (v)
			info("Masking %d contig%c", nsq, nsq > 1 ? 's' : ' ');
		size_t i, l, from, to;
		cgranges_t *cr = cr_init();
		uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
		for (i = l = 0, x[0] = x[1] = 0; i < ks->seq.l; ++i)
		{
			int c = seq_nt4_table[(uint8_t)ks->seq.s[i]];
			if (c < 4) // not an "N" base
			{
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= k)
				{
					uint64_t y = x[0] < x[1] ? x[0] : x[1];
					km_t *h1 = h[(uint8_t)(y & ~(~0U << 3))];
					if (km_get(h1, y) != kh_end(h1))
						cr_add(cr, ks->name.s, i - k + 1, i + 1, 0);
				}
			}
			else // if there is an "N", restart
				l = 0, x[0] = x[1] = 0;
		}
		cr_sort(cr);
		cr_merge_pre_index(cr);
		if (r)
		{
			for (i = 0; i < cr->n_r; ++i)
			{
				const cr_intv_t *r = &cr->r[i];
				from = (int32_t)r->x, to = (int32_t)r->y;
				fprintf(fo, "%s\t%zu\t%zu\n", ks->name.s, from, to);
			}
		}
		else
		{
			// masking
			char *s = ks->seq.s;
			for (i = 0; i < cr->n_r; ++i)
			{
				const cr_intv_t *r = &cr->r[i];
				size_t m = 0;
				from = (int32_t)r->x, to = (int32_t)r->y;
				for (m = from; m < to; ++m)
					s[m] = 'N';
			}
			dump_ks1(ks, w, fo);
		}
		cr_destroy(cr);
	}
	kseq_destroy(ks);
	gzclose(fp);
	fclose(fo);
}

static int is_dir(const char *path)
{
	struct stat statbuf;
	if (stat(path, &statbuf) != 0) return 0;
	return S_ISDIR(statbuf.st_mode);
}

static int mkdir_p(const char *path)
{
	// check for regular file
	if (!access(path, F_OK) && !is_dir(path))
		error("Not a directory [%s]\n", path);
	const size_t len = strlen(path);
	char _path[PATH_MAX];
	errno = 0;
	if (len > sizeof(_path)-1)
	{
		errno = ENAMETOOLONG;
		return EXIT_FAILURE;
	}
	strcpy(_path, path);
	char *p = 0;
	for (p = _path + 1; *p; p++)
	{
		if (*p == '/')
		{
			*p = '\0';
			if (mkdir(_path, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0)
				if (errno != EEXIST)
					return EXIT_FAILURE;
			*p = '/';
		}
	}
	if (mkdir(_path, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0)
		if (errno != EEXIST)
			return EXIT_FAILURE;
	return EXIT_SUCCESS;
}

static int ends_with(const char *str, const char *suffix)
{
	int ret = 0;
	int str_len = strlen(str);
	int suffix_len = strlen(suffix);
	if ((str_len >= suffix_len) && (0 == strcasecmp(str + (str_len-suffix_len), suffix)))
	ret = 1;
	return ret;
}

static void error(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	fputs(ERR, stderr);
	fputc(' ', stderr);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(EXIT_FAILURE);
}

static char *now(void)
{
	char buf[80];
	time_t _now = time(0);
	strftime(buf, sizeof(buf), "\e[34m%D %X\e[0m", localtime(&_now));
	return strdup(buf);
}

static timer *runtime (int diff_t)
{
	timer *my_timer = calloc(1, sizeof(timer));
	my_timer->d = diff_t / 86400;
	my_timer->h = diff_t % 86400 / 3600;
	my_timer->m = diff_t % 86400 % 3600 / 60;
	my_timer->s = diff_t % 86400 % 3600 % 60;
	return my_timer;
}

static void info(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	char *right_now = now();
	fprintf(stderr, "\r\e[K%s %s \e[33m", INF, right_now);
	vfprintf(stderr, format, ap);
	fputs("\e[0m", stderr);
	va_end(ap);
	free(right_now);
}

static void success(const int n)
{
	fprintf(stderr, "\r\e[A\r\e[%dC%s\n", n + 21, SUC);
}

static void horiz(const int _n)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	int i, n = (w.ws_col >= _n) ? _n : w.ws_col;
	for (i = 0; i < n; ++i)
		fputs(BAR, stdout);
	fputc('\n', stdout);
}

static void usage()
{
	char cwd[PATH_MAX];
    getcwd(cwd, sizeof(cwd));
	struct winsize win;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &win);
	int w = min(max(24 + strlen(cwd), 42), win.ws_col);
	horiz(w);
	puts("\e[1m Mask query fasta(s) by kmers in subject fasta\e[0m");
	horiz(w);
	fprintf(stdout, "Usage: \e[1;31m%s\e[0;0m \e[2m[options]\e[0m \e[33m<fna>\e[0m \e[2m[<fna>]\e[0m\n", __progname);
	putchar('\n');
	printf("  -s subject fasta file\n");
	printf("  -p output prefix \e[2m[%s]\e[0m\n", __progname);
	printf("  -o output directory \e[2m[%s]\e[0m\n", cwd);
	printf("  -r output region bed \e[2m[false]\e[0m\n");
	printf("  -k kmer size (max supported %d) \e[2m[%d]\e[0m\n", KMER, KMER);
	printf("  -t num of threads to mask fna \e[2m[%d]\e[0m\n", THREADS);
	printf("  -w wraping column for masked fna \e[2m[%d]\e[0m\n", WRAP);
	putchar('\n');
	puts("  -h display help message");
	puts("  -v display programme version");
	horiz(w);
	exit(1);
}
