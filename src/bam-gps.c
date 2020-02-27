#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <assert.h>
#include <libgen.h>
#include <ctype.h>
#include <stdint.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/stat.h>
#include <zlib.h>

#include <bam.h>
#include <bedidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/kfunc.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>
#include <htslib/thread_pool.h>
#include <pthread.h>

#define version "0.2.0"
#define MAX_DEP 10000
#define MIN_DEP 10
#define MIN_MAPQ 50
#define MIN_QLEN 30
#define THREAD 8
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);
#define header "CHROM\tPOS\tDEPTH\tREF\tA\tC\tG\tT\tN\tINS\tDEL"
static const unsigned int BAM_FILTER = (BAM_FQCFAIL | BAM_FUNMAP | BAM_FMUNMAP | BAM_FSECONDARY | BAM_FDUP);
pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
	int n, m;
	uint64_t *a;
	int *idx;
	int filter;
} bed_reglist_t;

KHASH_MAP_INIT_STR(reg, bed_reglist_t)
typedef khash_t(reg) reghash_t;
//KHASH_SET_INIT_STR(set)

#define EMPTY_HASH(type, hash) do                 \
{                                                 \
	khint_t k;                                    \
	khash_t(type) *h_map = (khash_t(type) *)hash; \
	for (k = 0; k < kh_end(h_map); ++k)           \
	{                                             \
		if (kh_exist(h_map, k))                   \
		{                                         \
			free((char *)kh_key(h_map, k));       \
			kh_del(type, h_map, k);               \
		}                                         \
	}                                             \
} while (0)

KHASH_MAP_INIT_STR(map, int)

typedef struct
{
	samFile *fp;
	bam_hdr_t *hdr;
	hts_itr_t *itr;
	hts_idx_t *idx;
	int min_mapq, min_qlen;
} aux_t;

typedef struct
{
	char *in, *out, *ref, *reg;
	int max_dep, min_dep, biallelic;
	int thread;
} arg_t;

typedef struct
{
	uint64_t id;
	char *str;
} buf_t;

typedef struct
{
	char *reg;
	arg_t *arg;
	faidx_t *fai;
	buf_t *buf;
} par_t;

struct option long_options[] =
{
	{"in", required_argument, 0, 'i'},
	{"out", required_argument, 0, 'o'},
	{"ref", required_argument, 0, 'r'},
	{"reg", required_argument, 0, 'R'},
	{"max_dep", required_argument, 0, 'M'},
	{"min_dep", required_argument, 0, 'm'},
	{"biallelic", no_argument, 0, 'b'},
	{"thread", no_argument, 0, 't'},
	{"version", no_argument, 0, 'v'},
	{"help", no_argument, 0, 'h'},
	{0, 0, 0, 0}
};

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void error(const char *format, ...);
// make sure required args are filled
void chk_arg(uint64_t ra, int na, char *s);
// check for invalid args
void val_arg(const arg_t *);
void _mkdir(const char *dir);
void usage(char *s);
int read_bam(void *data, bam1_t *b);
int buf_cmp (void const *a, void const *b);
void *gps(void *par);

int main(int argc, char *argv[])
{
	int i = 0;
	if (argc == 1) usage(argv[0]);
	arg_t *arg = calloc(1, sizeof(arg_t));
	arg->max_dep = MAX_DEP;
	arg->min_dep = MIN_DEP;
	arg->thread = THREAD;
	int c = 0, opt_idx = 0;
	uint64_t required_args = 0;
	while ((c = getopt_long(argc, argv,"i:o:r:R:M:m:t:bvh", long_options, &opt_idx)) != -1)
	{
		switch (c)
		{
			case 'i':
				arg->in = optarg;
				required_args |= 1<<0;
				break;
			case 'o':
				arg->out = optarg;
				required_args |= 1<<1;
				break;
			case 'r':
				arg->ref = optarg;
				required_args |= 1<<2;
				break;
			case 'R':
				arg->reg = optarg;
				required_args |= 1<<3;
				break;
			case 'M':
				arg->max_dep = atoi(optarg);
				break;
			case 'm':
				arg->min_dep = atoi(optarg);
				break;
			case 't':
				arg->thread = atoi(optarg);
				break;
			case 'b':
				arg->biallelic = 1;
				break;
			case 'v':
				fputs(version, stderr);
				fputc('\n', stderr);
				return(EXIT_SUCCESS);
			case 'h':
				usage(argv[0]);
				return EXIT_SUCCESS;
			case ':':
				fprintf(stderr, "Option -%c requires an argument\n", optopt);
				return EXIT_FAILURE;
				break;
			case '?':
				fprintf(stderr, "Option -%c is undefined\n", optopt);
				return EXIT_FAILURE;
				break;
		}
	}
	chk_arg(required_args, 4, argv[0]);
	val_arg(arg);
	char *output = strdup(arg->out);
	_mkdir(dirname(output));
	//free(output);
	void *bed = bed_read(arg->reg);
	if (!bed) error("Failed to read region file %s\n", arg->reg);
	bed_unify(bed);
	faidx_t *fai = fai_load(arg->ref); assert(fai);

	hts_tpool *proc = hts_tpool_init(arg->thread);
	hts_tpool_process *que = hts_tpool_process_init(proc, arg->thread * 2, 1);
	// dispatch jobs
	bed_reglist_t *p;
	reghash_t *h = (reghash_t *)bed;
	uint64_t m = 0, n = 0;
	for (khint_t k = kh_begin(h); k != kh_end(h); ++k)
		if (kh_exist(h,k) && (p = &kh_val(h,k)) != NULL && p->n > 0)
			m += p->n;
	buf_t *buf = calloc(m, sizeof(buf_t)); assert(buf);
	for (khint_t k = kh_begin(h); k != kh_end(h); ++k)
	{
		if (kh_exist(h,k) && (p = &kh_val(h,k)) != NULL && p->n > 0)
		{
			const char *chr = kh_key(h, k);
			for (i = 0; i < p->n; ++i)
			{
				par_t *par = calloc(1, sizeof(par_t));
				par->arg = arg;
				par->fai = fai;
				par->buf = buf + n++;
				asprintf(&par->reg, "%s:%d-%d", chr,  (uint32_t)(p->a[i]>>32) + 1, (uint32_t)(p->a[i]));
				hts_tpool_dispatch(proc, que, gps, par);
			}
		}
	}
	hts_tpool_process_flush(que);
	hts_tpool_process_destroy(que);
	hts_tpool_destroy(proc);
	fai_destroy(fai);
	bed_destroy(bed);

	BGZF *in = bgzf_open(arg->in, "r"); assert(in);
	bam_hdr_t *hdr = bam_hdr_read(in);
	qsort(buf, n, sizeof(buf_t), buf_cmp);
	FILE *fp = fopen(arg->out, "w");
	if (!fp) error("Error creating output file: %s\n", arg->out);
	fputs(header, fp);
	fputc('\n', fp);
	for (i = 0; i < n; ++i)
	{
		if ((buf + i)->id)
		{
			fprintf(fp, "%s\t%d\t%s", hdr->target_name[(buf + i)->id >> 32], (uint32_t)(buf + i)->id + 1, (buf + i)->str);
			free((buf + i)->str);
		}
	}
	bgzf_close(in);
	bam_hdr_destroy(hdr);
	free(buf);
	fclose(fp);
	free(arg);
	return 0;
}

void error(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(EXIT_FAILURE);
}


void chk_arg(uint64_t ra, int na, char *s)
{
	int required_args_set = __builtin_popcount(ra & ((1<<na) - 1));
	if (required_args_set != na)
	{
		fprintf(stderr, "\033[31m[ERROR]\033[0m %d required argument%s unspecified\n", na - required_args_set, (na - required_args_set > 1) ? "s" : "");
		usage(s);
		exit(1);
	}
}

static void chk_file(const char *fn)
{
	if (fn && access(fn, F_OK) == -1)
	{
		fprintf(stderr, "Error: Failed to read file: %s\n", fn);
		exit(1);
	}
}

void val_arg(const arg_t *arg)
{
	chk_file(arg->in);
	chk_file(arg->reg);
	chk_file(arg->ref);
	if (arg->max_dep <= 0 || arg->min_dep <= 0)
		error("Expecting positive depth value\n");
	if (arg->max_dep < arg->min_dep)
		error("Expecting larget max_dep (%d) than min_dep (%d)\n", arg->max_dep, arg->min_dep);
}

size_t count_chars(char* s, char c)
{
	return *s ? ((c == *s) + count_chars(s + 1, c)) : 0;
}

void _mkdir(const char *dir)
{
	char tmp[PATH_MAX];
	char *p = NULL;
	size_t len;
	snprintf(tmp, sizeof(tmp),"%s",dir);
	len = strlen(tmp);
	if(tmp[len - 1] == '/') tmp[len - 1] = 0;
	for(p = tmp + 1; *p; p++)
	{
		if(*p == '/')
		{
			*p = 0;
			mkdir(tmp, (S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH));
			*p = '/';
		}
	}
	mkdir(tmp, (S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH));
}

int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data;
	int ret;
	while (1)
	{
		ret = sam_itr_next(aux->fp, aux->itr, b);
		if (ret < 0) break;
		if (b->core.flag & BAM_FILTER) continue;
		break;
	}
	return ret;
}

/**
 * @brief get insertion or deletion sequence
 *
 * @param p pointer to bam_pileup1_t
 *
 * @return insertion or deletion sequence
 */
char *get_id_seq(const bam_pileup1_t *p)
{
	int i = 0, l = abs(p->indel);
	char *seq = calloc(l, sizeof(char)); assert(seq);
	for (i = 1; i <= l; ++i) // w/ leading anchor base
		seq[i - 1] = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos + i)];
	seq[i - 1] = '\0';
	return seq;
}

int buf_cmp (void const *a, void const *b)
{
	return (*(buf_t *)b).id > (*(buf_t *)a).id ? -1 : 1;
}

void *gps(void *par)
{
	khint_t k;
	int tid = 0, beg, end, pos = 0;
	int i, l = 0, n_plp = 0, r = 0;
	kstring_t iseq = {0, 0, 0}, dseq = {0, 0, 0};
	khash_t(map) *mi = kh_init(map), *md = kh_init(map);

	par_t *p = (par_t *)par;
	aux_t *data = calloc(1, sizeof(aux_t)); assert(data);
	data->fp = hts_open(p->arg->in, "r"); assert(data->fp);
	data->hdr = sam_hdr_read(data->fp); assert(data->hdr);
	data->idx = sam_index_load(data->fp, p->arg->in); assert(data->idx);
	data->itr = sam_itr_querys(data->idx, data->hdr, p->reg); assert(data->itr);
	hts_idx_destroy(data->idx);
	beg = data->itr->beg;
	end = data->itr->end;
	const bam_pileup1_t *plp = calloc(1, sizeof(bam_pileup1_t)); assert(plp);
	bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)&data);
	bam_mplp_set_maxcnt(mplp, p->arg->max_dep);
	while (bam_mplp_auto(mplp, &tid, &pos, &n_plp, &plp) > 0)
	{
		if (pos < beg) continue;
		if (pos >= end) break;
		if (tid >= data->hdr->n_targets) continue;
		if (n_plp < p->arg->min_dep) continue;
		// reference sequence
		int nt[5] = {0}, idl[2] = {0}, gap = 0;
		for (i = 0; i < n_plp; ++i)
		{
			char *seq = 0;
			const bam_pileup1_t *p = plp + i;
			if (p->indel) seq = get_id_seq(p);
			if (p->indel > 0)
			{
				if ((k = kh_get(map, mi, seq)) == kh_end(mi))
				{
					k = kh_put(map, mi, strdup(seq), &r);
					kh_val(mi, k) = 1;
				}
				else
					++kh_val(mi, k);
			}
			if (p->indel < 0)
			{
				if ((k = kh_get(map, md, seq)) == kh_end(md))
				{
					k = kh_put(map, md, strdup(seq), &r);
					kh_val(md, k) = 1;
				}
				else
					++kh_val(md, k);
			}
			if (p->is_del || p->is_refskip)
			{
				++gap;
				continue;
			}
			++nt[seq_nt16_int[seq_nt16_table[(int)bam_nt16_rev_table[bam1_seqi((char *)bam1_seq(p->b), p->qpos)]]]];
			if (seq) free (seq);
		}
		if (n_plp - gap < p->arg->min_dep) continue;
		// summarize indels
		if (kh_size(mi))
		{
			for (k = kh_begin(mi); k != kh_end(mi); ++k)
				if (kh_exist(mi, k))
					ksprintf(&iseq, "%s:%d,", kh_key(mi, k), kh_val(mi, k));
			iseq.s[--iseq.l] = '\0';
		}
		else
			kputs("0", &iseq);
		if (kh_size(md))
		{
			for (k = kh_begin(md); k != kh_end(md); ++k)
				if (kh_exist(md, k))
					ksprintf(&dseq, "%s:%d,", kh_key(md, k), kh_val(md, k));
			dseq.s[--dseq.l] = '\0';
		}
		else
			kputs("0", &dseq);
		if (p->arg->biallelic)
		{
			if ((nt[0]>0) + (nt[1]>0) + (nt[2]>0) + (nt[3]>0) + count_chars(iseq.s, ':') + count_chars(dseq.s, ':') > 2)
				goto cleanup;
		}
		// output info
		/*
		 * https://github.com/samtools/htslib/issues/404#issuecomment-251434613
		 */
		pthread_mutex_lock(&mtx);
		char *ref = faidx_fetch_seq(p->fai, data->hdr->target_name[tid], pos, pos, &l);
		pthread_mutex_unlock(&mtx);
		p->buf->id = (uint64_t)tid << 32 | pos;
		asprintf(&p->buf->str, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", n_plp - gap, ref,
				nt[0], nt[1], nt[2], nt[3], nt[4], iseq.s, dseq.s);
		free(ref);
cleanup:
		EMPTY_HASH(map, mi);
		EMPTY_HASH(map, md);
		iseq.l = dseq.l = 0;
	}
	bam_mplp_destroy(mplp);
	bam_hdr_destroy(data->hdr);
	hts_itr_destroy(data->itr);
	hts_close(data->fp);
	free(data);
	free(p->reg);
	free(p);
	return NULL;
}

void usage(char *str)
{
	fprintf(stderr, "\n\033[1;2mGet pileup summary for GRAIL's conta\033[0;0m\n");
	fprintf(stderr, "\n\033[1;2mBuild %s on %s at %s\033[0;0m\n", version, __DATE__, __TIME__);
	fprintf(stderr, "\n\033[1mUsage\033[0m: \033[1;31m%s\033[0;0m -i <bam> -o <out> -r <ref> -R <reg>\n", basename(str));
	fprintf(stderr, "\n  Required:\n");
	fprintf(stderr, "    -i --in         [STR]  Input bam\n");
	fprintf(stderr, "    -o --out        [STR]  Output pileup summary\n");
	fprintf(stderr, "    -r --ref        [STR]  Reference fa\n");
	fprintf(stderr, "    -R --reg        [STR]  Regions of interest\n");
	fprintf(stderr, "\n  Optional:\n");
	fprintf(stderr, "    -M --max_dep    [INT]  Max depth of coverage (10000)\n");
	fprintf(stderr, "    -m --min_dep    [INT]  Min depth of coverage (10)\n");
	fprintf(stderr, "    -b --biallelic  [BOOL] Restrict alleles to biallelic (false)\n");
	fprintf(stderr, "    -t --thread     [INT]  Number of threads (8)\n");
	fputc('\n', stderr);
	fprintf(stderr, "    -v --version          Show version\n");
	fprintf(stderr, "    -h --help             Show help message\n");
	fprintf(stderr, "\n\033[1mContact\033[0m: liujc@geneplus.org.cn\n\n");
	exit(EXIT_FAILURE);
}
