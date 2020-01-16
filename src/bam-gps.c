#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <assert.h>
#include <libgen.h>
#include <ctype.h>
#include <getopt.h>
#include <sys/stat.h>
#include <zlib.h>

#include <bam.h>
#include <htslib/faidx.h>
#include <htslib/kfunc.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>

#define version "0.1.0"
#define MAX_DEP 8000
#define MIN_DEP 30
#define MIN_MAPQ 20
#define MIN_QLEN 10
#define header "CHROM\tPOS\tDEPTH\tREF\tA\tC\tG\tT\tN\tINS\tDEL"

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
	bamFile fp;
	bam_hdr_t *hdr;
	int min_mapq, min_len;
} aux_t;

typedef struct
{
	char *in, *out, *ref, *reg;
	int max_dep, min_dep;
} arg_t;

struct option long_options[] =
{
	{"in", required_argument, 0, 'i'},
	{"out", required_argument, 0, 'o'},
	{"ref", required_argument, 0, 'r'},
	{"reg", required_argument, 0, 'R'},
	{"max_dep", required_argument, 0, 'M'},
	{"min_dep", required_argument, 0, 'm'},
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
void pileup(aux_t *data, FILE *fp, faidx_t *fai, void *bed, const int max_dep, const int min_dep);

int main(int argc, char *argv[])
{
	int i = 0;
	if (argc == 1) usage(argv[0]);
	arg_t *arg = calloc(1, sizeof(arg_t));
	arg->max_dep = MAX_DEP;
	arg->min_dep = MIN_DEP;
	int c = 0, opt_idx = 0;
	uint64_t required_args = 0;
	while ((c = getopt_long(argc, argv,"i:o:r:R:M:m:vh", long_options, &opt_idx)) != -1)
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
	FILE *fp = fopen(arg->out, "w");
	if (!fp) error("Error creating output file: %s\n", arg->out);
	fputs(header, fp);
	fputc('\n', fp);
	faidx_t *fai = fai_load(arg->ref); assert(fai);
	// prepare to read bam
	aux_t *data = calloc(1, sizeof(aux_t)); assert(data);
	data->fp = bam_open(arg->in, "r"); assert(data->fp);
	data->hdr = bam_hdr_read(data->fp); assert(data->hdr);
	data->min_len = MIN_QLEN;
	data->min_mapq = MIN_MAPQ;
	bam_hdr_t *h = data->hdr;
	pileup(data, fp, fai, bed, arg->max_dep, arg->min_dep);
	// infer depth
	bam_hdr_destroy(data->hdr);
	bam_close(data->fp);
	bed_destroy(bed);
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
		ret = bam_read1(aux->fp, b);
		if (ret < 0)
			break;
		if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))
			continue;
		if ((int)b->core.qual < aux->min_mapq)
			continue;
		if (bam_cigar2qlen(&b->core, bam_get_cigar(b)) < aux->min_len)
			continue;
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

void pileup(aux_t *data, FILE *fp, faidx_t *fai, void *bed, const int max_dep, const int min_dep)
{
	khint_t k;
	int i, l = 0, tid = 0, pos = 0, n_plp = 0, r = 0;
	khash_t(map) *mi = kh_init(map);
	khash_t(map) *md = kh_init(map);
	bam_hdr_t *h = data->hdr;
	const bam_pileup1_t *plp = calloc(1, sizeof(bam_pileup1_t)); assert(plp);
	bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)&data);
	bam_mplp_set_maxcnt(mplp, max_dep);
	kstring_t iseq = {0, 0, 0}, dseq = {0, 0, 0};
	while (bam_mplp_auto(mplp, &tid, &pos, &n_plp, &plp) > 0) // pos: 0-based
	{
		if (n_plp < min_dep) continue;
		if (!bed_overlap(bed, h->target_name[tid], pos, pos + 1))
			continue;
		// reference sequence
		unsigned nt[5] = {0}, idl[2] = {0}, gap = 0;
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
		// output info
		char *ref = faidx_fetch_seq(fai, h->target_name[tid], pos, pos, &l);
		fprintf(fp, "%s\t", h->target_name[tid]);
		fprintf(fp, "%d\t", pos + 1);
		fprintf(fp, "%d\t", n_plp - gap);
		fprintf(fp, "%s\t", ref);
		fprintf(fp, "%d\t", nt[0]);
		fprintf(fp, "%d\t", nt[1]);
		fprintf(fp, "%d\t", nt[2]);
		fprintf(fp, "%d\t", nt[3]);
		fprintf(fp, "%d\t", nt[4]);
		fprintf(fp, "%s\t", iseq.s);
		fprintf(fp, "%s\n", dseq.s);
		EMPTY_HASH(map, mi);
		EMPTY_HASH(map, md);
		iseq.l = dseq.l = 0;
		free(ref);
	}
	bam_mplp_destroy(mplp);
}

void usage(char *str)
{
	fprintf(stderr, "\n\033[1;2mGet pileup summary for GRAIL's conta\033[0;0m\n");
	fprintf(stderr, "\n\033[1;2mBuild %s on %s at %s\033[0;0m\n", version, __DATE__, __TIME__);
	fprintf(stderr, "\n\033[1mUsage\033[0m: \033[1;31m%s\033[0;0m -i <bam> -o <out> -r <ref> -R <reg>\n", basename(str));
	fprintf(stderr, "\n  Required:\n");
	fprintf(stderr, "    -i --in         [STR] Input bam\n");
	fprintf(stderr, "    -o --out        [STR] Output pileup summary\n");
	fprintf(stderr, "    -r --ref        [STR] Reference fa\n");
	fprintf(stderr, "    -R --reg        [STR] Region of interest\n");
	fprintf(stderr, "\n  Optional:\n");
	fprintf(stderr, "    -M --max_dep    [INT] Max depth of coverage (8000)\n");
	fprintf(stderr, "    -m --min_dep    [INT] Min depth of coverage (10)\n");
	fputc('\n', stderr);
	fprintf(stderr, "    -v --version          Show version\n");
	fprintf(stderr, "    -h --help             Show help message\n");
	fprintf(stderr, "\n\033[1mContact\033[0m: liujc@geneplus.org.cn\n\n");
	exit(EXIT_FAILURE);
}
