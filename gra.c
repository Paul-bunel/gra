/**
 * @file gra.c
 * @author Paul Bunel (paul.bunel@etu.umontpellier.fr)
 * @brief Modifie un fichier multifasta pour retirer les gènes n'étant pas des
 * gènes du goût, et ajoute dans le header des gènes restant la famille de gène
 * correspondante ainsi que son état (fonctionnel/pseudogène) 
 * @param file Le fichier fasta d'entrée contenant les gènes à tester
 * @return Un nouveau fichier fasta contenant les gènes validés annotés
 * @version 0.1
 * @date 2022-05-13
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

struct String {
    size_t len;
    char *str;
};

typedef struct Tokens {
    size_t size;
    char **tokens;
} Tokens;

typedef struct MapGeneProfile {
    char *gene;
    char *profile;
    char *e_value;
} MapGeneProfile;

Tokens split_string(const char *str, char *sep) {
    size_t count_tokens = 0;
    int c = 0;
 	for (int i = 0; i < strlen(str) - 1; i++) {
        if (!c && str[i] >= 33 && str[i] < 127) { c = 1; }
  		if (c && str[i] == sep[0] && str[i + 1] != sep[0] && str[i + 1] != '\0') {
  			count_tokens++;
 		}
	}
    char **tokens = malloc(++count_tokens * sizeof(char*));

    char *s = strdup(str);
    char *tok = s, *end = s;
    int i = 0;
    while (tok != NULL) {
        strsep(&end, sep);
        if (tok[0] != '\0') {
            tokens[i] = strdup(tok);
            i++;
        }
        tok = end;
    }

    Tokens t = (Tokens) { .size = count_tokens, .tokens = tokens };

    return t;
}

struct String getFiles(int n, char **filenames) {
    FILE **files = calloc(n, sizeof(FILE*));

    for (int i = 0; i < n; i++) {
        files[i] = fopen(filenames[i], "r");
        if (files[i] == NULL) {
            fprintf(stderr, "Impossible d'ouvrir le fichier %s\n", filenames[i]);
            exit(EXIT_FAILURE);
        }
    }

    int total_fsize = 0;
    int *sizes = calloc(n, sizeof(int));
    for (int i = 0; i < n; i++) {
        fseek(files[i], 0, SEEK_END);
        total_fsize += ftell(files[i]);
        sizes[i] = ftell(files[i]);
        fseek(files[i], 0, SEEK_SET);
    }

    printf("Taille du premier fichier : %d, taille totale : %d\n", sizes[0], total_fsize);

    char *fasta = calloc(total_fsize + 1, sizeof(char));
    char *buffer = calloc(total_fsize + 1, sizeof(char));

    for (int i = 0; i < n; i++) {
        fread(buffer, sizes[i], 1, files[i]);
        strcat(fasta, buffer);
        fclose(files[i]);
    }

    fasta[total_fsize] = '\0';

    struct String res = {total_fsize + 1, fasta};

    return res;
}

char complement(char c) {
    switch (c) {
        case 'A':
            return 'T';
            break;
        case 'a':
            return 't';
            break;
        case 'T':
            return 'A';
            break;
        case 't':
            return 'a';
            break;
        case 'G':
            return 'C';
            break;
        case 'g':
            return 'c';
            break;
        case 'C':
            return 'G';
            break;
        case 'c':
            return 'g';
            break;
        default:
            return '\0';
            break;
    }
}

char* reverse_complement(char* seq, int crlf) {
    size_t s = strlen(seq), j = 0, offset = 0;
    char* rev = calloc(s + s/60 + 2, sizeof(char));
    strcpy(rev, seq);
    char c;
    for (int i = s; i >= 0; i--) {
        c = complement(seq[i]);
        if (c != '\0') { rev[j++] = c; offset++; }
        if (offset == 60) {
            if (crlf) { rev[j++] = '\r'; rev[j++] = '\n'; }
            else { rev[j++] = '\n'; }
            offset = 0;
        }
    }
    if (rev[j-1] != '\n') {
        if (crlf) { rev[j++] = '\r'; rev[j++] = '\n'; }
        else { rev[j++] = '\n'; }
    }
    rev[j] = '\0';

    return rev;
}

const size_t count_gene(const char *filename) {
    FILE *fd = fopen(filename, "r");
    char buffer[256];
    size_t n_gene = 0;
    while (fgets(buffer, sizeof(buffer), fd) != NULL) {
        if (buffer[0] == '>') { n_gene++; }
    }

    fclose(fd);
    return n_gene;
}

FILE* hmmscan(const char *fasta, const char *profile, const char *e_value,
                                                           const char *cpu) {
    /**
        hmmscan -o "/dev/null" \
        --tblout "RESULTS/${name^^}/res_$DBname.txt" \
        --qformat $format \
        "HMM_PROFILE/${name^^}/$name.hmm" \
        $db
     */
    char *fun = "/bin/bash -c '"
        "hmmscan -o /dev/null "
        "--tblout /dev/stdout "
        "--qformat fasta "
        "-E ";

    size_t l_cmd = strlen(fun) + strlen(e_value) + strlen(cpu) +
        strlen(profile) + strlen(fasta) + 11;
    char *cmd = calloc(l_cmd, sizeof(char));
    snprintf(cmd, l_cmd, "%s%s --cpu %s %s %s'", fun, e_value, cpu, profile, fasta);
    cmd[l_cmd - 1] = '\0';

    printf("Command :\n%s\n", cmd);
    
    FILE *scan_res = popen(cmd, "r");
    if (scan_res == NULL) {
        fprintf(stderr, "Erreur lors du scan HMM");
    }

    free(cmd);
    return scan_res;
}

/**
 * @brief Fonction de parsage du fichier table résultat de hmm_scan. La fonction
 * va récupérer chaque référence de gène présent dans le fichier résultat ainsi
 * que le profil associé
 * @param scan_res 
 * @return Dictionnaire associant une référence de gène et son profil attribué
 */
size_t parse_scan_res(
    FILE *scan_res,
    const size_t n_gene,
    MapGeneProfile *retrievedGenes
) {
    char buffer[512];
    size_t cpt = 0;

    while (fgets(buffer, sizeof(buffer), scan_res) != NULL) {
        if (buffer[0] != '#') {
            Tokens tokens = split_string(buffer, " ");
            int retrieved = 0;

            for (int i = 0; i < cpt; i++) {
                if (strcmp(retrievedGenes[i].gene, tokens.tokens[2]) == 0) {
                    retrieved = 1;
                }
            }

            if (!retrieved) {
                retrievedGenes[cpt++] = (MapGeneProfile) {
                    .gene = tokens.tokens[2],
                    .profile = tokens.tokens[0],
                    .e_value = tokens.tokens[4]
                };
            }

            free(tokens.tokens);
        }
    }
    fclose(scan_res);

    return cpt;
}

char* get_filename(const char* filename) {
    Tokens tokens = split_string(filename, "/");
    int i;
    Tokens raw_filename = split_string(tokens.tokens[tokens.size - 1], ".");

    free(tokens.tokens);
    return raw_filename.tokens[0];
}

char* default_output_name(const char* fasta) {
    char *fasta_truncated = get_filename(fasta);
    size_t l_fa_out_name = strlen(fasta_truncated) + strlen("_GRA_OUTPUT.fasta") + 1;
    char *fa_out_name = calloc(l_fa_out_name, sizeof(char));
    snprintf(fa_out_name, l_fa_out_name, "%s_GRA_OUTPUT.fasta", fasta_truncated);

    return fa_out_name;
}

/**
 * @brief fonction écrivant dans un nouveau fichier fasta les gènes reconnus
 * par les profils HMM, annoté avec le nom du profil l'ayant reconnu.
 * @param fasta le fichier fasta d'entrée contenant tous les gènes à tester
 * @param retrievedGenes les gènes reconnus
 * @param n_gene le nombre de gènes dans le fichier fasta d'entrée
 */
void fasta_output(const char* fasta, MapGeneProfile *retrievedGenes,
                                    size_t n_retrieved_genes, char* fa_out_name) {
    FILE *fa_in = fopen(fasta, "r");

    FILE *fa_out = fopen(fa_out_name, "w");

    char seq[32768]; seq[0] = '\0';
    char buffer[256];
    int read = 0;
    int lol = 0;
    while (fgets(buffer, sizeof(buffer), fa_in) != NULL) {
        if (buffer[0] == '>') {
            fprintf(fa_out, "%s", seq);
            memset(seq, 0, sizeof(seq));
            read = 0;
            for (int i = 0; i < n_retrieved_genes; i++) {
                if (strstr(buffer, retrievedGenes[i].gene) != NULL) {
                    lol++;
                    read = 1;
                    
                    int j = 0;
                    while (buffer[j] != '\n') {
                        j++;
                    }
                    if (buffer[j-1] == '\r') { j--; }
                    buffer[j] = '|';
                    buffer[j+1] = '\0';
                    strcat(buffer, retrievedGenes[i].profile);
                    strcat(buffer, "|\0");
                    strcat(buffer, retrievedGenes[i].e_value);
                    strcat(buffer, "\n\0");

                    strcat(seq, buffer);
                }
            }
        } else if (read) { strcat(seq, buffer); }
    }
    printf("Nombre de gènes écrits : %i\n", lol);
    fprintf(fa_out, "%s", seq);
    fclose(fa_in);
    fclose(fa_out);
}

const char** get_parameters(int argc, char ** argv, size_t *j) {
    /**
     * @brief Parse les arguments entrés en ligne de commande afin de récupérer
     * les différents paramètres du programme
     * ./gra <seqfile>
     * -h help
     * -p hmm profile
     * -e e_value threshold
     * -o output file
     * -c number of thread 
     */
    const char *hmm = "HMM_PROFILE/TAS/TAS_ncbi.hmm";
    const char *e_value = "0.9e-30";
    const char *fa_out = default_output_name(argv[argc-1]);
    const char *cpu = "2";
    const char *rc = "n";
    *j = 1;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            // printf("Parsing argument : %s\n", argv[i]);
            switch (argv[i][1]) {
                case 104:
                    fprintf(stderr,
                        "=====   GRA   =====\n\n"
                        "Usage: gra [-options] <seqfile>\n\n"
                        "%-12s: list options.\n"
                        "%-12s: scan <seqfile> against <hmmdb>.\n"
                        "%-12s: take the reverse complement of each sequence.\n"
                        "%-12s: set <x> as the E-value threshold. GRA will "
                            "only report sequences with an E-value <= <x>.\n"
                            "%-14sDefaults is 0.9e-45.\n"
                        "%-12s: redirect output to the file <f>.\n%-14sDefault "
                            "is \"<seqfile>_GRA_OUTPUT.fasta\".\n%-14sIf there "
                            "are multiple fasta as input, default is used for "
                            "each.\n"
                        "%-12s: Set the number of parallel worker threads to"
                            " <n>. On multicore machines, the default is 2.\n",
                        "-h", "-p <hmmdb>", "-r", "-e <x>", "", "-o <f>", "", "", "-c <n>"
                    );
                    exit(EXIT_SUCCESS);
                    break;
                case 112:
                    hmm = argv[++i];
                    // printf("HMM profile : %s\n", hmm);
                    break;
                case 101:
                    e_value = argv[++i];
                    // printf("E-value threshold : %s\n", e_value);
                    break;
                case 111:
                    fa_out = argv[++i];
                    // printf("Output file : %s\n", fa_out);
                    break;
                case 99:
                    cpu = argv[++i];
                    // printf("Number of thread : %s\n", cpu);
                    break;
                case 114:
                    rc = "y";
                    break;
                default:
                    printf("Unknown option \"%s %s\", type gra -h to see options.\n",
                        argv[i], argv[i+1]);
                    exit(EXIT_FAILURE);
                    break;
            }
            *j = i + 1;
        }
    }

    return (const char*[]){hmm, e_value, fa_out, cpu, rc};
}

char* write_rc_fasta(const char *fasta) {
    size_t l_fa_rc_name = strlen(fasta) + strlen("_rc") + 1;
    char *fa_rc_name = calloc(l_fa_rc_name, sizeof(char));
    char *fasta_truncated = get_filename(fasta);
    snprintf(fa_rc_name, l_fa_rc_name, "%s_rc.fasta", fasta_truncated);

    FILE *fa_in = fopen(fasta, "r");
    FILE *fa_rc = fopen(fa_rc_name, "w");

    char seq[32768]; seq[0] = '\0';
    char buffer[256];
    size_t s;
    int crlf = 0, d = 1;
    while (fgets(buffer, sizeof(buffer), fa_in) != NULL) {
        if (buffer[0] == '>') {
            s = strlen(buffer);
            if (buffer[s-1] == '\n' && buffer[s-2] == '\r') { crlf = 1; d = 0; }
            if (seq[0] != '\0') {
                char *rev = reverse_complement(seq, crlf);
                fprintf(fa_rc, "%s", rev);
                memset(seq, 0, sizeof(seq));
            }
            buffer[s-2+d] = '/';
            buffer[s-1+d] = 'r';
            buffer[s+d] = 'c';
            if (crlf) { buffer[s+1] = '\r'; }
            buffer[s+2] = '\n';
            buffer[s+3] = '\0';
            fprintf(fa_rc, "%s", buffer);
        } else { strcat(seq, buffer); }
    }
    char *rev = reverse_complement(seq, crlf);
    fprintf(fa_rc, "%s", rev);
    free(rev);

    return fa_rc_name;
}

void run_gra(
    const char *fasta,
    const char *hmm,
    const char *e_value,
    char *fa_out,
    const char *cpu,
    const char *rc
) {
    printf("Args :\nfasta : %s\nhmm : %s\ne_value : %s\nfa_out : %s\ncpu : %s\n",
        fasta, hmm, e_value, fa_out, cpu);

    const size_t n_gene = count_gene(fasta);
    printf("nombre de gène dans le fichier d'entrée : %lu\n", n_gene);

    char* filename = fasta;
    if (rc[0] == 'y') {
        filename = write_rc_fasta(fasta);
    }

    FILE *scan_res = hmmscan(filename, hmm, e_value, cpu);
    printf("Fin hmmscan\n");
    // FILE *f = fopen("vrac/aug_test.txt", "r");
    MapGeneProfile *retrievedGenes = malloc(sizeof(MapGeneProfile) * n_gene);
    size_t n_retrieved_genes = parse_scan_res(scan_res, n_gene, retrievedGenes);
    printf("Fin parsage\n");

    printf("nombre de gène reconnus : %lu\n", n_retrieved_genes);

    fasta_output(filename, retrievedGenes, n_retrieved_genes, fa_out);
    for (int i = 0; i < n_retrieved_genes; i++) {
        free(retrievedGenes[i].gene);
        free(retrievedGenes[i].profile);
        free(retrievedGenes[i].e_value);
    }
    free(retrievedGenes);
    // memset(retrievedGenes, 0, n_gene);
    retrievedGenes = NULL;
}

int main(int argc, char **argv) {
    if (argc < 2 || (argv[argc-2][0] == '-' && argv[argc-2][1] != 'r') ||
        (argv[argc-1][0] == '-' && argv[argc-1][1] != 'h')) {
        fprintf(stderr, "Incorrect number of command line arguments.\n"
        "Usage: gra [-options] <seqfiles>\nType gra -h to see options.\n");
        exit(EXIT_FAILURE);
    }

    size_t j;
    const char **params = get_parameters(argc, argv, &j);
    const char *hmm = params[0];
    const char *e_value = params[1];
    char *fa_out = params[2];
    const char *cpu = params[3];
    const char* rc = params[4];
    char** seqfiles = &argv[j];

    int free_fa_out = 1;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'o') {
            free_fa_out = 0;
        }
    }

    for (int i = 0; i < argc - j; i++) {
        if (argc - j > 1) { fa_out = default_output_name(seqfiles[i]); }
        run_gra(seqfiles[i], hmm, e_value, fa_out, cpu, rc);
        if (free_fa_out) { free(fa_out); }
    }

    return EXIT_SUCCESS;
}
