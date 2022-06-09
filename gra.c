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

FILE* hmmscan(const char *profile, const char *fasta) {
    /**
        hmmscan -o "/dev/null" \
        --tblout "RESULTS/${name^^}/res_$DBname.txt" \
        --qformat $format \
        "HMM_PROFILE/${name^^}/$name.hmm" \
        $db
     */
    char *fun = "/bin/bash -c '"
        "hmmscan -o /dev/null "
        "-E 0.9e-45 "
        "--tblout /dev/stdout "
        "--qformat fasta ";

    size_t l_cmd = strlen(fun) + strlen(profile) + strlen(fasta) + 3;
    char *cmd = calloc(l_cmd, sizeof(char));
    snprintf(cmd, l_cmd, "%s%s %s'", fun, profile, fasta);
    cmd[l_cmd - 1] = '\0';

    
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
MapGeneProfile* parse_scan_res(FILE *scan_res, const size_t n_gene) {
    char buffer[512];
    MapGeneProfile *retrievedGenes = malloc(sizeof(MapGeneProfile) * n_gene);

    size_t cpt = 0;
    while (fgets(buffer, sizeof(buffer), scan_res) != NULL) {
        if (buffer[0] != '#') {
            Tokens tokens = split_string(buffer, " ");
            int retrieved = 0;
            
            for (int i = 0; i < n_gene; i++) {
                if (retrievedGenes[i].gene != NULL &&
                    strcmp(retrievedGenes[i].gene, tokens.tokens[2]) == 0) {
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
    return retrievedGenes;
}

char* get_filename(const char* filename) {
    Tokens tokens = split_string(filename, "/");
    int i;
    Tokens raw_filename = split_string(tokens.tokens[tokens.size - 1], ".");

    free(tokens.tokens);
    return raw_filename.tokens[0];
}

/**
 * @brief fonction écrivant dans un nouveau fichier fasta les gènes reconnus
 * par les profils HMM, annoté avec le nom du profil l'ayant reconnu.
 * @param fasta le fichier fasta d'entrée contenant tous les gènes à tester
 * @param retrievedGenes les gènes reconnus
 * @param n_gene le nombre de gènes dans le fichier fasta d'entrée
 */
void fasta_output(const char* fasta, MapGeneProfile *retrievedGenes, size_t n_gene) {
    FILE *fa_in = fopen(fasta, "r");

    char *fasta_truncated = get_filename(fasta);
    size_t l_fa_out_name = strlen(fasta_truncated) + strlen("_GRA_OUTPUT.fasta") + 1;
    char *fa_out_name = calloc(l_fa_out_name, sizeof(char));
    snprintf(fa_out_name, l_fa_out_name, "%s_GRA_OUTPUT.fasta", fasta_truncated);
    FILE *fa_out = fopen(fa_out_name, "w");

    char seq[32768]; seq[0] = '\0';
    char buffer[256];
    int read = 0;
    while (fgets(buffer, sizeof(buffer), fa_in) != NULL) {
        if (buffer[0] == '>') {
            fprintf(fa_out, "%s", seq);
            memset(seq, 0, sizeof(seq));
            read = 0;
            for (int i = 0; i < n_gene; i++) {
                if (retrievedGenes[i].gene != NULL &&
                    strstr(buffer, retrievedGenes[i].gene) != NULL) {
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
    fprintf(fa_out, "%s", seq);
    free(fa_out_name);
    free(retrievedGenes);
    fclose(fa_in);
    fclose(fa_out);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Utilisation : %s <fichier.fasta> <fichier.hmm> ...\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    const char* fasta = argv[1];
    const char* hmm = argv[2];
    const size_t n_gene = count_gene(fasta);
    printf("nombre de gène dans le fichier d'entrée : %lu\n", n_gene);

    FILE *scan_res = hmmscan(hmm, fasta);
    MapGeneProfile *retrievedGenes = parse_scan_res(scan_res, n_gene);

    size_t n_retrieved_genes = 0;
    for (int i = 0; i < n_gene; i++) {
        if (retrievedGenes[i].gene != NULL) { n_retrieved_genes++; }
    }
    printf("nombre de gène reconnus : %lu\n", n_retrieved_genes);

    fasta_output(fasta, retrievedGenes, n_gene);

    return EXIT_SUCCESS;
}
