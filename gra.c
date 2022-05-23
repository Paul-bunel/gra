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

typedef struct MapGeneProfile {
    char *gene;
    char *profile;
    char *e_value;
} MapGeneProfile;

char** split_string(char *str, char *sep) {
    size_t count_tokens = 0;
    int c = 0;
 	for (int i = 0; i < strlen(str) - 1; i++) {
        if (!c && str[i] >= 33 && str[i] < 127) { c = 1; }
  		if (c && str[i] == sep[0] && str[i + 1] != sep[0] && str[i + 1] != '\0') {
  			count_tokens++;
 		}
	}
    char **tokens = malloc(++count_tokens * sizeof(char*));
    // printf("Nombre token = %d\n", count_tokens);

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

    // for (int i = 0; i < count_tokens; i++) {
    //     printf("tokens n°%i : %s\n", i, tokens[i]);
    // }

    free(s);

    return tokens;
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
            printf("Buffer : %s", buffer);
            char **tokens = split_string(buffer, " ");
            int retrieved = 0;
            
            for (int i = 0; i < n_gene; i++) {
                if (retrievedGenes[i].gene != NULL &&
                    strcmp(retrievedGenes[i].gene, tokens[2]) == 0) {
                        retrieved = 1;
                }
            }

            if (!retrieved) {
                retrievedGenes[cpt++] = (MapGeneProfile) {
                    .gene = tokens[2],
                    .profile = tokens[0],
                    .e_value = tokens[4]
                };
            }
        }
    }

    fclose(scan_res);
    return retrievedGenes;
}

void fasta_output(const char* fasta, MapGeneProfile *retrievedGenes, size_t n_gene) {
    FILE *fd = fopen(fasta, "r");
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), fd) != NULL) {
        if (buffer[0] == '>') {
            for (int i = 0; i < n_gene; i++) {
                if (retrievedGenes[i].gene != NULL &&
                    strstr(buffer, retrievedGenes[i].gene) != NULL) {
                    printf("gene %s retrouvé sur le header : %s\n", retrievedGenes[i].gene, buffer);
                    // TODO: Récupérer la séquence et l'écrire dans un nouveau fichier fasta
                }
            }
        }
    }

    fclose(fd);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Utilisation : %s fichier1.fasta fichier2.fasta ...\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // struct String fasta = getFiles(argc - 1, &argv[1]);
    const char* fasta = argv[1];
    const size_t n_gene = count_gene(fasta);

    FILE *scan_res = hmmscan("HMM_PROFILE/TAS/TAS_ncbi_nuc.hmm", fasta);
    MapGeneProfile *retrievedGenes = parse_scan_res(scan_res, n_gene);

    for (int i = 0; i < n_gene; i++) {
        if (retrievedGenes[i].gene != NULL) {
            printf("gène n°%d :%s, profile : %s, e_value : %s\n",
                i, retrievedGenes[i].gene, retrievedGenes[i].profile,
                retrievedGenes[i].e_value);
        }
    }

    fasta_output(fasta, retrievedGenes, n_gene);

    return EXIT_SUCCESS;
}
