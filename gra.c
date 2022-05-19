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

char** split_string(char *str, char *sep) {
    size_t count_tokens = 0;
    int c = 0;
 	for (int i = 0; i < strlen(str) - 1; i++) {
        if (str[i] >= 33 && str[i] < 127) { c = 1; }
  		if (c && str[i] == sep[0] && str[i + 1] != sep[0] && str[i + 1] != '\0') {
  			count_tokens++;
 		}
	}
    char **tokens = malloc(++count_tokens * sizeof(char*));
    printf("Nombre token = %i\n", count_tokens);

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

FILE* hmmscan(char *profile, char *fasta) {
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
void parse_scan_res(FILE *scan_res) {
    char buffer[512];
    while (fgets(buffer, sizeof(buffer), scan_res) != NULL) {
        if (buffer[0] != '#') {
            printf("%s", buffer);
            char **tokens = split_string(buffer, " ");
            
            // int i = 0;
            char *token;
            // while (tokens[i] != NULL) {
            for (int i = 0; i < 19; i++) {
                token = strdup(tokens[i]);
                printf("Token n°%i = %s\n", i, token);
            }
            // printf("allo\n");
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Utilisation : %s fichier1.fasta fichier2.fasta ...\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // struct String fasta = getFiles(argc - 1, &argv[1]);
    char* fasta = argv[1];

    FILE *scan_res = hmmscan("HMM_PROFILE/TAS/TAS_ncbi_nuc.hmm", fasta);
    parse_scan_res(scan_res);

    // printf("Réalisation de la commande bash...\n");
    // int l_cmd = strlen("/bin/bash -c 'head <<< \"") + fasta.len + 2;
    // char *cmd = calloc(l_cmd, sizeof(char));
    // strcpy(cmd, "/bin/bash -c 'head <<< \"");
    // strcat(cmd, fasta.str);

    // cmd[l_cmd - 3] = '"';
    // cmd[l_cmd - 2] = '\'';
    // cmd[l_cmd - 1] = '\0';

    // printf("Commande : \n");

    // for (int i = 10; i >= 1; i--) {
    //     printf("cmd[l_cmd - %d] = %i\n", i, cmd[l_cmd - i]);
    // }
    // int status = system(cmd);

    return EXIT_SUCCESS;
}

// caca\0 -> 4
// '"caca"'\0

