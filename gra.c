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

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Utilisation : %s fichier1.fasta fichier2.fasta ...\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    struct String fasta = getFiles(argc - 1, &argv[1]);

    printf("Réalisation de la commande bash...\n");
    int l_cmd = strlen("/bin/bash -c 'head <<< \"") + fasta.len + 3;
    char *cmd = calloc(l_cmd, sizeof(char));
    strcpy(cmd, "/bin/bash -c 'head <<< \"");
    strcat(cmd, fasta.str);

    cmd[l_cmd - 3] = '\'';
    cmd[l_cmd - 2] = '"';
    cmd[l_cmd - 1] = '\0';

    printf("Commande : \n");
    int status = system(cmd);

    return EXIT_SUCCESS;
}

// caca\0 -> 4
// '"caca"'\0

