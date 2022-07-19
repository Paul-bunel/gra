import sys
import re
from statistics import *
from Bio import SeqIO

def stats_evalue(res_file):
    genes = ["TAS1R1", "TAS1R2", "TAS1R3", "TAS2R"]
    genes_evalue = {gene: [] for gene in genes}
    
    for record in SeqIO.parse(res_file, "fasta"):
        header = record.description.split('|')
        genes_evalue[header[-2]].append(float(header[-1]))

    for gene in genes:
        print(gene, ':')
        print("\tMean :", mean(genes_evalue[gene]))
        print("\tMin :", min(genes_evalue[gene]))
        print("\tMax :", max(genes_evalue[gene]))


def stats_results(res_file, *in_files):
    """Vérifie que les gènes sont biens attribués au bon profil HMM
    Pour cela : parcours les fichier d'entrée et vérifie que chaque gène est
    retrouvé et annoté correctement dans la sortie"""
    gene_res = {}
    fp = 0
    fn = 0
    tp = tn = 0
    for record in SeqIO.parse(res_file, "fasta"):
        accession = record.id
        profile, evalue = record.description.split('|')[-2:]
        if accession in gene_res:
            print("%s est déjà dans le dictionnaire !" % accession)
        gene_res[accession] = (profile, evalue)

    # [print(f"{key :<30}: {value}") for key, value in gene_res.items()]        
    dico = {}
    for file in in_files:
        dico[file] = 0
        profile = re.search("TAS[12]R[123]*", file)
        if profile:
            profile = profile[0]
        for record in SeqIO.parse(file, "fasta"):
            accession = record.id
            if profile:
                if accession in gene_res:# and float(gene_res[accession][1]) < 0.9e-45:
                    if gene_res[accession][0] != profile:
                        print("%s de profil %s a été reconnu au mauvais profil : %s" %
                            (accession, profile, gene_res[accession]))
                        fn += 1
                    else:
                        tp += 1
                else:
                    print("%s de profil %s n'a pas été reconnu" % (accession, profile))
                    fn += 1
            else:
                if accession in gene_res:# and float(gene_res[accession][1]) < 0.9e-45:
                    # print(record.description, gene_res[accession])
                    # print("Le gene", accession, ",ne devant pas être reconnu, à",
                        # "été identifié au profil", gene_res[accession], "par erreur")
                    dico[file] += 1
                    fp += 1
                else:
                    tn += 1

    print(dico)
    print("Nombre de faux positifs : %s, nombre de faux négatifs : %s" % (fp, fn),
          "Nombre de vrais positifs : %s, nombre de vrais négatifs : %s" % (tp, tn))



with open(sys.argv[1]) as file:
    stats_results(
        file,
        "SEQUENCES/TEST/TAS1R1_RefSeq_Predicted_Pseudo.fasta",
        "SEQUENCES/TEST/TAS1R2_RefSeq_Predicted_Pseudo.fasta",
        "SEQUENCES/TEST/TAS1R3_RefSeq_Predicted_Pseudo.fasta",
        "SEQUENCES/TEST/TAS2R_RefSeq_Predicted_Pseudo.fasta",
        "SEQUENCES/TEST/TAS1R1_INSDC_clean.fasta",
        "SEQUENCES/TEST/TAS1R2_INSDC_clean.fasta",
        "SEQUENCES/TEST/TAS1R3_INSDC_clean.fasta",
        "SEQUENCES/TEST/TAS2R_INSDC_clean.fasta",
        "SEQUENCES/TEST/Vomeronasal_RefSeq_Mammalia.fasta",
        "SEQUENCES/TEST/OR_cetaceans.fasta",
        "SEQUENCES/TEST/Megaptera_novaeangliae_seqs_nTAS.fasta",
        "SEQUENCES/TEST/Tursiop_truncatus_seqs_nTAS.fasta",
        "SEQUENCES/TEST/Globicephala_melas_seqs_nTAS.fasta"
    )