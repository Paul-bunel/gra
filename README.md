# gra
## Gustatory Receptor family Assigner  

Utilisation : gra \[-options\] \<seqfile\>  
  
-h          : liste des options.  
-p \<hmmdb\>  : scannne \<seqfile\> contre \<hmmdb\>.  
-r          : prend le complément inverse de chaque séquence de \<seqfile\>.  
-e \<x\>      : pose \<x\> comme seuil de E-value. GRA reportera uniquement les  
              séquences avec une E-value <= \<x\>. La valeur par défault est 0.9e-45.  
-o \<f\>      : redirige la sortie vers le fichier \<f\>.  
              Par défaut, écrit dans "\<seqfile\>_GRA_OUTPUT.fasta".  
              S'il y a plusieurs fichier \<seqfile\>, le défaut est alors automatiquement utilisé.  
-c \<n\>      : Place le nombre de thread parallèles à \<n\>. Sur les machines multicoeurs, la valeur  
              par défaut est deux.  

## Description des différents sous-dossiers du dossier gra

* documentation : documentation générée automatiquement par doxygen
* HMM_PROFILE : contient les différents profils HMM utilisés par gra
* scripts : contient les différents scripts utilisés pendant la réalisation de gra.
* SEQUENCES : les séquences ayant servi à générer les profils HMM, et celles ayant servi à les tester
* vrac : différents fichiers et dossier n'étant pas utile pour gra. (Notamment des résultats de test, d'anciens profils, etc.)

