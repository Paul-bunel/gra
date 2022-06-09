#!/bin/bash
# utilisation : ./make_profil.sh <msafile_TAS.fasta>
# sortie : HMM_PROFILE/TAS/TASx/<TASx.hmm>
SRC=$1

name=$(sed 's|.*/||' <<< $SRC | sed 's/msafile_//' | sed 's/\..*//')
printf "File : $name\n"
gene=$(grep -o 'TAS[12]R[123]*' <<< $name)
printf "Gene : $gene\n\n"

TAS=''
[[ $name == *TAS* ]] && TAS='/TAS'

[ ! -d "HMM_PROFILE$TAS/${gene^^}" ] && mkdir "HMM_PROFILE$TAS/${gene^^}" &&
    printf "Creating directory 'HMM_PROFILE$TAS/${gene^^}'\n\n"

printf "Building HMM profile...\n"
printf "hmmbuild -n $gene -o \"/dev/null\" --dna \"HMM_PROFILE$TAS/${gene^^}/$name.hmm\" $SRC\n\n"
hmmbuild -n $gene -o "/dev/null" --dna "HMM_PROFILE$TAS/${gene^^}/$name.hmm" $SRC

printf "Converting to binaries files\n"
printf "hmmpress -f \"HMM_PROFILE$TAS/${gene^^}/$name.hmm\" > /dev/null\n"
hmmpress -f "HMM_PROFILE$TAS/${gene^^}/$name.hmm" > /dev/null
