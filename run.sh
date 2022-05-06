#!/bin/bash
ALIGNMENTS='SEQUENCES/ALIGNMENTS'

for file in $ALIGNMENTS/*
do
    name=$(sed 's|.*/||' <<< $file | sed 's/msafile_//' | sed 's/\..*//')
    echo "Database : $name"
    
    [ ! -d "HMM_PROFILE/${name^^}" ] && mkdir "HMM_PROFILE/${name^^}" &&
        echo "Creating directory 'HMM_PROFILE/${name^^}'"

    echo "Building HMM profile..."
    hmmbuild -o "/dev/null" --amino "HMM_PROFILE/${name^^}/$name.hmm" $file

    echo "Converting to binaries files"
    hmmpress -f "HMM_PROFILE/${name^^}/$name.hmm" > /dev/null

    [ ! -d "RESULTS/${name^^}" ] && mkdir "RESULTS/${name^^}" &&
        echo "Creating directory 'RESULTS/${name^^}'"


    echo "Scanning sources databases on the profile..."
    for db in $ALIGNMENTS/*
    do
        [[ $db == *.fasta ]] && format='fasta'
        [[ $db == *.afa ]] && format='afa'
        [[ $db == *.sto ]] && format='stockholm'

        DBname=$(sed 's|.*/||' <<< $db | sed 's/msafile_//' | sed 's/\..*//')
        hmmscan -o "/dev/null" \
            --tblout "RESULTS/${name^^}/res_$DBname.txt" \
            --qformat $format \
            "HMM_PROFILE/${name^^}/$name.hmm" \
            $db
        # hmmscan -o "/dev/null" --pfamtblout "RESULTS/${name^^}/res_pfam.txt" \
        #     "HMM_PROFILE/${name^^}/$name.hmm" \
        #     $file
    done
done