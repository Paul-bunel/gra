# gra
Gustatory Receptor family Assigner

Usage: gra [-options] <seqfile>  
  
-h          : list options.  
-p <hmmdb>  : scan <seqfile> against <hmmdb>.  
-r          : take the reverse complement of each sequence.  
-e <x>      : set <x> as the E-value threshold. GRA will only report sequences with an E-value <= <x>.  
              Defaults is 0.9e-45.  
-o <f>      : redirect output to the file <f>.  
              Default is "<seqfile>_GRA_OUTPUT.fasta".  
              If there are multiple fasta as input, default is used for each.  
-c <n>      : Set the number of parallel worker threads to <n>. On multicore machines, the default is 2.  

