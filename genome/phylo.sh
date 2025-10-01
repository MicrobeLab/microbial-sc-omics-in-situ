#!/usr/bin/env bash

phylophlan -i /path/to/genomes \
    -d /path/to/phylophlan \
    --diversity low \
    -f supermatrix_aa.cfg \
    -o /path/to/output 
    --nproc 16  
