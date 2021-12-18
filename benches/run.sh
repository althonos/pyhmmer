#!/bin/sh

N=$(nproc)

hyperfine -P threads 1 $N \
    'python -m pyhmmer.hmmer --jobs {threads} hmmsearch Pfam.v33-1.hmm tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    'python -m pyhmmer.hmmer --jobs {threads} hmmsearch Pfam.v33-1.pressed.hmm tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    'python -m pyhmmer.hmmer --jobs {threads} hmmsearch Pfam.v33-1.pressed.hmm.h3m tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    'hmmsearch --cpu {threads} Pfam.v33-1.hmm tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    'hmmsearch --cpu {threads} Pfam.v33-1.pressed.hmm tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    'hmmsearch --cpu {threads} Pfam.v33-1.pressed.hmm.h3m tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    'hmmscan --cpu {threads} Pfam.v33-1.pressed.hmm tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    --export-json benches/hmmsearch8.json -r 4 -w 1

hyperfine -P threads 1 $N \
    'python -m pyhmmer.hmmer --jobs {threads} phmmer tests/data/seqs/938293.PRJEB85.HG003687.faa tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    'phmmer --cpu {threads} tests/data/seqs/938293.PRJEB85.HG003687.faa tests/data/seqs/938293.PRJEB85.HG003687.faa' \
    --export-json benches/phmmer.json -r 4 -w 1

hyperfine -P threads 1 $N \
    'python -m pyhmmer.hmmer --jobs {threads} nhmmer tests/data/seqs/CP040672.1.genes_100.fna tests/data/seqs/CP000560.2.fna' \
    'nhmmer --cpu {threads} tests/data/seqs/CP040672.1.genes_100.fna tests/data/seqs/CP000560.2.fna' \
    --export-json benches/nhmmer.json -r 4 -w 1
