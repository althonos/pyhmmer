#!/bin/sh

N=$(nproc)

cp pyhmmer/tests/data/seqs/CP040672.1.genes_100.fna /tmp/CP040672.1.genes_100.fna
cp pyhmmer/tests/data/seqs/CP000560.2.fna /tmp/CP000560.2.fna
cp pyhmmer/tests/data/seqs/938293.PRJEB85.HG003687.faa /tmp/938293.PRJEB85.HG003687.faa

if [ ! -e "/tmp/Pfam.v33-1.pressed.hmm" ]; then
	wget "http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz" -O- | gunzip > /tmp/Pfam.v33-1.pressed.hmm
fi
if [ ! -e "/tmp/Pfam.v33-1.pressed.hmm.h3i" ]; then
	hmmpress /tmp/Pfam.v33-1.pressed.hmm
fi
if [ ! -e "/tmp/Pfam.v33-1.hmm" ]; then
	ln -s /tmp/Pfam.v33-1.pressed.hmm /tmp/Pfam.v33-1.hmm
fi

hyperfine -P threads 1 $N \
    'python -m pyhmmer.hmmer --jobs {threads} hmmsearch /tmp/Pfam.v33-1.hmm /tmp/938293.PRJEB85.HG003687.faa' \
    'python -m pyhmmer.hmmer --jobs {threads} hmmsearch /tmp/Pfam.v33-1.pressed.hmm /tmp/938293.PRJEB85.HG003687.faa' \
    'python -m pyhmmer.hmmer --jobs {threads} hmmsearch /tmp/Pfam.v33-1.pressed.hmm.h3m /tmp/938293.PRJEB85.HG003687.faa' \
    'hmmsearch --cpu {threads} /tmp/Pfam.v33-1.hmm /tmp/938293.PRJEB85.HG003687.faa' \
    'hmmsearch --cpu {threads} /tmp/Pfam.v33-1.pressed.hmm /tmp/938293.PRJEB85.HG003687.faa' \
    'hmmsearch --cpu {threads} /tmp/Pfam.v33-1.pressed.hmm.h3m /tmp/938293.PRJEB85.HG003687.faa' \
    'hmmscan --cpu {threads} /tmp/Pfam.v33-1.pressed.hmm /tmp/938293.PRJEB85.HG003687.faa' \
    --export-json benches/hmmsearch.json -r 4 -w 1

hyperfine -P threads 1 $N \
    'python -m pyhmmer.hmmer --jobs {threads} phmmer /tmp/938293.PRJEB85.HG003687.faa /tmp/938293.PRJEB85.HG003687.faa' \
    'phmmer --cpu {threads} /tmp/938293.PRJEB85.HG003687.faa /tmp/938293.PRJEB85.HG003687.faa' \
    --export-json benches/phmmer.json -r 4 -w 1

hyperfine -P threads 1 $N \
    'python -m pyhmmer.hmmer --jobs {threads} nhmmer /tmp/CP040672.1.genes_100.fna /tmp/CP000560.2.fna' \
    'nhmmer --cpu {threads} /tmp/CP040672.1.genes_100.fna /tmp/CP000560.2.fna' \
    --export-json benches/nhmmer.json -r 4 -w 1
