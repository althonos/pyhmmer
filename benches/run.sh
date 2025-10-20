#!/bin/sh

N=$(nproc)

cp src/pyhmmer/tests/data/seqs/CP040672.1.genes_100.fna /tmp/CP040672.1.genes_100.fna
cp src/pyhmmer/tests/data/seqs/CP000560.2.fna /tmp/CP000560.2.fna

if [ ! -e "/tmp/562.PRJEB4685.faa" ]; then
	wget "https://progenomes.embl.de/dumpSequence.cgi?p=562.PRJEB4685&t=p&a=562" -O- | gunzip > /tmp/562.PRJEB4685.faa
fi
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
    'python -m pyhmmer.hmmer --jobs {threads} hmmsearch /tmp/Pfam.v33-1.hmm /tmp/562.PRJEB4685.faa' \
    'python -m pyhmmer.hmmer --jobs {threads} hmmsearch /tmp/Pfam.v33-1.pressed.hmm /tmp/562.PRJEB4685.faa' \
    'python -m pyhmmer.hmmer --jobs {threads} hmmscan /tmp/Pfam.v33-1.pressed.hmm /tmp/562.PRJEB4685.faa' \
    'hmmsearch --cpu {threads} /tmp/Pfam.v33-1.hmm /tmp/562.PRJEB4685.faa' \
    'hmmsearch --cpu {threads} /tmp/Pfam.v33-1.pressed.hmm /tmp/562.PRJEB4685.faa' \
    'hmmscan --cpu {threads} /tmp/Pfam.v33-1.pressed.hmm /tmp/562.PRJEB4685.faa' \
    --export-json benches/hmmsearch.json -r 3

hyperfine -P threads 1 $N \
    'python -m pyhmmer.hmmer --jobs {threads} phmmer /tmp/562.PRJEB4685.faa /tmp/562.PRJEB4685.faa' \
    'phmmer --cpu {threads} /tmp/562.PRJEB4685.faa /tmp/562.PRJEB4685.faa' \
    --export-json benches/phmmer.json -r 3

hyperfine -P threads 1 $N \
    'python -m pyhmmer.hmmer --jobs {threads} nhmmer /tmp/CP040672.1.genes_100.fna /tmp/CP000560.2.fna' \
    'nhmmer --cpu {threads} /tmp/CP040672.1.genes_100.fna /tmp/CP000560.2.fna' \
    --export-json benches/nhmmer.json -r 3
