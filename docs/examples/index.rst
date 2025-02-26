Examples
========

Analyses
--------

This section of the documentation shows several examples adapted from real
examples, and running with the latest version of the PyHMMER interface.

.. toctree::
   :maxdepth: 2
   :caption: Analyses

   Build an HMM from an multiple sequence alignment <msa_to_hmm>
   Analyse the active site of an enzymatic domain <active_site>
   Fetch marker genes from a genome <fetchmgs>
   Run an iterative search to build a HMM for rhodopsins <iterative_search>


Code & Data
-----------

This section shows more practical tutorials about how you can use the
PyHMMER API in combination with Python and other Python libraries:

.. toctree::
    :maxdepth: 2
    :caption: Code & Data

    Recipes for extending PyHMMER <recipes>
    Create a Python package with embedded HMMs <embed_hmms>
    Improve performance with some tips <performance_tips>


Dependent projects
------------------

PyHMMER is being used in several projects, including:

- `GECCO <https://gecco.embl.de>`_, a tool for detecting Biosynthetic Gene
  Clusters in genomic data, uses PyHMMER to annotate proteins with Pfam domains
  to use as sequence features for a machine learning model.
- `BiG-SLICE <https://github.com/medema-group/bigslice>`_, an interactive tool 
  for the large scale analysis of Biosynthetic Gene Clusters data, uses
  PyHMMER for the Pfam domain annotation step.
- `duomolog <https://github.com/tijeco/duomolog>`_, a method to identify the
  best set of homologous sequences from two homology searching approaches,
  uses PyHMMER to search for homologous sequences before comparing the results
  to BLAST hits.
- `cinful <https://github.com/wilkelab/cinful>`_, a pipeline to identify
  microcins along with their associated immunity proteins and export machinery,
  uses PyHMMER to build HMMs from MAFFT-derived alignments, and to run
  `~pyhmmer.hmmer.hmmsearch`.
- `BUSCOlite <https://github.com/nextgenusfs/buscolite>`_, a package for
  running simplified `BUSCO <https://busco.ezlab.org/>`_ analysis for gene
  prediction, uses PyHMMER to discover Conserved Orthologous Genes.
- `FastAAI <https://github.com/cruizperez/FastAAI>`_, a package providing
  fast estimation of Average Amino Acid Identities (AAI) for bacterial and
  viral genomes, uses PyHMMER to detect bacterial and archaeal domains.
- `SADIE <https://github.com/jwillis0720/sadie>`_, the Sequencing Analysis and
  Data library for Immunoinformatics Exploration, uses PyHMMER as an aligner
  backend to re-number sequence alignment columns.
- `GSub <https://github.com/FlorianCHA/Gsub>`_, a GUI for submitting viral
  sequences to `GenBank <https://www.ncbi.nlm.nih.gov/genbank/>`_, uses
  PyHMMER to detect potential polymerase-encoding ORFs.
- `Bakta <https://github.com/oschwengers/bakta>`_, a standardised annotation
  tool for bacterial genomes and plasmids, uses PyHMMER to annotate genes
  with different HMM libraries such as
  `AntiFam <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3308159/>`_
  or `Pfam <https://academic.oup.com/nar/article/49/D1/D412/5943818>`_.
- `Pharokka <https://github.com/gbouras13/pharokka>`_, a rapid standardised
  annotation tool for bacteriophage genomes and metagenomes, uses PyHMMER
  to annotate genes with the `PHROGS <https://phrogs.lmge.uca.fr/>`_ HMMs.
- `VEBA <https://github.com/jolespin/veba>`_, A modular end-to-end suite
  for in silico recovery, clustering, and analysis of prokaryotic, microeukaryotic,
  and viral genomes from metagenomes.  Uses PyHMMER to annotate proteins with
  homology to `Pfam`_, KOfam, NCBIfam-AMR, and `AntiFam`_.
- `PyHMMSearch <https://github.com/jolespin/pyhmmsearch>`_, Fast implementation of
  HMMSEARCH optimized for high-memory systems using PyHMMER.
- `PyKOfamSearch <https://github.com/jolespin/pyhmmsearch>`_, Fast implementation of
  KofamScan optimized for high-memory systems using PyHMMER.

.. hint::

    If you use PyHMMER in a scientific context, feel free to open a
    `pull request <https://github.com/althonos/pyhmmer/>`_ to add yourself
    to the list!


External resources
------------------

The `pytrimal <https://pypi.org/project/pytrimal>`_ documentation contains an
example for building an HMM with PyHMMER from a trimmed alignment, check it
`here <https://pytrimal.readthedocs.io/en/stable/examples/hmmer.html>`_.

For a more detailed explanation of HMMER features, you should also check
the `HMMER User's Guide <http://eddylab.org/software/hmmer/Userguide.pdf>`_.
