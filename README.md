#  pyHMMER [![Stars](https://img.shields.io/github/stars/althonos/pyhmmer.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyhmmer/stargazers)


<!-- ## 🟡♦️🟦 Overview -->


## ⏱️ Benchmarks

Benchmarks were run on a i7-8550U CPU running at 1.80GHz, using a FASTA file
containing 2100 protein sequences (see `tests/data/seqs/938293.PRJEB85.HG003687.faa`)
and a subset of the PFam HMM library containing 2873 domains. Both commands
were set to use 8 threads, and were run 20 times.

| Command                     | mean (s) | σ (ms) | min (s) | max (s) |
|-----------------------------|----------|--------|---------|---------|
| `hmmsearch`                 | 40.723   | 223    | 40.402  | 41.321  |
| `python -m hmmer.hmmsearch` | 40.128   | 946    | 37.706  | 42.457  |


## 📜 License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The HMMER3 and Easel code is available under the BSD 3-clause license. See
`vendor/hmmer/LICENSE` and `vendor/easel/LICENSE` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the [original HMMER authors](http://hmmer.org/). It was developed by
[Martin Larralde](https://github.com/althonos/pyhmmer) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
