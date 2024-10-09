# Test Data

This directory contains test data that were collected from different sources.
This file attempts to track provenance to facilitate update and licensing.

## HMMs

### `txt`

This folder contains HMMs in ASCII format (either from HMMER2 or HMMER3):

- `PF02826.hmm` was extracted manually from Pfam v33.1.
- `PKSI-AT.hmm2` was obtained from the AntiSMASH repository ([`/detection/hmm_detection/data/PKS_AT.hmm`](https://github.com/antismash/antismash/blob/master/antismash/detection/hmm_detection/data/PKS_AT.hmm))
- `t2pks.hmm` was obtained from the AntiSMASH repository ([`modules/t2pks/data/t2pks.hmm`](https://github.com/antismash/antismash/blob/master/antismash/modules/t2pks/data/t2pks.hmm))

### `bin`

This folder was obtained by converting files from the ASCII folder using
`hmmconvert -b`, provided with HMMER `v3.3.1`.

### `db`

This folder was obtained by pressing files from the ASCII folder using
`hmmpress`, provided with HMMER `v3.3.1`.
