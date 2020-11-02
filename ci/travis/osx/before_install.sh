#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh


# --- Fix HMMER makefile -----------------------------------------------------

sed -i.bak "s/= @CFLAGS@/= @CFLAGS@ @PIC_CFLAGS@/g" vendor/hmmer/src/impl_sse/Makefile.in
