#!/bin/sh -e

. $(dirname $(dirname $0))/functions.sh

# --- Launch manylinux container ---------------------------------------------

IMG="pypywheels/manylinux2010-pypy_x86_64"
log Launching \`manylinux\` docker container
docker pull $IMG
docker run -d -e TERM=$TERM -v $TRAVIS_BUILD_DIR:/io --name manylinux --rm -it $IMG sh


# --- Fix HMMER makefile -----------------------------------------------------

sed -i "s/= @CFLAGS@/= @CFLAGS@ @PIC_CFLAGS@/g" vendor/hmmer/src/impl_sse/Makefile.in
