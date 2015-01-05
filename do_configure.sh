#!/bin/bash

./configure \
  --enable-validation-tools \
  --enable-test \
  --enable-debug \
  --enable-gsl \
  --with-optimiz-level=O0 \
  --enable-numi \
  --enable-rwght \
  --with-log4cpp-inc=$LOG4CPP_INC \
  --with-log4cpp-lib=$LOG4CPP_LIB \
  --with-libxml2-inc=/usr/include/libxml2 \
  --with-libxml2-lib=/usr/lib64 \
  >& log.config

  # --enable-vle-extension \
  # --enable-doxygen-doc \
  # --with-doxygen-path=/usr/bin/doxygen \
  # --disable-lhapdf \
  # --enable-cernlib \
