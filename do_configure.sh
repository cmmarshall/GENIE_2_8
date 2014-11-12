#!/bin/bash

./configure \
  --enable-doxygen-doc \
  --enable-debug \
  --enable-test \
  --enable-numi \
  --enable-gsl \
  --with-optimiz-level=O0 \
  --with-doxygen-path=/usr/bin/doxygen \
  --with-log4cpp-inc=$LOG4CPP_INC \
  --with-log4cpp-lib=$LOG4CPP_LIB \
  --with-libxml2-inc=/usr/include/libxml2 \
  --with-libxml2-lib=/usr/lib64 

