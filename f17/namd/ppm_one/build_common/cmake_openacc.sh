#!/bin/bash

rm -rf CMakeCache.txt
rm -rf CMakeFiles

SOURCE=../common

cmake \
    -DCMAKE_C_COMPILER:STRING=pgcc \
    -DCMAKE_CXX_COMPILER:STRING=pgc++ \
    -DCMAKE_C_FLAGS:STRING="-Minfo=accel -tp=haswell -ta=multicore" \
    -DCMAKE_CXX_FLAGS:STRING="-Minfo=accel -tp=haswell -ta=multicore" \
    $SOURCE
