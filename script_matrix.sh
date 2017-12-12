#!/bin/sh
qrsh -l cuda=TRUE
nvcc cudamatrix2.c -o cudamatrix2
cudamatrix2
