#!/usr/bin/env bash

mkdir -p build
gcc vcm.c -std=c99 -Wall -Ofast -march=native -fopenmp -o build/vcm -lm
