#!/bin/bash

gcc -O2 -g characteristic.c lambertw.c vector.c -lm -o compiled/characteristic
time ./compiled/characteristic
python visualize.py psi.txt
