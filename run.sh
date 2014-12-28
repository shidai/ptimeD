#!/bin/sh

gcc -lm -lgsl -lgslcblas -lfftw3 -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio fitDm.c fdt.c ptimeD.c readfits.c readPred.c T2toolkit.c tempo2pred.c cheby2d.c t1polyco.c -o ptimeD 
