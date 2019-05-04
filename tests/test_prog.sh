#!/bin/bash

../cancer_gillespie_simulation -o output/test1_ -s 123456
../cancer_gillespie_simulation -o output/test2_ -s 123456 -P20,20
../cancer_gillespie_simulation -o output/test3_ -s 123456 -P20,20 -M100
../cancer_gillespie_simulation -o output/test4_ -s 123456 -X2 -M100 -V0.05 -R5 -D200
