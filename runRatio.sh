#!/bin/bash

rm -f ./dat/summaryRatio.dat
echo "#       Ex        0_1^+        0_2^+        2_1^+        0_6^+" > ./dat/summaryRatio.dat

for i in 19.000 19.100 19.200 19.300 19.400 19.500\
		19.600 19.700 19.800 19.900 20.000\
		20.100 20.200 20.300 20.400 20.500\
		20.600 20.700 20.800 20.900 21.000\
		21.100 21.200 21.300 21.400 21.500\
		21.600 21.700 21.800 21.900 22.000\
		22.100 22.200 22.300 22.400 22.500\
		22.600 22.700 22.800 22.900 23.000\
		23.100 23.200 23.300 23.400 23.500\
		23.600 23.700 23.800 23.900 24.000\
		24.100 24.200 24.300 24.400 24.500\
		24.600 24.700 24.800 24.900 25.000\
		25.100 25.200 25.300 25.400 25.500\
		25.600 25.700 25.800 25.900 26.000\
		26.100 26.200 26.300 26.400 26.500\
		26.600 26.700 26.800 26.900 27.000\
		27.100 27.200 27.300 27.400 27.500\
		27.600 27.700 27.800 27.900 28.000\
		28.100 28.200 28.300 28.400 28.500\
		28.600 28.700 28.800 28.900 29.000\
		29.100 29.200 29.300 29.400 29.500\
		29.600 29.700 29.800 29.900 30.000\
		30.100 30.200 30.300 30.400 30.500\
		30.600 30.700 30.800 30.900 31.000\
		31.100 31.200 31.300 31.400 31.500\
		31.600 31.700 31.800 31.900 32.000\
		32.100 32.200 32.300 32.400 32.500\
		32.600 32.700 32.800 32.900 33.000\
		33.100 33.200 33.300 33.400 33.500\
		33.600 33.700 33.800 33.900 34.000\
		34.100 34.200 34.300 34.400 34.500\
		34.600 34.700 34.800 34.900 35.000
do

    root -l -b<<EOF
.x ./macro/smeared.C($i)
.q
EOF
    
done

