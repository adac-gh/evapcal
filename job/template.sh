#!/bin/bash
#PBS -q AM

cd $PBS_O_WORKDIR
cd ../

./evapcal 100000 23.500
