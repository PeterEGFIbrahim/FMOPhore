#!/bin/bash

#Create directories for each input and move the input there
for i in $(ls *inp); do
mkdir -p $(echo ${i} | awk -F '[.]' '{print $1}') 
#mkdir $(echo ${i} | awk -F '[.]' '{print $1}') 
mv ${i} $(echo ${i} | awk -F '[.]' '{print $1}')/.
cp DFTB $(echo ${i} | awk -F '[.]' '{print $1}')/.
cd $(echo ${i} | awk -F '[.]' '{print $1}')/
sed -i "s/-N /-N x$(echo ${i} | awk -F '[.]' '{print $1}')/g" DFTB
sed -i "s/input.inp/${i}/g" DFTB
sed -i "s/output.log/$(basename ${i} .inp).log/g" DFTB
qsub -jc long DFTB
cd ../
done
