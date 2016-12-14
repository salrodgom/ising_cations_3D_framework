#!/bin/bash -x
function go {
ni=$(ps aux | grep 'ising' | sed '/grep/d' | wc -l)
while [ $(echo "$ni >= 4" | bc -l) == 1 ] ; do
 sleep 10
 ni=$(ps aux | grep 'ising' | sed '/grep/d' | wc -l)
done
./ising_frameworks < input > ${molar_fraction}.${repetition}.output
}
make install
if [ ! -d CALCS ] ; then mkdir CALCS ; fi
for molar_fraction in $(seq 1 60) ; do
 for repetition in $(seq 1 4) ; do
  echo ${molar_fraction}.${repetition}
  echo "${molar_fraction}" > input
  go
  name=$(grep 'filename:' ${molar_fraction}.${repetition}.output | awk '{print $2}')
  cat ${name} src/oxygen.gin > ${molar_fraction}.${repetition}.gin
  mv ${name} CALCS/.
  mv ${molar_fraction}.${repetition}.output ${molar_fraction}.${repetition}.gin CALCS
 done
done
