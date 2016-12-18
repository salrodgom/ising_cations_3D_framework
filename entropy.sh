#!/bin/bash
for i in $(seq 1 59) ; do
 sed "s/N_SUBSTITUTIONS/$i/g" INSOD.back > INSOD
 ./comsod_mod > output.entropy
 entropy=$(grep 'Maximum entropy for this composition:' output.entropy | head -n1 | awk '{print $6}')
 entropy_sgo=$(grep 'Maximum entropy for this composition:' output.entropy | tail -n1 | awk '{print $9}')
 echo $i $entropy $entropy_sgo
 rm output.entropy INSOD
done
