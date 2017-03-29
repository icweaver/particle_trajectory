#!/bin/bash
#for i in `seq 5.0 0.1 7.1`;
#do
#  sed -i.bak "116s/.*/            v_frac = $i/" sho.f90
#  make
#done
sed -i.bak "116s/.*/            v_frac = 1/" sho.f90
make
