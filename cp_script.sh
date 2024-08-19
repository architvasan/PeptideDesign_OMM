#!/bin/bash
for i in $(seq 1 4); do
    let rep1_cp=1+${i}*4
    cp -r replica1 replica$rep1_cp 
    let rep2_cp=2+${i}*4
    cp -r replica2 replica$rep2_cp 
    let rep3_cp=3+${i}*4
    cp -r replica3 replica$rep3_cp 
    let rep4_cp=4+${i}*4
    cp -r replica4 replica$rep4_cp 
done
