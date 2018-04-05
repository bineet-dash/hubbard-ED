#!/bin/bash
for i in {0..48}
do
 ./parallel_Ht.out $(($1*100)) $(($i*100)) 100 8 0
done