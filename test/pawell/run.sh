#!/bin/bash

for i in 1 2
do
./ram < testCS2.re$i >& testcs2out$i
cp testCS2.spec testCS2.spec$i
done