#!/usr/bin/env bash
for i in `seq 1 100';
do
    bsub < jobscript
done
