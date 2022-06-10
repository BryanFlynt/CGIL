#!/bin/bash

NUM_SRC= 10
NUM_TRG= 10

./make_random_sources 2 ${NUM_SRC} 2 > source.txt
./make_random_targets 2 ${NUM_TRG} > target.txt
./cgil2
