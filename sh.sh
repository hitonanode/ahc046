#!/bin/bash

set -eux

CASE_ID="$1"

# g++ -std=c++20 -O2 -DHITONANODE_LOCAL main.cpp
g++ -O2 -std=c++20 -DHITONANODE_LOCAL main.cpp

cat in/$CASE_ID.txt | ./a.out | pbcopy
