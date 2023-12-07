#!/usr/bin/env bash

# Step 1: Send file to remote computer
scp -r intopt.c ka3617pe-s@power.cs.lth.se:

# Step 2: SSH into remote, compile and get size
ssh -q -T ka3617pe-s@power.cs.lth.se 'mv intopt.c intopt/; cd intopt/; make;'
#size --common poly.o' | tail -1 | awk '{ print $4; }'
