#!/bin/bash

for f in $(grep -l "CANCELLED AT" slurm-*); do
    if [[ -f restart_1/$f ]]; then
        if [[ -z $(grep "CANCELLED AT" restart_1/$f) ]]; then
            cp $f restart_2
        fi
    else
        cp $f restart_2
    fi
done

