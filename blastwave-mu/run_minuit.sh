#!/usr/bin/env bash
export LD_LIBRARY_PATH=/snap/root-framework/current/usr/local/lib:$LD_LIBRARY_PATH
cd "$HOME/Documents/blastwave-mu"
./build/bwmu "$@"
