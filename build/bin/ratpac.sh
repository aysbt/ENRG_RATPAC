#!/bin/bash

export RATROOT=/usr/local
export PATH=$RATROOT/bin:$PATH
export LD_LIBRARY_PATH=$RATROOT/lib:$LD_LIBRARY_PATH
export RATSHARE=/usr/local/share/RAT
export GLG4DATA=$RATSHARE/data
export PYTHONPATH=$RATSHARE/python:$PYTHONPATH
export ROOT_INCLUDE_PATH=$RATROOT/include
# For Mac OS X
export DYLD_LIBRARY_PATH=$RATROOT/lib:$DYLD_LIBRARY_PATH
