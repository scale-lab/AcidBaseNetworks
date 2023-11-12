#!/bin/bash
set -e
source venv/bin/activate
export PYTHONPATH=$(pwd):$(pwd)/chemcpupy

for sz in 8 12 16 28 
do
  python simulation/acidbase_writer_acc.py bin ${sz} 3 > acc_bin_3class_${sz}.txt &
done
wait