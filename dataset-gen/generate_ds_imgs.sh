#!/bin/bash
set -e

source venv/bin/activate
echo "=================================="
python mnist_generate_echo_ds.py bin 8 2 true &
echo "=================================="
python mnist_generate_echo_ds.py bin 12 2 true &
echo "=================================="
python mnist_generate_echo_ds.py bin 16 2 true &
echo "=================================="
python mnist_generate_echo_ds.py bin 28 2 true &
echo "=================================="

echo "=================================="
python mnist_generate_echo_ds.py bin 8 3 true &
echo "=================================="
python mnist_generate_echo_ds.py bin 12 3 true &
echo "=================================="
python mnist_generate_echo_ds.py bin 16 3 true &
echo "=================================="
python mnist_generate_echo_ds.py bin 28 3 true &
echo "=================================="


echo "=================================="
python mnist_generate_echo_ds.py 3bit 8 2 true &
echo "=================================="
python mnist_generate_echo_ds.py 3bit 12 2 true &
echo "=================================="
python mnist_generate_echo_ds.py 3bit 16 2 true &
echo "=================================="
python mnist_generate_echo_ds.py 3bit 28 2 true &
echo "=================================="
wait
echo "DONE!"