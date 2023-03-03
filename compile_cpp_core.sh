#!/bin/bash

rm -rf bin
mkdir bin

for i in 256 384 512 768 1024 1280 1536 2048 3072 4096 6144 8192 10240 12288 16384 24576 32768 49152 65536
do
	g++ -fPIC  -O3 -march=native -funroll-loops -std=c++14 -c -Dmaxn=$i coderedlib.cpp -o bin/coderedlib-$i.o # -march=native might not be supported on arm 
	g++ -shared -O3 -march=native -funroll-loops -std=c++14 bin/coderedlib-$i.o -o bin/coderedlib-$i.so
done
