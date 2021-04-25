#!/bin/bash

mkdir -p auxiliar
cd auxiliar
wget https://github.com/Jaimenms/dasslc2py/archive/refs/heads/master.zip
unzip master.zip
rm master.zip
cd dasslc2py-master
python setup.py install
cd ../..
rm -rf auxiliar/
