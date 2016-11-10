#!/bin/bash

wget http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz -O -|tar xzvf -
cd Minuit2-5.34.14
mkdir build
cd build
../configure --prefix=$HOME
make -j 2
make install
cd ../../
rm -fr Minuit2-5.34.14

echo 'Add: export CPATH=$CPATH:$HOME/include to the .bashrc'
