#!/bin/bash
mkdir Minuit2-5.34.14
curl https://scisoft.fnal.gov/scisoft/packages/minuit2/v5_34_14/minuit2-5.34.14-source.tar.bz2|tar -C Minuit2-5.34.14 xjvf -
#curl http://project-mathlibs.web.cern.ch/project-mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz|tar xzvf -
cd Minuit2-5.34.14
mkdir build
cd build
../configure --prefix=$HOME
make -j 2
make install
cd ../../
rm -fr Minuit2-5.34.14

echo 'Add: export CPATH=$CPATH:$HOME/include to the .bashrc'
