#!/bin/bash
echo "removing old files"
rm -f libHHKinFit.so run

echo "creating shared library"
g++ -fPIC -shared src/*.cpp `root-config --cflags --glibs` -D HHKINFIT_STANDALONE -I ./include -o libHHKinFit.so

echo "creating executable"
g++ example/example.C `root-config --cflags --glibs` -D HHKINFIT_STANDALONE -I ./include -L . -lHHKinFit  -o runHHKinFit
