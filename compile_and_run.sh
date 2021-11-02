#!/bin/bash

testName=$1
echo "Test = $testName"

dim=$2
echo "Dim = $dim"

snapshot=$3

if [ snapshot = 1 ]
then
  echo "Taking snapshot"
fi



if [[ $dim = 2 ]]
then
  make
  ./fluidsolver2d.exe $testName 10000 $snapshot
elif [[ $dim = 3 ]]
then
  make sol3d
  ./fluidsolver3d.exe $testName 10000 $snapshot
fi