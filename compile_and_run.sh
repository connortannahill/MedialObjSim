#!/bin/bash
make
./fluidsolver3d.exe Q2ManyObject3D 0
python plot3D.py Q2ManyObject3D 1   