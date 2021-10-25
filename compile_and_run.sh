#!/bin/bash
make sol3d
./fluidsolver3d.exe Q2ManyObject3D 10000
python plot3D.py Q2ManyObject3D 1   