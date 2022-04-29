# MedialObjSim

Requires:
- Eigen [https://eigen.tuxfamily.org/dox/GettingStarted.html](https://eigen.tuxfamily.org/dox/GettingStarted.html)
- nanoflann.hpp [https://github.com/jlblancoc/nanoflann/blob/master/include/nanoflann.hpp](https://github.com/jlblancoc/nanoflann/blob/master/include/nanoflann.hpp)
- Python libraries: matplotlib, numpy, itertools, seaborn, plotly (for 3D)

To build, run `make`

Before running, make sure to create an `output` folder in repo if it doesn't already exist.

To run fluid simulation (in two dimensions) after building, run generated exe file. Example: 
```sh
./fluidsolver2d.exe <max_steps>
```

To plot:
```sh
python plot.py <test_name> <mode> <num_objects>
```
where `test_name` is from `output/<test_name>` file generated from running C++ code


### Test Driver Input Files

In folder `TestDrivers/2DDrivers/`, see `template.txt` for the format of the input file and `test_1.txt` for basic usage.

To run, change Makefile `sol2d` to `main_from_txt.cpp`. Example command line call:

```sh
./fluidsolver2d.exe TestDrivers/2DDrivers/test_1.txt 100
```

