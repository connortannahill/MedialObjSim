# BloodSim

Requires:
- Eigen [https://eigen.tuxfamily.org/dox/GettingStarted.html](https://eigen.tuxfamily.org/dox/GettingStarted.html)
- nanoflann.hpp [https://github.com/jlblancoc/nanoflann/blob/master/include/nanoflann.hpp](https://github.com/jlblancoc/nanoflann/blob/master/include/nanoflann.hpp)
- Python libraries: matplotlib, numpy, itertools, seaborn, plotly (for 3D)

To build, run `make`

Before running, make sure to create an `output` folder in repo if it doesn't already exist.

To run fluid simulation after building, run generated exe file. Example: 
```sh
./fluidsolver2d.exe <max_steps>
```

To plot:
```sh
python plot.py <test_name> <mode> <num_objects>
```
where `test_name` is from `output/<test_name>` file generated from running C++ code
