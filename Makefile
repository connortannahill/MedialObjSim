CC=g++
IDIRS=../src
CFLAGS=-Wall -std=c++11 -I$(IDIRS) -O3
DRIVERS2D = ./TestDrivers/2DDrivers
DRIVERS3D = ./TestDrivers/3DDrivers

SRC2D=$(wildcard ./src/2DSolver/*.cpp) $(wildcard ./src/LASolver/*.cpp) $(wildcard ./src/Utils/*.cpp)
SRC3D=$(wildcard ./src/3DSolver/*.cpp) $(wildcard ./src/LASolver/*.cpp) $(wildcard ./src/Utils/*.cpp)

# Basic drivers
###############

sol2d : main.cpp $(SRC2D) $(LIBS)
	$(CC) $(CFLAGS) $^ -o fluidsolver2d.exe

sol3d : main.cpp $(SRC3D) $(LIBS)
	$(CC) $(CFLAGS) $^ -o fluidsolver3d.exe

# All Test Drivers
#################

all: test_circle_flow_2d

# 2D Test Drivers
#################

test_circle_flow_2d : $(DRIVERS2D)/circle_flow.cpp $(SRC2D) $(LIBS)
	$(CC) $(CFLAGS) $^ -o circle_flow_2d.exe



# 3D Test Drivers
#################


# fluidsolver2d.exe : $(SRC2D) $(LIBS)
# 	$(CC) $(CFLAGS) $^ -o $@
# fluidsolver3d.exe : $(SRC3D) $(LIBS)
# 	$(CC) $(CFLAGS) $^ -o $@

print-% : ;@echo $* = $($*)
