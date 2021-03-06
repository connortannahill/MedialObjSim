CC=g++
IDIRS=-I../src -I./lib/eigen/ -I./lib/nanoflann/include/
CFLAGS=-Wall -std=c++17 $(IDIRS) -O3
DRIVERS2D = ./TestDrivers/2DDrivers
DRIVERS3D = ./TestDrivers/3DDrivers

SRC2D=$(wildcard ./src/2DSolver/*.cpp) $(wildcard ./src/LASolver/*.cpp) $(wildcard ./src/Utils/*.cpp)
SRC3D=$(wildcard ./src/3DSolver/*.cpp) $(wildcard ./src/LASolver/*.cpp) $(wildcard ./src/Utils/*.cpp)

# Basic drivers
###############

# sol2d : main.cpp $(SRC2D) $(LIBS)
sol2d : main_from_txt.cpp $(SRC2D) $(LIBS)
	$(CC) $(CFLAGS) $^ -o fluidsolver2d.exe

sol3d : main_from_txt_3d.cpp $(SRC3D) $(LIBS)
	$(CC) $(CFLAGS) $^ -o fluidsolver3d.exe

# All Test Drivers
#################

all: test_circle_flow_2d test_circle_flow_display_steps_2d

# 2D Test Drivers
#################

test_circle_flow_2d : $(DRIVERS2D)/circle_flow.cpp $(SRC2D) $(LIBS)
	$(CC) $(CFLAGS) $^ -o circle_flow_2d.exe

test_circle_flow_display_steps_2d : $(DRIVERS2D)/circle_flow_display_steps.cpp $(SRC2D) $(LIBS)
	$(CC) $(CFLAGS) $^ -o circle_flow_display_steps.exe

col_test_2d : $(DRIVERS2D)/col_test.cpp $(SRC2D) $(LIBS)
	$(CC) $(CFLAGS) $^ -o col_test_2d.exe



# 3D Test Drivers
#################


# fluidsolver2d.exe : $(SRC2D) $(LIBS)
# 	$(CC) $(CFLAGS) $^ -o $@
# fluidsolver3d.exe : $(SRC3D) $(LIBS)
# 	$(CC) $(CFLAGS) $^ -o $@

print-% : ;@echo $* = $($*)
