"""
Script for generating some bulk testing scripts
"""

import numpy as np
import random
from itertools import product

testName = input('Enter the test name: ')
dim = int(input('Dimension = '))

objs = []

with open('TestDrivers/{0}DDrivers/{1}.txt'.format(dim, testName), 'w') as f:
    f.write('Many object test generic, can be from many sources\n')
    f.write(testName + '\n\n')

    print('Getting the parameters for the simulation')
    line = input('xa <space> xb = ')
    f.write('{0} {1}\n'.format(*[float(i) for i in line.split()]))

    line = input('ya <space> yb = ')
    f.write('{} {}\n'.format(*[float(i) for i in line.split()]))

    if dim == 3:
        line = input('za <space> zb = ')
        f.write('{} {}\n'.format(*[float(i) for i in line.split()]))


    nx, ny, nz = 0, 0, 0
    if dim == 2:
        nx, ny = [int(i) for i in input('nx <space> ny = ').split()]
        f.write('{} {}\n'.format(nx, ny))
    else :
        nx, ny, nz = [int(i) for i in input('nx <space> ny <space> nz = ').split()]
        f.write('{} {} {}\n'.format(nx, ny, nz))

    u_cons, v_cons, w_cons = 0, 0, 0
    if dim == 2:
        u_cons, v_cons = [float(i) for i in input('cons_u <space> cons_v = ').split()]
        f.write('{} {}\n'.format(u_cons, v_cons))
    else:
        u_cons, v_cons, w_cons = [int(i) for i in input('cons_u <space> cons_v <space> cons_w = ').split()]
        f.write('{} {} {}\n'.format(u_cons, v_cons, w_cons))

    g_x, g_y, g_z = 0, 0, 0
    if dim == 2:
        g_x, g_y = [float(i) for i in input('g_x <space> g_y = ').split()]
        f.write('{} {}\n'.format(g_x, g_y))
    else:
        g_x, g_y, g_z = [int(i) for i in input('g_x <space> g_y <space> g_z = ').split()]
        f.write('{} {} {}\n'.format(g_x, g_y, g_z))


    f.write('{}\n'.format(float(input('Re = '))))
    f.write('{}\n'.format(int(bool(input('useEno = ')))))
    f.write('{}\n'.format(input('BC type = ')))
    f.write('{}\n'.format(float(input('tEnd = '))))
    f.write('\n')

    print('The subset of the domain where we place objects')
    XA, XB = [float(x) for x in input('XA <space> XB = ').split()]
    YA, YB = [float(x) for x in input('YA <space> YB = ').split()]

    ZA, ZB = 0.0, 0.0
    if dim == 3:
        ZA, ZB = [float(x) for x in input('ZA <space> ZB = ')]
    
    obj_type = input('Type of object = ')

    r = float(input('Radius of objects = '))
    eps = float(input('Packing tolerance (how far apart at minimum) = '))

    hx = abs(XB - XA)/nx
    hy = abs(YB - YA)/ny
    hz = 0
    if dim == 3:
        hz = abs(ZB - ZA)/nz

    # Compute the tolerance for the embedding 5h
    # eps = 0
    if dim == 2:
        epsMin = np.linalg.norm([5*hx, 5*hy])
    else:
        epsMin = np.linalg.norm([5*hx, 5*hy, 5*hz])
    
    while eps < epsMin:
        eps = float(input('This packing tolerance {0} will be too low for the chosen grid density, please choose another > {1}: '.format(eps, epsMin)))
    

    b = r+eps
    
    # Compute the maximum
    iMax = int(np.floor((XB - XA)/(2*(r + eps))))
    jMax = int(np.floor((YB - YA)/(2*(r + eps))))
    kMax = 1
    if dim == 3:
        kMax = np.floor((ZB - ZA)/(2*(r + eps)))
    
    numObjs = iMax*jMax*kMax
    f.write('{}\n'.format(numObjs))
    # f.write('\n')

    print('The parameter list (rotations are generated randomly) and added to all objects, enter "None -1" to stop')
    paramList = {}
    name, val = input('[paramName] [paramVal] = ').split()
    while (val != '-1'):
        print('{0} {1}'.format(name, val))
        paramList[name] = float(val)
        name, val = input('[paramName] [paramVal] = ').split()
    
    print(paramList)

    """ Now generate the centers for each object """
    centerList = []
    if dim == 2:
        centerList = [(XA+tup[0]*2*b+r, YA+tup[1]*2*b+r) for tup in product(range(0, iMax), range(0, jMax))]
    else:
        centerList = [(XA+tup[0]*2*(b)+(r), YA+tup[1]*2*(b)+(r), ZA+tup[2]*2*(b)+(r)) for tup in product(range(0, iMax), range(0, jMax), range(0, kMax+1))]
    
    """ For each of these centers, give the output for the file! """

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
    for center in centerList:
        f.write('\n{}\n'.format(obj_type))
        f.write('1\n')
        f.write('{}\n'.format(u_cons))
        f.write('{}\n'.format(v_cons))
        if dim == 3:
            f.write('{}\n'.format(w_cons))
        f.write('cx {}\n'.format(center[0]))
        f.write('cy {}\n'.format(center[1]))
        f.write('r {}\n'.format(r))
        if dim == 3:
            f.write('cz {}\n'.format(center[2]))
        
        for key, val in paramList.items():
            print('{0} {1}'.format(key, val))
            f.write('{0} {1}\n'.format(key, val))
        
        # Add a random rotation between 0 and 360
        f.write('deg {}\n'.format(random.randint(0, 360)))
        
        f.write('.\n')

        circle = plt.Circle(center, r)
        ax.add_patch(circle)
    
    plt.xlim((0, 5))
    plt.ylim((0, 5))


    # (or if you have an existing figure)
    plt.show()
