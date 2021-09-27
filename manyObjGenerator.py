"""
Script for generating some bulk testing scripts
"""

import numpy as np
from string import Template
import random
import subprocess
from matplotlib.pyplot import cm
from itertools import product
import matplotlib.pyplot as plt
import glob
import time
from pathlib import Path
import json
import re

"""
TEMPLATES AND MISC
"""

FUNS = """
create_test()
grid_scale_test_2d()
obj_scale_test_2d()
full_scale_test_2d()
run_mesh_scale_experiment()
create_scale_plot()
exit()
"""

def exit():
    import sys
    sys.exit()

out_dir_template = './TestDrivers/{0}DDrivers/{1}'

FUNS_LIST = FUNS.split()

def getTemplate2D():
    TEMPLATE_2D = """
        $xa
        $xb
        $ya
        $yb
        $nx
        $ny
        $g_x
        $g_y
        $Re
        $useEno
        $bcType
        $tEnd
        """
    return TEMPLATE_2D

def getTemplate3D():
    TEMPLATE_3D = """
        $xa
        $xb
        $ya
        $yb
        $za
        $zb
        $nx
        $ny
        $nz
        $cons_u
        $cons_v
        $cons_w
        $g_x
        $g_y
        $g_z
        $Re
        $useEno
        $bcType
        $tEnd
        """
    return TEMPLATE_3D

def getObjTemplate2D():
    OBJ_TEMPLATE = """
        $objType
        $objMove
        $cons_u
        $cons_v
        mass $massVal
        density $densityVal
        r $r
        E $Eval
        eta $etaVal
        a $aVal
        c $cVal
        """
    return OBJ_TEMPLATE

def plot_scale_experiment(testName, nums, times, ax, color, label):
    
    import statsmodels.stats.api as sms

    # Build each of the confidence intervals
    mean = []
    low_int = []
    high_int = []
    mean_max = 0
    low_max = 0
    high_max = 0
    i = 0

    denom = 1

    for num, time_list in zip(num, times):
        time_list = [time/denom for time in time_list]
        interval = sms.DescrStatsW(time_list).tconfint_mean()
        mean.append(np.mean(time_list))
        low_int.append(interval[0])
        high_int.append(interval[1])
        if i == 0:
            mean_max = mean[0]
            low_max = low_int[0]
            high_max = high_int[0]
        i += 1
    
    # ax.plot(np.arange(len(nums)), mean / mean_max, label=label, color=color)
    # ax.fill_between(np.arange(len(nums)), low_int / low_max, high_int / high_max, color=color, alpha=.1)
    ax.plot(np.arange(len(nums)), mean, label=label, color=color)
    ax.fill_between(np.arange(len(nums)), low_int, high_int, color=color, alpha=.1)


def create_input_from_dict(in_dict, template, f):
    s = Template(template)
    inFile = s.substitute(in_dict)

    f.write(inFile)

def outputPacking(dim, in_dict, obj_dict, f, lims, eps):
    # Compute the locations
    XA = XB = YA = YB = ZA = ZB = 0
    print(lims)
    if dim == 2:
        XA, XB, YA, YB = lims
    else:
        XA, XB, YA, YB, ZA, ZB = lims

    
    r = float(obj_dict['r'])
    b = r+eps

    create_input_from_dict(in_dict, (getTemplate2D() if dim == 2 else getTemplate3D()), f)
    f.write('\n')
    
    # Compute the maximum
    iMax = int(np.floor((XB - XA)/(2*(r + eps))))
    jMax = int(np.floor((YB - YA)/(2*(r + eps))))
    kMax = 1
    if dim == 3:
        kMax = np.floor((ZB - ZA)/(2*(r + eps)))
    
    numObjs = iMax*jMax*kMax

    if numObjs == 0:
        return

    f.write('\n{}\n'.format(numObjs))

    """ Now generate the centers for each object """
    centerList = []
    if dim == 2:
        centerList = [(XA+tup[0]*2*b+r, YA+tup[1]*2*b+r) for tup in product(range(0, iMax), range(0, jMax))]
    else:
        centerList = [(XA+tup[0]*2*(b)+(r), YA+tup[1]*2*(b)+(r), ZA+tup[2]*2*(b)+(r)) for tup in product(range(0, iMax), range(0, jMax), range(0, kMax+1))]
    
    """ For each of these centers, give the output for the file! """

    # fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
    for center in centerList:
        
        create_input_from_dict(obj_dict, getObjTemplate2D(), f)
        
        # Add a random rotation between 0 and 360
        f.write('deg {}\n'.format(random.randint(0, 360)))

        f.write('cx {0}\n'.format(center[0]))
        f.write('cy {0}\n'.format(center[1]))
        if dim == 3:
            f.write('cz = {0}\n'.format(center[2]))
        
        f.write('.\n')

def full_scale_test_2d():
    nVals = [(2**i)*10 for i in range(6)]

    outFileName = input('test name = ')
    dim = 2

    paramList = [s[1:] for s in (getTemplate2D() if dim == 2 else getTemplate3D()).split()]

    paramList.remove('nx')
    paramList.remove('ny')

    in_dict = {s: input('{} = '.format(s)) for s in paramList}

    nTups = []

    XA = float(input('XA = '))
    XB = float(input('XB = '))

    YA = float(input('YA = '))
    YB = float(input('YB = '))
    if dim == 3:
        ZA = float(input('ZA = '))
        ZB = float(input('ZB = '))
        lims = [XA, XB, YA, YB, ZA, ZB]
    else:
        lims = [XA, XB, YA, YB]

    hx = abs(in_dict['xb'] - in_dict['xa'])/in_dict['nx']
    hy = abs(in_dict['yb'] - in_dict['ya'])/in_dict['ny']

    epsMin = np.linalg.norm([5*hx, 5*hy])

    xa, xb = float(in_dict['xa']), float(in_dict['xb'])
    ya, yb = float(in_dict['ya']), float(in_dict['yb'])
    x_diff = xb - xa
    y_diff = yb - ya
    for n in nVals:
        if (min(x_diff, y_diff) == x_diff):
            scal = y_diff / x_diff
            nTups.append((n, int(scal*n)))
        else:
            scal = x_diff / y_diff
            nTups.append((int(scal*n), n))

    objParamList = [s[1:] for s in getObjTemplate2D().split() if s[0] == '$']

    obj_dict = {s: input('{} = '.format(s)) for s in objParamList}

    for n in nTups:
        testName = outFileName + str(min(n))

        with open(out_dir_template.format(dim, str(testName)+'.txt'), 'w') as f:
            f.write('A grid-based scaling test for performance\n')
            f.write(testName + '\n\n')

            in_dict['nx'] = str(n[0])
            in_dict['ny'] = str(n[1])

            outputPacking(dim, in_dict, obj_dict, f, lims, epsMin)

def obj_scale_test_2d():
    outFileName = input('test name = ')
    dim = 2

    paramList = [s[1:] for s in (getTemplate2D() if dim == 2 else getTemplate3D()).split()]
    print(paramList)

    in_dict = {s: input('{} = '.format(s)) for s in paramList}

    XA = float(input('XA = '))
    XB = float(input('XB = '))

    YA = float(input('YA = '))
    YB = float(input('YB = '))
    if dim == 3:
        ZA = float(input('ZA = '))
        ZB = float(input('ZB = '))
        lims = [XA, XB, YA, YB, ZA, ZB]
    else:
        lims = [XA, XB, YA, YB]

    hx = abs(float(in_dict['xb']) - float(in_dict['xa']))/int(in_dict['nx'])
    hy = abs(float(in_dict['yb']) - float(in_dict['ya']))/int(in_dict['ny'])

    epsMin = np.linalg.norm([5*hx, 5*hy])

    # Get the initial r value
    rBase = min((XA + XB)/2.0, (YA + YB)/2.0)
    
    rVals = [rBase / (2**i) - epsMin/2.0 for i in range(6)]

    objParamList = [s[1:] for s in getObjTemplate2D().split() if s[0] == '$']
    print(objParamList)

    obj_dict = {s: input('{} = '.format(s)) for s in objParamList}
    del obj_dict['r']

    for r in rVals:
        testName = outFileName + 'Radius{0:.2}'.format(float(r))

        with open(out_dir_template.format(dim, str(testName)+'.txt'), 'w') as f:
            f.write('A fixed-grid scaling test for the number of objects\n')
            f.write(testName + '\n\n')

            obj_dict['r'] = r

            outputPacking(dim, in_dict, obj_dict, f, lims, epsMin)

def grid_scale_test_2d():
    nVals = [(2**i)*10 for i in range(2, 8)]

    outFileName = input('test name = ')
    dim = 2

    paramList = [s[1:] for s in (getTemplate2D() if dim == 2 else getTemplate3D()).split()]

    paramList.remove('nx')
    paramList.remove('ny')

    in_dict = {s: input('{} = '.format(s)) for s in paramList}

    nTups = []

    XA = float(input('XA = '))
    XB = float(input('XB = '))

    YA = float(input('YA = '))
    YB = float(input('YB = '))
    if dim == 3:
        ZA = float(input('ZA = '))
        ZB = float(input('ZB = '))
        lims = [XA, XB, YA, YB, ZA, ZB]
    else:
        lims = [XA, XB, YA, YB]


    xa, xb = float(in_dict['xa']), float(in_dict['xb'])
    ya, yb = float(in_dict['ya']), float(in_dict['yb'])
    x_diff = xb - xa
    y_diff = yb - ya
    for n in nVals:
        if (min(x_diff, y_diff) == x_diff):
            scal = y_diff / x_diff
            nTups.append((n, int(scal*n)))
        else:
            scal = x_diff / y_diff
            nTups.append((int(scal*n), n))

    objParamList = [s[1:] for s in getObjTemplate2D().split() if s[0] == '$']

    obj_dict = {s: input('{} = '.format(s)) for s in objParamList}

    eps = 0

    for i, n in enumerate(nTups):
        testName = outFileName + str(min(n))

        with open(out_dir_template.format(dim, str(testName)+'.txt'), 'w') as f:
            f.write('A grid-based scaling test for performance\n')
            f.write(testName + '\n\n')

            in_dict['nx'] = str(n[0])
            in_dict['ny'] = str(n[1])

            if i == 0:
                nx = n[0]
                ny = n[1]

                hx = abs(float(in_dict['xb']) - float(in_dict['xa']))/nx
                hy = abs(float(in_dict['yb']) - float(in_dict['ya']))/ny

                hz = 0
                if dim == 3:
                    nz = n[2]
                    hz = abs(float(in_dict['zb']) - float(in_dict['za']))/nz

                # Compute the tolerance for the embedding 5h
                epsMin = 0
                if dim == 2:
                    epsMin = np.linalg.norm([5*hx, 5*hy])
                else:
                    epsMin = np.linalg.norm([5*hx, 5*hy, 5*hz])
                print('epsMin = {}'.format(epsMin))

                eps = float(input('eps = '))
                while eps < epsMin:
                    eps = float(input('This packing tolerance {0} will be too low for the chosen grid density, please choose another > {1}: '.format(eps, epsMin)))


            outputPacking(dim, in_dict, obj_dict, f, lims, eps)

def create_test():
    testName = input('Enter the test name: ')
    dim = int(input('Dimension = '))

    objs = []

    with open(out_dir_template.format(dim, str(testName)+'.txt'), 'w') as f:
        f.write('Many object test generic, can be from many sources\n')
        f.write(testName + '\n')

        paramList = [s[1:] for s in (getTemplate2D() if dim == 2 else getTemplate3D()).split()]

        in_dict = {s: input('{} = '.format(s)) for s in paramList}

        objParamList = [s[1:] for s in getObjTemplate2D().split() if s[0] == '$']

        obj_dict = {s: input('{} = '.format(s)) for s in objParamList}

        XA = float(input('XA = '))
        XB = float(input('XB = '))

        YA = float(input('YA = '))
        YB = float(input('YB = '))
        if dim == 3:
            ZA = float(input('ZA = '))
            ZB = float(input('ZB = '))
            lims = [XA, XB, YA, YB, ZA, ZB]
        else:
            lims = [XA, XB, YA, YB]

        outputPacking(dim, in_dict, obj_dict, f, lims)

def run_mesh_scale_experiment():

    testName = input('test name = ')
    dim = int(input('dimension = '))

    # Get all input file names
    inputFiles = [file for file in glob.glob('./TestDrivers/{0}DDrivers/{1}*'.format(str(dim), testName))]
    inputFiles_dirgone = [file for file[file.rfind('/')+1:] in inputFiles]
    num_list = np.argsort([int(file[len(testName):]) for file in inputFiles_dirgone])

    inputFiles = list(np.array(inputFiles)[num_list])
    inputFiles_dirgone = list(np.array(inputFiles_dirgone)[num_list])

    subprocess.run('make')

    for i in range(len(inputFiles)):
        times = []
        num_runs = 10

        for run in num_runs:
            start = time.time()

            subprocess.run('./mesh.exe {0}'.format(inputFiles[i]).split())
            times.append(time.time() - start)
        
        # Dump the data file
        Path("output/ScaleTest/{1}".format(testName)).mkdir(parents=True, exist_ok=True)

        with open('output/ScaleTest/{0}/{1}.out'.format(testName, inputFiles_dirgone), 'w+') as f:
            f.write('{}\n'.format(num_runs))
            for time in times():
                f.write(str(time) + ' ')

def create_scale_plot():
    testName = input('test name = ')
    dim = int(input('dimension = '))

    # Get all input file names
    dataFiles= [file for file in glob.glob('output/ScaleTest/{0}/*'.format(str(dim), testName))]
    dataFiles_dirgone = [file for file[file.rfind('/')+1:] in dataFiles]
    num_list = np.argsort([int(file[len(testName):]) for file in dataFiles_dirgone])
    num_list = [int(file[len(testName):]) for file in dataFiles]

    sort_list = np.argsort(num_list)
    num_list = list(np.array(num_list)[sort_list])
    dataFiles = list(np.array(dataFiles)[sort_list])

    color=cm.rainbow(np.linspace(0,1,len(dataFiles)))
    fig, ax = plt.subplots()

    for i, dataFile in enumerate(dataFiles):
        times = {}
        print(dataFile)
        num_runs = 0
        times = []
        with open('output/ScaleTest/{0}/{1}'.format(testName, dataFile)) as f:
            input = [float(num) for num in f.readlines()]
            num_runs = int(input[0])
            times = input[1:]
        
        fig, ax = plt.subplots()
        
        # Dump the data file
        label = str(num_list[i])
        plot_scale_experiment(dataFiles[i], num_list, times, ax, color[i], label)

        # Make modified ticks
        ticks = ['${}$'.format(num) for num in num_list]
        ticks.insert(0, '$1$')
        ax.set_xticklabels(ticks)

        plt.xlabel('$N^{}$'.format(dim))
        plt.ylabel('Average CPU Time')
        plt.title('{}'.format(testName))
        plt.legend()

        pName = "output/ScaleTest/{0}/{1}/".format(testName, num_list[i])
        print(pName)
        Path(pName).mkdir(parents=True, exist_ok=True)
        plt.savefig('{0}ParTest{1}{2}.png'.format(pName, testName, num_list[i]))
        

while(True):
    print('====================================')
    print("Functions available:")
    print(FUNS)
    print('====================================')

    mode = input('call function() = ')
    if mode in FUNS_LIST:
        eval(mode)
    else:
        print('Function {} does not exist!'.format(mode))