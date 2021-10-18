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
MAX_STEPS = 100000

FUNS = """
create_test()
grid_scale_test()
obj_scale_test()
full_scale_test()
run_scale_experiment()
create_scale_plot()
run_solver()
exit()
"""

def removePadding(s):
    return ''.join([str.strip()+'\n' for str in s.split('\n')])

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
        $cons_u
        $cons_v
        $g_x
        $g_y
        $Re
        $useEno
        $bcType
        $tEnd"""
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
        $tEnd"""
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
        c $cVal"""
    return OBJ_TEMPLATE

def getObjTemplate3D():
    OBJ_TEMPLATE = """
        $objType
        $objMove
        $cons_u
        $cons_v
        $cons_w
        mass $massVal
        density $densityVal
        r $r
        E $Eval
        eta $etaVal
        a $aVal
        c $cVal"""
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
    template = removePadding(template)
    s = Template(template)
    inFile = s.substitute(in_dict)

    f.write(inFile)

def outputPacking(dim, in_dict, obj_dict, f, lims, eps):
    # Compute the locations
    XA = XB = YA = YB = ZA = ZB = 0
    print('lims = {}'.format(lims))
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
        kMax = int(np.floor((ZB - ZA)/(2*(r + eps))))
    
    numObjs = iMax*jMax*kMax

    if numObjs == 0:
        return

    f.write('\n{}\n'.format(numObjs))

    """ Now generate the centers for each object """
    centerList = []
    if dim == 2:
        centerList = [(XA+tup[0]*2*b+r, YA+tup[1]*2*b+r) for tup in product(range(0, iMax), range(0, jMax))]
    else:
        centerList = [(XA+tup[0]*2*(b)+(r), YA+tup[1]*2*(b)+(r), ZA+tup[2]*2*(b)+(r)) for tup in product(range(0, iMax), range(0, jMax), range(0, kMax))]
    
    print('size of centerList = {}'.format(len(centerList)))
    
    """ For each of these centers, give the output for the file! """

    # fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
    for center in centerList:
        
        create_input_from_dict(obj_dict, (getObjTemplate2D() if dim == 2 else getObjTemplate3D()), f)
        
        # Add a random rotation between 0 and 360
        f.write('deg {}\n'.format(random.randint(0, 360)))

        f.write('cx {0}\n'.format(center[0]))
        f.write('cy {0}\n'.format(center[1]))
        if dim == 3:
            f.write('cz {0}\n'.format(center[2]))
        
        f.write('.\n')

def full_scale_test():
    nVals = [(2**i)*10 for i in range(2, 8)]

    outFileName = input('test name = ')
    dim = int(input(('dim = ')))

    paramList = [s[1:] for s in (getTemplate2D() if dim == 2 else getTemplate3D()).split()]

    paramList.remove('nx')
    paramList.remove('ny')
    if dim == 3:
        paramList.remove('nz')

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

    hx = abs(float(in_dict['xb']) - float(in_dict['xa']))/float(in_dict['nx'])
    hy = abs(float(in_dict['yb']) - float(in_dict['ya']))/float(in_dict['ny'])
    hz = 1
    if dim == 3:
        hz = abs(float(in_dict['zb']) - float(in_dict['za']))/float(in_dict['zy'])

    if dim == 2:
        epsMin = np.linalg.norm([5*hx, 5*hy])
    else:
        epsMin = np.linalg.norm([5*hx, 5*hy, 5*hz])

    xa, xb = float(in_dict['xa']), float(in_dict['xb'])
    ya, yb = float(in_dict['ya']), float(in_dict['yb'])
    za, zb = float(in_dict['za']), float(in_dict['zb'])

    diffs = [xb - xa, yb - ya, zb - za]

    for n in nVals:
        min_edge = min(diffs)

        min_idx = diffs.index(min_edge)

        if min_idx == 0:
            nTups.append((n, int((diffs[1]/min_edge)*n), int((diffs[2]/min_edge)*n)))
        elif min_idx == 1:
            nTups.append((int((diffs[0]/min_edge)*n), n, int((diffs[2]/min_edge)*n)))
        else:
            nTups.append((int((diffs[0]/min_edge)*n), int((diffs[1]/min_edge)*n), n))

    objParamList = [s[1:] for s in (getObjTemplate2D().split() if dim == 2 else getObjTemplate3D().split()) if s[0] == '$']

    obj_dict = {s: input('{} = '.format(s)) for s in objParamList}

    for n in nTups:
        testName = outFileName + str(min(n))

        with open(out_dir_template.format(dim, str(testName)), 'w') as f:
            f.write('A grid-based scaling test for performance\n')
            f.write(testName + '\n\n')

            in_dict['nx'] = str(n[0])
            in_dict['ny'] = str(n[1])

            if dim == 3:
                in_dict['nz'] = str(n[2])

            outputPacking(dim, in_dict, obj_dict, f, lims, epsMin)

def obj_scale_test():
    outFileName = input('test name = ')
    dim = int(input)

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

    objParamList = [s[1:] for s in (getObjTemplate2D().split() if dim == 2 else getObjTemplate3D().split()) if s[0] == '$']
    print(objParamList)

    obj_dict = {s: input('{} = '.format(s)) for s in objParamList}
    del obj_dict['r']

    for r in rVals:
        testName = outFileName + 'Radius{0:.2}'.format(float(r))

        with open(out_dir_template.format(dim, str(testName)), 'w') as f:
            f.write('A fixed-grid scaling test for the number of objects\n')
            f.write(testName + '\n\n')

            obj_dict['r'] = r

            outputPacking(dim, in_dict, obj_dict, f, lims, epsMin)

def grid_scale_test():
    nVals = [(2**i)*10 for i in range(2, 8)]

    outFileName = input('test name = ')
    dim = int(input('dim = '))

    paramList = [s[1:] for s in (getTemplate2D() if dim == 2 else getTemplate3D()).split()]

    paramList.remove('nx')
    paramList.remove('ny')
    if dim == 3:
        paramList.remove('nz')

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
    za = zb = 0
    if dim == 3:
        za, zb = float(in_dict['za']), float(in_dict['zb'])

    diffs = [xb - xa, yb - ya, zb - za]

    for n in nVals:
        min_edge = min(diffs[:dim])

        min_idx = diffs.index(min_edge)

        if dim == 2:
            if min_idx == 0:
                nTups.append((n, int((diffs[1]/min_edge)*n)))
            elif min_idx == 1:
                nTups.append((int((diffs[0]/min_edge)*n), n))
        else:
            if min_idx == 0:
                nTups.append((n, int((diffs[1]/min_edge)*n), int((diffs[2]/min_edge)*n)))
            elif min_idx == 1:
                nTups.append((int((diffs[0]/min_edge)*n), n, int((diffs[2]/min_edge)*n)))
            else:
                nTups.append((int((diffs[0]/min_edge)*n), int((diffs[1]/min_edge)*n), n))

    objParamList = [s[1:] for s in (getObjTemplate2D().split() if dim == 2 else getObjTemplate3D().split()) if s[0] == '$']

    obj_dict = {s: input('{} = '.format(s)) for s in objParamList}

    eps = 0

    for i, n in enumerate(nTups):
        testName = outFileName + str(min(n))

        with open(out_dir_template.format(dim, str(testName)), 'w') as f:
            f.write('A grid-based scaling test for performance\n')
            f.write(testName + '\n\n')

            in_dict['nx'] = str(n[0])
            in_dict['ny'] = str(n[1])
            if dim == 3:
                in_dict['nz'] = str(n[2])

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

    with open(out_dir_template.format(dim, str(testName)), 'w') as f:
        f.write('Many object test generic, can be from many sources\n')
        f.write(testName + '\n')

        paramList = [s[1:] for s in (getTemplate2D() if dim == 2 else getTemplate3D()).split()]

        in_dict = {s: input('{} = '.format(s)) for s in paramList}

        objParamList = [s[1:] for s in (getObjTemplate2D().split() if dim == 2 else getObjTemplate3D().split()) if s[0] == '$']

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
        
        n = None
        if dim == 2:
            n = (int(in_dict['nx']), int(in_dict['ny']))
        else:
            n = (int(in_dict['nx']), int(in_dict['ny']), int(in_dict['nz']) )
        
        print('n = {}'.format(n))
        print('xa = {0} xb = {1}'.format(in_dict['xa'], in_dict['xb']))
        print('ya = {0} yb = {1}'.format(in_dict['ya'], in_dict['yb']))
        print('za = {0} zb = {1}'.format(in_dict['za'], in_dict['zb']))

        nx = n[0]
        ny = n[1]

        hx = abs(float(in_dict['xb']) - float(in_dict['xa']))/nx
        print('h_x = {}'.format(hx))
        hy = abs(float(in_dict['yb']) - float(in_dict['ya']))/ny
        print('h_y = {}'.format(hy))

        hz = 0
        if dim == 3:
            nz = n[2]
            hz = abs(float(in_dict['zb']) - float(in_dict['za']))/nz
            print('h_z = {}'.format(hz))

        # Compute the tolerance for the embedding 5h
        epsMin = 0
        if dim == 2:
            epsMin = 5*np.linalg.norm([hx, hy])
        else:
            epsMin = 5*np.linalg.norm([hx, hy, hz])
        print('epsMin = {}'.format(epsMin))

        eps = float(input('eps = '))
        while eps < epsMin:
            eps = float(input('This packing tolerance {0} will be too low for the chosen grid density, please choose another > {1}: '.format(eps, epsMin)))

        outputPacking(dim, in_dict, obj_dict, f, lims, eps)

def run_scale_experiment():

    testName = input('test name = ')
    dim = int(input('dimension = '))

    # Get all input file names
    inputFiles = [file for file in glob.glob('./TestDrivers/{0}DDrivers/{1}*'.format(str(dim), testName))]
    print(inputFiles)
    inputFiles = [file[file.rfind('/')+1:]  for file in inputFiles]
    print(inputFiles)
    num_list = np.argsort([int(file[len(testName):]) for file in inputFiles])

    inputFiles = list(np.array(inputFiles)[num_list])

    if dim == 2:
        subprocess.run('make')
    else:
        subprocess.run('make sol3d')

    for i in range(len(inputFiles)):
        times = []
        num_runs = 1
        import time

        for run in range(num_runs):
            start = time.time()

            if dim == 2:
                subprocess.run('./fluidsolver2d.exe {0} {1} 1'.format(inputFiles[i], MAX_STEPS).split())
            else:
                subprocess.run('./fluidsolver3d.exe {0} {1} 1'.format(inputFiles[i], MAX_STEPS).split())
            times.append(time.time() - start)
        
        # Dump the data file
        Path("output/ScaleTest/{0}".format(testName)).mkdir(parents=True, exist_ok=True)

        with open('output/ScaleTest/{0}/{1}.out'.format(testName, inputFiles[i]), 'w+') as f:
            f.write('{}\n'.format(num_runs))
            for time in times:
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
        ax.set_xticklabels(ticks)

        plt.xlabel('$N^{}$'.format(dim))
        plt.ylabel('Average CPU Time')
        plt.title('{}'.format(testName))
        plt.legend()

        pName = "output/ScaleTest/{0}/{1}/".format(testName, num_list[i])
        print(pName)
        Path(pName).mkdir(parents=True, exist_ok=True)
        plt.savefig('{0}ParTest{1}{2}.png'.format(pName, testName, num_list[i]))

def run_solver():
    test_name = input('test name = ')
    dim = int(input('dim = '))

    if dim == 2:
        subprocess.run('make')
    else:
        subprocess.run('make sol3d')

    start = time.time()

    if dim == 2:
        subprocess.run('./fluidsolver2d.exe {0} {1} 1'.format(inputFiles[i], MAX_STEPS).split())
    else:
        subprocess.run('./fluidsolver3d.exe {0} {1} 1'.format(inputFiles[i], MAX_STEPS).split())
    times = (time.time() - start)

    print('took {} secs'.format(times))

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