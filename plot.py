import matplotlib.pyplot as plt
import numpy as np
import itertools
# import seaborn as sb
import sys
import os, glob

if len(sys.argv) == 1:
    print('Must provide a plotting mode!')
    sys.exit()

testNum = -1
if (len(sys.argv) == 4):
    testNum = int(sys.argv[3])

test = sys.argv[1]
mode = int(sys.argv[2])
f_name = None
if testNum == -1:
    f_name = './output/{0}/'.format(test)
else:
    f_name = './output/{0}/{1}/'.format(test, testNum)

# numObj = 1
# if (len(sys.argv) == 4):
#     numObj = int(sys.argv[3])



# Get the number of objects
# cur_dir = os.getcwd()
# os.chdir(f_name)
print('looking for numObj in {}'.format(f_name+'MSSEdges*'))
# print('glob = {}'.format(glob.glob(f_name+'MSSEdges*')))
numObj = len(glob.glob(f_name+'MSSEdges*'))
# os.chdir(cur_dir)


"""
For plotting the Pool isocontour
"""

if mode == -2:
    str_base = ''
    from random import randint
    with open('./TestDrivers/2DDrivers/{}'.format(test)) as f:
        for ln in f.readlines():
            str_temp = ln.strip()
            print(str_temp)
            if str_temp.startswith('deg'):
                str_temp = str_temp[:str_temp.rfind('g')+2] + ('90' if randint(0, 1) == 0 else '0')
                # str_temp = str_temp[:str_temp.rfind('g')+2] + str(randint(1, 4) * 45)
                # val = randint(0, 360) 
                # print(val)
                # str_temp = "deg " + str(val)

            str_base += str_temp + '\n'

    with open('./TestDrivers/2DDrivers/Q3ManyObject', 'w') as f:
        f.write(str_base)

# assert(False)

if mode == -1:

    t = np.arange(0.0, 1.0 + 0.01, 0.01)
    s = np.cos(2*2*np.pi*t)
    plt.plot(t, s, '-', lw=2)

    plt.xlabel('time (s)')
    plt.ylabel('voltage (mV)')
    plt.title('About as simple as it gets, folks')
    plt.grid(True)

    plt.axes().set_aspect('equal', 'datalim')


    # plt.show()

    # f_name += 'colOut'

    # out = np.genfromtxt(f_name, delimiter=',')
    # x = out[:,0]
    # y = out[:,1]

    # plt.scatter(x, y)
    # plt.show()


if mode == 1:
    f_name += 'out'
    nskip = 1
    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0][::nskip]
    y = out[:,1][::nskip]
    u = out[:,2][::nskip]
    v = out[:,3][::nskip]
    p = out[:,4][::nskip]

    # Create a quiver plot for this velocity field
    fig, ax1 = plt.subplots()
    q = ax1.quiver(x, y, u, v)
    plt.show()
elif mode == 2:
    f_name += 'poolOut'
    out = np.genfromtxt(f_name, delimiter=',')

    x = out[:,0]
    y = out[:,1]
    phi = out[:,2]

    n = int(np.sqrt(phi.size))

    fig, ax = plt.subplots(1)
    img = ax.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)), levels=[0])
    plt.show()
elif mode == 3:
    # Plot the pool velocity field
    f_name += 'poolVel'
    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]
    u = out[:,2]
    v = out[:,3]

    fig, ax1 = plt.subplots()
    q = ax1.quiver(x, y, u, v)
    plt.show()
elif mode == 4:
    f_name += 'poolOut'
    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]
    phi = out[:,2]

    n = int(np.sqrt(phi.size))
    fig, ax = plt.subplots(1)
    img = ax.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)))
    CB = fig.colorbar(img, shrink=0.8, extend='both')
    plt.show()
elif mode == 5:
    # Plot the pool velocity field
    f_name += 'poolVel'
    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]
    u = out[:,2]
    v = out[:,3]

    fig, ax1 = plt.subplots()
    q = ax1.quiver(x, y, u, v)

    f_name += 'poolOut'
    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]
    phi = out[:,2]

    n = int(np.sqrt(phi.size))
    q = ax1.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)))
    plt.show()
elif mode == 6:
    f_name += 'MSSEdges0'

    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]

    for i in range(0, x.size, 2):
        plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=0.5)
    # plt.xlim((0, 1))
    # plt.ylim((0, 1))
    plt.show()
elif mode == 7:
    f_temp = f_name
    for i in range(numObj):
        f_name += 'MSSVels{0}'.format(i)

        out = np.genfromtxt(f_name, delimiter=',')
        x = out[:,0]
        y = out[:,1]
        u = out[:,2]
        v = out[:,3]

        fig, ax1 = plt.subplots()
        q = ax1.quiver(x, y, u, v)

        f_name = f_temp


    # f_name = 'poolOut.txt'
    # out = np.genfromtxt(f_name, delimiter=',')
    # x = out[:,0]
    # y = out[:,1]
    # phi = out[:,2]

    # n = int(np.sqrt(phi.size))
    # q = ax1.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)))
    plt.show()
    # CB = fig.colorbar(img, shrink=0.8, extend='both')
elif mode == 8:
    f_name += 'MSSNodes'

    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]

    plt.plot(x, y, marker='o', linestyle='None', color='r')

    plt.show()
elif mode == 9:
    f_temp = f_name
    f_name += 'poolOut'
    out = np.genfromtxt(f_name, delimiter=',')

    x = out[:,0]
    y = out[:,1]
    phi = out[:,2]

    n = int(np.sqrt(phi.size))

    fig, ax = plt.subplots(1)
    img = ax.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)), levels=[0], colors='b')

    f_name = f_temp

    for i in range(numObj):
        f_name += 'MSSEdges{0}'.format(i)

        out = np.genfromtxt(f_name, delimiter=',')
        x = out[:,0]
        y = out[:,1]

        for i in range(0, x.size, 2):
            plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=0.5)
        plt.xlim((0, 4))
        plt.ylim((0, 2))

        f_name = f_temp


    plt.show()
elif mode == 10:
    f_temp = f_name
    f_name += 'poolOut'
    out = np.genfromtxt(f_name, delimiter=',')

    # nx = 100
    # ny = 100

    x = out[:,0]
    y = out[:,1]
    phi = out[:,2]

    fig, ax = plt.subplots(1)

    # ax.contour(x, y, phi)
    
    for i in range(numObj):
        f_name = f_temp
        f_name += 'MSSEdges{0}'.format(i)

        out = np.genfromtxt(f_name, delimiter=',')
        x_mss = out[:,0]
        y_mss = out[:,1]

        # x_list = [x[i:i+2] for i in range(x.size, 2)]
        # y_list = [y[i:i+2] for i in range(y.size, 2)]
        lines = []
        # print(x.size)
        # print(list(range(0, x.size, 2)))
        # for j in range(0, x.size, 2):
        #     # print('hi')
        #     lines.append([(x[j],y[j]), (x[j+1],y[j+1])])
        # # print('done for')
        # # lines = [[(x[i],y[i]), (x[i+1],y[i+1])] for i in range(x.size, 2)]

        # # print(lines[0])
        # # assert(False)
        # from matplotlib.collections import LineCollection
        # lc = LineCollection(lines, colors='r')

        # ax.add_collection(lc)

        # for i in range(0, x.size, 2):
        #     plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=1, lw=0.5)
        plt.plot(x_mss, y_mss, 'ro-', ms=1, lw=0.5)
        # plt.xlim((0.25, 1.75))
        # plt.ylim((1, 2))
    
    f_name = f_temp
    # f_name += 'poolVel'
    out = np.genfromtxt(f_name+'out', delimiter=',')
    x = out[:,0]
    y = out[:,1]
    u = out[:,2]
    v = out[:,3]
    # p = out[:,4]

    # v = v.flatten()



    temp = np.sqrt(u**2 + v**2)

    q = ax.quiver(x, y, u, v)
    # plt.imshow(np.reshape(u, (n, n)))
    # heat_map = sb.heatmap(np.reshape(temp, (n, n)))

    out = np.genfromtxt(f_name+'medialAxis', delimiter=',')

    if len(out) > 0:
        x = out[:,0]
        y = out[:,1]

        plt.scatter(x, y)
    
    # plt.title('$t = 1.73$')

    plt.gca().set_aspect('equal')
    # plt.axis('off')
    if testNum == -1:
        plt.savefig("{0}.png".format(test), bbox_inches='tight')
    else:
        plt.savefig("{0}{1}.png".format(test, testNum), bbox_inches='tight')
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    plt.show()
elif mode == 11:
    f_temp = f_name
    f_name += 'medialAxis'
    out = np.genfromtxt(f_name, delimiter=',')

    x = out[:,0]
    y = out[:,1]

    plt.scatter(x, y)

    for i in range(numObj):
        f_name = f_temp + 'MSSEdges{0}'.format(i)

        out = np.genfromtxt(f_name, delimiter=',')
        x = out[:,0]
        y = out[:,1]

        for i in range(0, x.size, 2):
            plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=0.5)
        # plt.xlim((0, 1)) plt.ylim((0, 1))

        f_name = f_temp

    plt.gca().set_aspect('equal')
    plt.show()

elif mode == 12:  # output both fluid velocity and object velocity plots
    plots = [('Fluid Velocity', 'out'), ('Object Velocity', 'poolVel')]

    cur_f_name = f_name + 'poolOut'
    out = np.genfromtxt(cur_f_name, delimiter=',')

    x_pool = out[:,0]
    y_pool = out[:,1]
    phi = out[:,2]

    x_min = min(x_pool)
    x_max = max(x_pool)

    y_min = min(y_pool)
    y_max = max(y_pool)

    x_low_ind = np.where(x_pool == x_min)
    x_high_ind = np.where(x_pool == x_max)
    
    y_low_ind = np.where(y_pool == y_min)
    y_high_ind = np.where(y_pool == y_max)

    nx = x_high_ind - x_low_ind + 1
    ny = y_high_ind - y_low_ind + 1


    # n = int(np.sqrt(phi.size))

    for idx, (name, file_ext) in enumerate(plots, start=1):
        fig, ax = plt.subplots(1)
        # img = ax.contour(np.reshape(x_pool, (n, n)), np.reshape(y_pool, (n, n)), np.reshape(phi, (n, n)), levels=[0], colors='b')

        f = plt.figure(idx)
        for i in range(numObj):
            mss_edges_f_name = f_name + 'MSSEdges{0}'.format(i)

            out = np.genfromtxt(mss_edges_f_name, delimiter=',')
            x = out[:,0]
            y = out[:,1]

            # for i in range(0, x.size, 2):
            #     plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=2, lw=0.5)
            # plt.xlim((0, 1))
            # plt.ylim((0, 1))
            x_list = [x[i:i+2] for i in range(x.size, 2)]
            y_list = [y[i:i+2] for i in range(y.size, 2)]
            lines = []
            # print(x.size)
            # print(list(range(0, x.size, 2)))
            for j in range(0, x.size, 2):
                # print('hi')
                lines.append([(x[j],y[j]), (x[j+1],y[j+1])])
            # print('done for')
            # lines = [[(x[i],y[i]), (x[i+1],y[i+1])] for i in range(x.size, 2)]

            # print(lines[0])
            # assert(False)
            from matplotlib.collections import LineCollection
            lc = LineCollection(lines, colors='r')

            ax.add_collection(lc)

        field_vels_f_name = f_name + file_ext
        plt.title(test + ' - ' + name)

        out = np.genfromtxt(field_vels_f_name, delimiter=',')
        x = out[:,0]
        y = out[:,1]
        u = out[:,2]
        v = out[:,3]

        q = ax.quiver(x, y, u, v)
        # plt.imshow(np.reshape(u, (n, n)))
        # heat_map = sb.heatmap(np.reshape(temp, (n, n)))
        # plt.title('Extrapolated Speed Field')

        plt.gca().set_aspect('equal')
        plt.axis('off')
    
    plt.show()
elif mode == 13:  # plot snapshots of steps through simulation
    print('mode 13')
    tracers = []

    # code for displaying plot
    step_num = 1
    def displayPlot():
        for tracer in tracers:
            tracer = (tracer[0], tracer[1])
            plt.scatter(tracer[0], tracer[1], c='b')

        x,y,x_pool,y_pool,x_med,y_med,u,v,steps = plots[curr_pos]
        print('steps = {}'.format(steps))
        plt.scatter(x, y, c='r', marker='o', s=1)
        # if (len(x_med) > 0):
        #     plt.scatter(x_med, y_med, s=0.5)
        # q = ax.quiver(x_pool, y_pool, u, v)
        # plt.title(test + ', nstep = ' + steps)
        plt.gca().set_aspect('equal')
        # plt.savefig('{0}step{1}.png'.format(f_name, step_num))
        fig.canvas.draw()
    

    # code for navigating using left, right arrow keys
    curr_pos = 0
    def key_event(e):
        global curr_pos

        print(curr_pos)
        if e.key == "right":
            curr_pos = curr_pos + 1
        elif e.key == "left":
            curr_pos = curr_pos - 1
        else:
            return
        curr_pos = curr_pos % len(plots)

        ax.cla()
        displayPlot()


    # get data from each step
    steps = [x for x in os.listdir(f_name) if x.isdigit()]
    steps.sort(key=float)
    steps = steps[::1]
    print(steps)

    plots = []
    max_x = -np.inf
    min_x = np.inf
    max_y = -np.inf
    min_y = np.inf
    for step in steps:


        f_temp = f_name + str(step)
        f_temp += '/MSSTracers'
        out = np.genfromtxt(f_temp, delimiter=',')

        for tracer in out:
            tracers.append((tracer[0], tracer[1]))

        curr_name = f_name + str(step) + '/poolOut'
        print('currName = {0}'.format(curr_name))
        out = np.genfromtxt(curr_name, delimiter=',')

        numObj = len(glob.glob(f_name + str(step) + '/MSSEdges*'))
        # numObj -= 10

        x_objs = []
        y_objs = []

        for i in range(numObj):
            mss_edges_f_name = f_name + str(step) + '/MSSNodes{0}'.format(i)

            out = np.genfromtxt(mss_edges_f_name, delimiter=',')
            x_temp = out[:,0]
            y_temp = out[:,1]
            x_objs = np.concatenate((x_objs, x_temp))
            y_objs = np.concatenate((y_objs, y_temp))

        out = np.genfromtxt(f_name+str(step)+'/out', delimiter=',')
        # out = np.genfromtxt(f_name+str(step)+'/poolVel', delimiter=',')
        # Generate array of indices (for the velocities)
        n = out.shape[0]
        inds = np.arange(out.shape[0])
        np.random.shuffle(inds)

        off = 2#int(n/10)
        # off = int(1)
        # # nx = 840
        # # ny = 420
        # nx = 320
        # ny = 640
        # nx = 320
        # ny = 640
        nx = 640
        ny = 1280
        # print(out.shape)

        mask = np.arange(nx) % off == 0
        mask = np.tile(mask, ny)

        out_fluid = out[mask,:]
        # out_fluid = out

        
        # out_fluid = np.reshape(out, (5, ny, nx))[:,::off,::off].reshape(-1, 5)

        x_pool = out_fluid[:,0]
        y_pool = out_fluid[:,1]
        max_x = max(max(x_pool), max_x)
        min_x = min(min(x_pool), min_x)
        max_y = max(max(y_pool), max_y)
        min_y = min(min(y_pool), min_y)
        u = out_fluid[:,2]
        v = out_fluid[:,3]

        # u /= np.sqrt(u**2 + v**2)
        # v /= np.sqrt(u**2 + v**2)

        out = np.genfromtxt(f_name+str(step)+'/medialAxis', delimiter=',')

        if (out.size == 0) :
            x_med = np.array([])
            y_med = np.array([])
        elif (out.size == 3):
            x_med = np.array([out[0]])
            y_med = np.array([out[1]])
        else:
            x_med = out[:,0]
            y_med = out[:,1]

        plots.append((x_objs,y_objs,x_pool,y_pool,x_med,y_med,u,v,step))

    fig = plt.figure()
    fig.canvas.mpl_connect('key_press_event', key_event)
    ax = fig.add_subplot(111)
    displayPlot()
    
    plt.show()
elif mode == 14:
    f_temp = f_name
    f_name += 'poolOut'
    out = np.genfromtxt(f_name, delimiter=',')

    x = out[:,0]
    y = out[:,1]
    phi = out[:,2]


    n = int(np.sqrt(phi.size))

    x *= n
    y *= n

    fig, ax = plt.subplots(1)
    img = ax.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)), levels=[0], colors='b')
    
    f_name = f_temp
    f_name += 'out'
    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]
    u = out[:,2]
    v = out[:,3]
    p = out[:,4]

    # plt.imshow(np.reshape(p, (n, n)))
    # heat_map = sb.heatmap(np.flip(np.reshape(p, (n, n)), axis=0))
    heat_map = sb.heatmap(np.reshape(p, (n, n)), axis=0)
    # heat_map = sb.heatmap(np.flip(np.reshape(u, (n, n)), axis=0))
    # heat_map = sb.heatmap(np.flip(np.reshape(v, (n, n)), axis=0))
    plt.title('pressure field')

    plt.gca().set_aspect('equal')
    plt.axis('off')

    plt.show()
elif mode == 15:
    f_temp = f_name
    f_temp += 'poolOut'
    print('PoolOut from {}'.format(f_temp))
    outPool = np.genfromtxt(f_temp, delimiter=',')
    xPool = outPool[:,0]
    yPool = outPool[:,1]
    phi = outPool[:,2]


    xMin = min(outPool[:,0])
    xMax = max(outPool[:,0])
    yMin = min(outPool[:,1])
    yMax = max(outPool[:,1])


    f_temp = f_name
    f_temp += 'MSSTracers'
    out = np.genfromtxt(f_temp, delimiter=',')

    tracers = []
    for tracer in out:
        tracers.append((tracer[0], tracer[1]))

    f_temp = f_name
    f_temp += 'domain'
    out = np.genfromtxt(f_temp, delimiter=',')
    nx = out.shape[1]
    ny = out.shape[0]

    print('nx = {0} ny = {1}'.format(nx, ny))

    # plt.scatter(0.1*nx, 0.9*ny, c='b')




    

    fig, ax = plt.subplots(1)
    # for tracer in tracers:
    #     print('adding tracer {0}'.format(tracer))
    #     tracer = ((nx/xMax)*tracer[0], (ny/yMax)*tracer[1])
    #     # tracer = (tracer[0], tracer[1])
    #     print('adding tracer {0}'.format(tracer))
    #     plt.scatter(tracer[0], ny - tracer[1], c='y')
    plt.imshow(out[::-1, :])
    xPool *= nx/(xMax - xMin)
    yPool *= ny/(yMax - yMin)
    yPool = ny - yPool
    ax.contour(np.reshape(xPool, (ny, nx)), np.reshape(yPool, (ny, nx)), np.reshape(phi, (ny, nx)), levels=[0], colors='b')

    # for i in range(numObj):
    #     f_temp = f_name
    #     f_temp += 'MSSEdges{0}'.format(i)

    #     out = np.genfromtxt(f_temp, delimiter=',')
    #     x_mss = (nx/(xMax-xMin)) * out[:,0]
    #     y_mss = (ny/(yMax-yMin)) * out[:,1]
    #     # print(x_mss)

    #     plt.scatter(x_mss, ny - y_mss, c='y', s=0.5)

    # for tracer in tracers:
    #     print('adding tracer {0}'.format(tracer))
    #     tracer = ((nx/(xMax-xMin))*tracer[0], (ny/(yMax-yMin))*tracer[1])
    #     # tracer = (tracer[0], tracer[1])
    #     print('adding tracer {0}'.format(tracer))
    #     plt.scatter(tracer[0], ny - tracer[1], c='y', s=0.5)

    plt.show()
    
else:
    print('Invalid mode!')
    sys.exit()

"""
end pool
"""
