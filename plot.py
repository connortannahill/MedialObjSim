import matplotlib.pyplot as plt
import numpy as np
import itertools
import seaborn as sb
import sys
import os

if len(sys.argv) == 1:
    print('Must provide a plotting mode!')
    sys.exit()

test = sys.argv[1]
mode = int(sys.argv[2])
f_name = './output/{0}/'.format(test)

numObj = 1
if (len(sys.argv) == 4):
    numObj = int(sys.argv[3])


"""
For plotting the Pool isocontour
"""

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
        plt.xlim((0, 1))
        plt.ylim((0, 1))

        f_name = f_temp


    plt.show()
elif mode == 10:
    # plt.add_subplot(111,aspect='equal')
    # plt.axes().set_aspect('equal', 'datalim')
    f_temp = f_name
    f_name += 'poolOut'
    out = np.genfromtxt(f_name, delimiter=',')

    nx = 100
    ny = 100

    x = out[:,0]
    y = out[:,1]
    phi = out[:,2]

    n = int(np.sqrt(phi.size))

    # x *= n
    # y *= n

    fig, ax = plt.subplots(1)
    img = ax.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)), levels=[0], colors='b')

    # f_name = f_temp
    # f_name += 'MSSTracers'
    # out = np.genfromtxt(f_name, delimiter=',')
    # x = out[:,0]
    # y = out[:,1]

    # for i in range(out.shape[0]):
    #     plt.scatter(x[i], y[i], color='y')
    
    for i in range(numObj):
        f_name = f_temp
        f_name += 'MSSEdges{0}'.format(i)

        out = np.genfromtxt(f_name, delimiter=',')
        x = out[:,0]
        y = out[:,1]

        for i in range(0, x.size, 2):
            plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=4, lw=0.5)
        plt.xlim((0, 1))
        plt.ylim((0, 1))
    
    f_name = f_temp
    f_name += 'poolVel'
    # f_name += 'out'
    out = np.genfromtxt(f_name, delimiter=',')
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
    # plt.title('Extrapolated Speed Field')

    plt.gca().set_aspect('equal')
    plt.axis('off')
    # plt.savefig("test.png", bbox_inches='tight')


    plt.show()
elif mode == 11:
    f_name += 'medialAxis'
    out = np.genfromtxt(f_name, delimiter=',')

    x = out[:,0]
    y = out[:,1]

    plt.scatter(x, y)
    plt.show()
elif mode == 12:  # output both fluid velocity and object velocity plots
    plots = [('Fluid Velocity', 'out'), ('Object Velocity', 'poolVel')]

    cur_f_name = f_name + 'poolOut'
    out = np.genfromtxt(cur_f_name, delimiter=',')

    nx = 100
    ny = 100

    x_pool = out[:,0]
    y_pool = out[:,1]
    phi = out[:,2]

    n = int(np.sqrt(phi.size))

    for idx, (name, file_ext) in enumerate(plots, start=1):
        fig, ax = plt.subplots(1)
        img = ax.contour(np.reshape(x_pool, (n, n)), np.reshape(y_pool, (n, n)), np.reshape(phi, (n, n)), levels=[0], colors='b')

        f = plt.figure(idx)
        for i in range(numObj):
            mss_edges_f_name = f_name + 'MSSEdges{0}'.format(i)

            out = np.genfromtxt(mss_edges_f_name, delimiter=',')
            x = out[:,0]
            y = out[:,1]

            for i in range(0, x.size, 2):
                plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=2, lw=0.5)
            plt.xlim((0, 1))
            plt.ylim((0, 1))

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

    # code for displaying plot
    def displayPlot():
        x,y,steps = plots[curr_pos]
        # for i in range(0, x.size, 2):
        #     plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=2, lw=0.5)
        plt.plot(x, y, marker='o', linestyle='None', color='r', ms=2)
        plt.title(test + ', nstep = ' + steps)
        plt.xlim((0, 1))
        plt.ylim((0, 1))
        plt.gca().set_aspect('equal')
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

    # x = out[:,0]
    # y = out[:,1]
    # phi = out[:,2]

    # n = int(np.sqrt(phi.size))

    # for idx, (name, file_ext) in enumerate(plots, start=1):
    #     fig, ax = plt.subplots(1)
    #     img = ax.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)), levels=[0], colors='b')
    
    # get data from each step
    steps = os.listdir(f_name)
    steps.sort(key=float)
    print(steps)

    plots = []
    for step in steps:
        curr_name = f_name + step + '/poolOut'
        out = np.genfromtxt(curr_name, delimiter=',')

        # x = out[:,0]
        # y = out[:,1]

        x = []
        y = []

        for i in range(numObj):
            mss_edges_f_name = f_name + step + '/MSSNodes{0}'.format(i)

            out = np.genfromtxt(mss_edges_f_name, delimiter=',')
            x_temp = out[:,0]
            y_temp = out[:,1]
            x = np.concatenate((x, x_temp))
            y = np.concatenate((y, y_temp))

            # for i in range(0, x.size, 2):
            #     plt.plot(x[i:i+2], y[i:i+2], 'ro-', ms=2, lw=0.5)
            # plt.xlim((0, 1))
            # plt.ylim((0, 1))

        # print(x)

        plots.append((x,y,step))

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

else:
    print('Invalid mode!')
    sys.exit()

"""
end pool
"""
