import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob

# base_dir = './output/3D'

if len(sys.argv) == 1:
    print('Must provide a plotting mode!')
    sys.exit()

testNum = -1

test = sys.argv[1]
mode = int(sys.argv[2])
print("mode = {}".format(mode))

if (mode != 4):
    if (len(sys.argv) == 4):
        testNum = int(sys.argv[3])
elif mode == 4:
    print(len(sys.argv))
    if (len(sys.argv) == 6):
        testNum = int(sys.argv[3])


print('testNum = {}'.format(testNum))

f_name = None
if testNum == -1:
    f_name = './output/{0}/'.format(test)
else:
    f_name = './output/{0}/{1}/'.format(test, testNum)


nMSS= len(glob.glob(f_name+'MSS3DEdges*'))
print(glob.glob(f_name+'MSS3DEdges*'))
# nMss = int(sys.argv[3])
# narg = 1
# if len(sys.argv) == 3:
#     narg = int(sys.argv[2])


# fluid_fname = 'out3D'
# pool_fname = 'poolVals3D'
# pool_fname = 'poolVel3D'

if mode == -1:
    # f_name = 'colOut.txt'
    f_name += 'medialAxis3D'

    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]
    z = out[:,2]

    fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z)])

    fig.show()

if mode == -2:
    f_name += 'colOut'

    out = np.genfromtxt(f_name, delimiter=',')
    x = out[:,0]
    y = out[:,1]
    z = out[:,2]
    u = out[:,3]
    v = out[:,4]
    w = out[:,5]
    

    fig = go.Figure(data=[go.Cone(x=x, y=y, z=z,
                                  u=u, v=v, w=w)])

    fig.show()

def plot_3D(num):
    f_temp = f_name

    pool_iso_fname = f_temp + 'pool3DOut'
    # pool_vel_fname = 'pool3DVel'
    pool_vel_fname = f_temp + 'out3D'

    out_pool = np.genfromtxt(pool_iso_fname, delimiter=',')
    out_pool_vel = np.genfromtxt(pool_vel_fname, delimiter=',')

    # Generate array of indices (for the velocities)
    n = out_pool_vel.shape[0]
    inds = np.arange(out_pool_vel.shape[0])
    np.random.shuffle(inds)

    # off = n#int(n/10)
    off = int(n/100)
    out_fluid = out_pool_vel[inds[:off],:]

    x = out_fluid[:,0]
    y = out_fluid[:,1]
    z = out_fluid[:,2]
    u = out_fluid[:,3]
    v = out_fluid[:,4]
    w = out_fluid[:,5]
    p = out_fluid[:,6]

    x_pool = out_pool[:,0]
    y_pool = out_pool[:,1]
    z_pool = out_pool[:,2]
    phi    = out_pool[:,3]
    # cx = 1
    # cy = 1
    # cz = 1
    # a = 0.3
    # c = 0.105
    # r = 1
    # b = 2.25 * r

    # for i, pnt in enumerate(zip(x_pool, y_pool, z_pool)):
    #     x, y, z = pnt
    #     xsq = ((x-cx)/b)**2
    #     ysq = ((y-cy)/b)**2
    #     zsq = ((z-cz)/b)**2
    #     phi[i] = (xsq + ysq + zsq + a**2)**2 - 4*(a**2)*(xsq + ysq) - c**2


    # fig = go.Figure(data=[
    #     go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, sizemode="absolute"),
    #     go.Isosurface(x=x_pool, y=y_pool, z=z_pool, value=phi, isomin=-0.0001,
    #                     isomax=0.0001, colorscale=[(0, 'blue'), (1, 'red')])])
    fig = go.Figure(data=[
        go.Isosurface(x=x_pool, y=y_pool, z=z_pool, value=phi, isomin=-0.0001,
                        isomax=0.0001, colorscale=[(0, 'blue'), (1, 'red')])])

    fig.show()

    # import plotly.io as pio
    # fig.write_image("{0}{1}.png".format(f_temp, test))


if mode == 0:
    if len(sys.argv) == 4:
        plot_3D(int(sys.argv[3]))
    else:
        plot_3D(-1)
elif mode == 1:
    pool_fname = f_name + 'pool3DOut'

    # Read in the pool data
    out = np.genfromtxt(pool_fname, delimiter=',')

    x = out[:,0]
    y = out[:,1]
    z = out[:,2]
    phi = out[:,3]

    vel_fname = f_name + 'out3D'

    # Read in the pool data

    out_pool_vel = np.genfromtxt(vel_fname, delimiter=',')

    n = out_pool_vel.shape[0]
    inds = np.arange(out_pool_vel.shape[0])
    np.random.shuffle(inds)

    # off = n#int(n/10)
    off = int(n/2000)
    out = out_pool_vel[inds[:off],:]


    x_vel = out[:,0]
    y_vel = out[:,1]
    z_vel = out[:,2]
    u = out[:,3]
    v = out[:,4]
    w = out[:,5]

    # velData = [go.Cone(x=x_vel, y=y_vel, z=z_vel, u=u, v=v, w=w, sizemode="scaled")]
    velData = []

    mssList = []

    print(nMSS)
    
    for i in range(nMSS):
        f_name_t = f_name + 'MSS3DNodes{0}'.format(i)

        out = np.genfromtxt(f_name_t, delimiter=',')
        x_mss = out[:,0]
        y_mss = out[:,1]
        z_mss = out[:,2]

        mssList.append(go.Scatter3d(x=x_mss, y=y_mss, z=z_mss, mode='markers'))

        # Now append the object velocities

        # f_name_t = f_name + 'MSS3DNodes{0}'.format(i)

        # out = np.genfromtxt(f_name_t, delimiter=',')
        # x_mss = out[:,0]
        # y_mss = out[:,1]
        # z_mss = out[:,2]
    # from matplotlib import pyplot
    # from mpl_toolkits.mplot3d import Axes3D
    # import random


    # fig = pyplot.figure()
    # ax = Axes3D(fig)

    # ax.scatter(x_mss, y_mss, z_mss)
    # pyplot.show()


        # for i in range(0, x.size, 3):
        # mssList.append(go.Scatter3d(x=x[i:i+3], y=y[i:i+3], z=z[i:i+3]))
            # ax.plot(x[i:i+2], y[i:i+2], z[i:i+2], 'ro-', ms=0.5)
    # for i in range(1):``
    #     f_name_t = f_name + 'MSS3DEdges{0}'.format(i)

    #     out = np.genfromtxt(f_name_t, delimiter=',')
    #     x = out[:,0]
    #     y = out[:,1]
    #     z = out[:,2]

    #     for i in range(0, x.size, 2):
    #         ax.plot(x[i:i+2], y[i:i+2], z[i:i+2], 'ro-', ms=0.5)

    # fig = go.Figure(data=[go.Isosurface(
    #     x=x,
    #     y=y,
    #     z=z,
    #     value=phi,
    #     isomin=-0.0001,
    #     isomax=0.0001,
    #     colorscale=[(0,"blue"), (1,"red")],
    #     opacity=0.3), go.Scatter3d(x=x_mss, y=y_mss, z=z_mss)])
    fig = go.Figure(data=[go.Isosurface(
        x=x,
        y=y,
        z=z,
        value=phi,
        isomin=-0.0001,
        isomax=0.0001,
        colorscale=[(0,"blue"), (1,"red")],
        opacity=0.3)]+mssList+velData)
    # fig = go.Figure(data=mssList)

    fig.show()
elif mode == 2:
    pool_vel_fname = f_name + 'pool3DVel'
    out = np.genfromtxt(pool_vel_fname, delimiter=',')

    x = out[:,0]
    y = out[:,1]
    z = out[:,2]
    u = out[:,3]
    v = out[:,4]
    w = out[:,5]

    fig = go.Figure(
        data=go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, sizemode="absolute"))

    fig.show()
elif mode == 3:
    """ Use MPL 3D plotting to plot line segments of MSS """
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Read in input data

    for i in range(1):
        f_name_t = f_name + 'MSS3DEdges{0}'.format(i)

        out = np.genfromtxt(f_name_t, delimiter=',')
        x = out[:,0]
        y = out[:,1]
        z = out[:,2]

        for i in range(0, x.size, 2):
            ax.plot(x[i:i+2], y[i:i+2], z[i:i+2], 'ro-', ms=0.5)
        # ax.xlim((0, 1))
        # ax.ylim((0, 1))
    
    # Scatter plot of the centroids
    f_name_t = f_name + 'MSS3DNodes{}'.format(0)
    out = np.genfromtxt(f_name_t, delimiter=',')

    x = out[:,0]
    y = out[:,1]
    z = out[:,2]

    ax.scatter(x, y, z, c='b')
    
    plt.show()
elif mode == 4:
    """
    Plotting a 2D slice of the velocity field

    INPUT:
    python plot3D.py 4 <axis=(x|y|z)> <val>

    For the given axis and val, we interpolate along this axis with
    constant val
    """

    axis = sys.argv[4]
    if (axis != 'x' and axis != 'y' and axis != 'z') :
        raise ValueError("axis {} does not exist!".format(axis))

    val = float(sys.argv[5])

    # Read in the input values
    pool_vel_fname = f_name + 'out3D'
    # pool_vel_fname = 'out3D'
    out = np.genfromtxt(pool_vel_fname, delimiter=',')

    # The coordinates
    x = out[:,0]
    y = out[:,1]
    z = out[:,2]

    # The velocity vector components
    u = out[:,3]
    v = out[:,4]
    w = out[:,5]
    p = out[:,6]

    print('max u = {}'.format(np.max(np.abs(u))))
    print('max v = {}'.format(np.max(np.abs(v))))
    print('max w = {}'.format(np.max(np.abs(w))))

    # Get the x, y, z domain of the point cloud.
    x_bounds = [np.amin(x), np.amax(x)]
    y_bounds = [np.amin(y), np.amax(y)]
    z_bounds = [np.amin(z), np.amax(z)]

    # Extract the phi value
    pool_vel_fname = f_name + 'pool3DOut'
    out = np.genfromtxt(pool_vel_fname, delimiter=',')

    x_pool = out[:,0]
    y_pool = out[:,1]
    z_pool = out[:,2]

    x_min = min(x_pool)
    x_max = max(x_pool)

    y_min = min(y_pool)
    y_max = max(y_pool)

    z_min = min(z_pool)
    z_max = max(z_pool)

    x_low_ind = np.where(x_pool == x_min)[0][0]
    x_high_ind = np.where(x_pool == x_max)[0][0]
    
    y_low_ind = np.where(y_pool == y_min)[0][0]
    y_high_ind = np.where(y_pool == y_max)[0][0]

    z_low_ind = np.where(z_pool == z_min)[0][0]
    z_high_ind = np.where(z_pool == z_max)[0][0]

    nx = x_high_ind - x_low_ind + 1
    ny = int((y_high_ind - y_low_ind + 1)/(nx)+1)
    nz = int((z_high_ind - z_low_ind + 1)/((nx)*(ny))+1)

    print('nx = {0} ny = {1} nz = {2}'.format(nx, ny, nz))
    print('tot = {}'.format(x_pool.size))

    phi = out[:,3]


    # Find the closest value that is in the array to val.
    closest = -1.0

    fig, ax1 = plt.subplots()

    if axis == 'x':
        closest = np.argmin(np.abs(x - val))
        closest_pnt = x[closest]

        inds = (x == closest_pnt)
        plt.imshow(np.reshape(p[inds], (nz, ny)))
        plt.colorbar()
        # q = ax1.quiver(y[inds], z[inds], v[inds], w[inds])


        # img = ax1.contour(np.reshape(y[inds], (nz, ny)), np.reshape(z[inds], (nz, ny)), np.reshape(phi[inds], (nz, ny)), levels=[0], colors='b')
        # img = ax1.contour(np.reshape(y[inds], (n, n)), np.reshape(z[inds], (n, n)), np.reshape(phi[inds], (n, n)), colors='b')
    elif axis == 'y':
        closest = np.argmin(np.abs(y - val))
        closest_pnt = y[closest]

        inds = (y == closest_pnt)
        q = ax1.quiver(x[inds], z[inds], u[inds], w[inds])

        n = int(np.sqrt(phi[inds].size))

        img = ax1.contour(np.reshape(x[inds], (nz, nx)), np.reshape(z[inds], (nz, nx)), np.reshape(phi[inds], (nz, nx)), levels=[0], colors='b')
        # img = ax1.contour(np.reshape(x[inds], (n, n)), np.reshape(z[inds], (n, n)), np.reshape(phi[inds], (n, n)), colors='b')
    elif axis == 'z':
        closest = np.argmin(np.abs(z - val))
        print(closest)
        closest_pnt = z[closest]
        print(closest_pnt)

        inds = (z == closest_pnt)
        print(np.sum(inds))
        q = ax1.quiver(x[inds], y[inds], u[inds], v[inds])


        n = int(np.sqrt(phi[inds].size))
        img = ax1.contour(np.reshape(x[inds], (ny, nx)), np.reshape(y[inds], (ny, nx)), np.reshape(phi[inds], (ny, nx)), levels=[0], colors='b')
        # img = ax1.contour(np.reshape(x[inds], (n, n)), np.reshape(y[inds], (n, n)), np.reshape(phi[inds], (n, n)), colors='b')
    
    plt.show()
elif mode == 5:
    mssNum = 1
    mss_vel_fname = f_name +  'MSS3DVels{0}'.format(mssNum)

    out = np.genfromtxt(mss_vel_fname, delimiter=',')

    # The coordinates
    x = out[:,0]
    y = out[:,1]
    z = out[:,2]

    # The velocity vector components
    u = out[:,3]
    v = out[:,4]
    w = out[:,5]


    pool_fname = f_name + 'pool3DOut'

    # Read in the pool data
    out = np.genfromtxt(pool_fname, delimiter=',')

    xPool = out[:,0]
    yPool = out[:,1]
    zPool = out[:,2]
    phi = out[:,3]

    print('Max u = {0} Max v = {1}'.format(np.amax(np.abs(u)), np.amax(np.abs(v))))

    # Create the cone plot for the nodes
    fig = go.Figure(
        data=go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, sizemode="absolute"))

    fig.show()
