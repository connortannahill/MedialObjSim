import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
import sys

base_dir = './output/3D'

if len(sys.argv) == 1:
    print('Must provide a plotting mode!')
    sys.exit()

test = sys.argv[1]
mode = int(sys.argv[2])
f_name = './output/{0}/'.format(test)

nMss = int(sys.argv[3])
narg = 1
if len(sys.argv) == 3:
    narg = int(sys.argv[2])


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

if mode == 0:
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
    off = int(n/10)
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


    fig = go.Figure(data=[
        go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, sizemode="absolute"),
        go.Isosurface(x=x_pool, y=y_pool, z=z_pool, value=phi, isomin=-0.0001,
                        isomax=0.0001, colorscale=[(0, 'blue'), (1, 'red')])])

    fig.show()
elif mode == 1:
    pool_fname = f_name + 'pool3DOut'

    # Read in the pool data
    out = np.genfromtxt(pool_fname, delimiter=',')

    x = out[:,0]
    y = out[:,1]
    z = out[:,2]
    phi = out[:,3]

    mssList = []
    
    for i in range(nMss):
        f_name_t = f_name + 'MSS3DEdges{0}'.format(i)

        out = np.genfromtxt(f_name_t, delimiter=',')
        x_mss = out[:,0]
        y_mss = out[:,1]
        z_mss = out[:,2]
        mssList.append(go.Scatter3d(x=x_mss, y=y_mss, z=z_mss))

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
        opacity=0.3)]+mssList)

    fig.show()
elif mode == 2:
    pool_vel_fname = f_name +  'pool3DVel'.format(0)
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
    f_name_t = f_name + 'MSS3DCentroids{}'.format(0)
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

    axis = sys.argv[2]
    if (axis != 'x' and axis != 'y' and axis != 'z') :
        raise ValueError("axis {} does not exist!".format(axis))

    val = float(sys.argv[3])

    # Read in the input values
    pool_vel_fname = f_name + 'pool3DVel'
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

    # Get the x, y, z domain of the point cloud.
    x_bounds = [np.amin(x), np.amax(x)]
    y_bounds = [np.amin(y), np.amax(y)]
    z_bounds = [np.amin(z), np.amax(z)]

    # Extract the phi value
    pool_vel_fname = f_name + 'pool3DOut'
    out = np.genfromtxt(pool_vel_fname, delimiter=',')

    phi = out[:,3]


    # Find the closest value that is in the array to val.
    closest = -1.0

    fig, ax1 = plt.subplots()

    if axis == 'x':
        closest = np.argmin(np.abs(x - val))
        closest_pnt = x[closest]

        inds = (x == closest_pnt)
        q = ax1.quiver(y[inds], z[inds], v[inds], w[inds])


        img = ax1.contour(np.reshape(y[inds], (n, n)), np.reshape(z[inds], (n, n)), np.reshape(phi[inds], (n, n)), levels=[0], colors='b')
        # img = ax1.contour(np.reshape(y[inds], (n, n)), np.reshape(z[inds], (n, n)), np.reshape(phi[inds], (n, n)), colors='b')
    elif axis == 'y':
        closest = np.argmin(np.abs(y - val))
        closest_pnt = y[closest]

        inds = (y == closest_pnt)
        q = ax1.quiver(x[inds], z[inds], u[inds], w[inds])

        n = int(np.sqrt(phi[inds].size))

        img = ax1.contour(np.reshape(x[inds], (n, n)), np.reshape(z[inds], (n, n)), np.reshape(phi[inds], (n, n)), levels=[0], colors='b')
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
        img = ax1.contour(np.reshape(x[inds], (n, n)), np.reshape(y[inds], (n, n)), np.reshape(phi[inds], (n, n)), levels=[0], colors='b')
        # img = ax1.contour(np.reshape(x[inds], (n, n)), np.reshape(y[inds], (n, n)), np.reshape(phi[inds], (n, n)), colors='b')
    
    plt.show()
elif mode == 5:
    mss_vel_fname = f_name +  'MSS3DVels{0}'.format(0)

    out = np.genfromtxt(mss_vel_fname, delimiter=',')

    # The coordinates
    x = out[:,0]
    y = out[:,1]
    z = out[:,2]

    # The velocity vector components
    u = out[:,3]
    v = out[:,4]
    w = out[:,5]

    print('Max u = {0} Max v = {1}'.format(np.amax(np.abs(u)), np.amax(np.abs(v))))

    # Create the cone plot for the nodes
    fig = go.Figure(
        data=go.Cone(x=x, y=y, z=z, u=u, v=v, w=w, sizemode="absolute"))

    fig.show()
