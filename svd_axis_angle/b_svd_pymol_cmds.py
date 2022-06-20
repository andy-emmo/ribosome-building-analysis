import numpy as np


def Rx(a):
    return np.matrix([[1, 0, 0],
                      [0, np.cos(a), -np.sin(a)],
                      [0, np.sin(a), np.cos(a)]])


def Ry(a):
    return np.matrix([[np.cos(a), 0, np.sin(a)],
                      [0, 1, 0],
                      [-np.sin(a), 0, np.cos(a)]])


def Rz(a):
    return np.matrix([[np.cos(a), -np.sin(a), 0],
                      [np.sin(a), np.cos(a), 0],
                      [0, 0, 1]])


def euler_svd2ax(vec1, vec2):
    # what is the angle between yz and z=0
    theta = np.arctan(vec1[1] / vec1[2])
    # so I know I want to rotate around the x-axis by this amount
    psi1 = theta  # save this angle for later
    R = Rz(0) * Ry(0) * Rx(psi1)
    rot_Rzyz = np.matmul(R, vec1)
    #print(psi1)
    #print(rot_Rzyz)

    # now rotate the vector to align with the z-axis
    # we need to find the arctan of x/z
    theta = np.arctan(rot_Rzyz[0, 0] / rot_Rzyz[0, 2])
    # we need to move in the opposite direction
    theta1 = -theta
    R = Rz(0) * Ry(theta1) * Rx(psi1)
    rot_Rzyz = np.matmul(R, vec1)
    #print(theta1)
    #print(rot_Rzyz)

    # The rotation matrix
    R = Rz(0) * Ry(theta1) * Rx(psi1)
    #print(R)
    # now using the secondary (orthogonal) vector, rotate around z axis
    n2_rot = np.copy(np.matmul(R, vec2)[0]).squeeze()
    #print(n2_rot)
    #print(np.linalg.norm(n2_rot))
    #print(n2_rot.shape)
    phi1 = -np.arctan(n2_rot[1]/n2_rot[0])
    #print(phi1)
    R = Rz(phi1) * Ry(theta1) * Rx(psi1)
    #print(np.matmul(R, vec1))
    #print(np.matmul(R, vec2))
    return psi1, theta1, phi1


def do_svd(fn, show_plot=False, verbose=False):
    import json
    from la_functions import planeFit2, planeFit
    
    # open the coordinates we dumped previously
    # they are saved as a list
    f = open(fn)
    data = f.read()
    f.close()
    a = json.loads(data)
    # convert to numpy array
    a = np.array(a)
    if verbose:
        print('Shape of the input array: ', a.shape)
    # fit a plane to the points (SVD) 
    point, normals = planeFit2(a.T)
    if verbose:
        print('SVD point/plane fit:')
        print(point)
        print(normals)
    n1 = normals[-1]
    # then we require the secondary axis to get the proper rotation around z
    n2 = normals[1]
    # determine the Euler angles to align the least variability in our point cloud to the (0,0,1) axis
    # why do this? It makes the transformation matrix generated later when comparing state transitions
    # to have the centre of rotation close to the origin and the rotation axis aligned closely to the (0,0,1) axis
    angles = euler_svd2ax(n1, n2)
    point *= -1  # to subtract COM from coords
    
    if verbose:
        for r in angles:
            print('Euler angles (radians, degrees):')
            print(r, np.rad2deg(r))
    
    # output the required pymol commands
    print('\nPymol commands: ')
    print('cmd.translate({0}, "all", camera=0)'.format(point.tolist()))
    print('cmd.rotate("x", {0}, "all", camera=0, origin=[0,0,0])'.format(np.rad2deg(angles[0])))
    print('cmd.rotate("y", {0}, "all", camera=0, origin=[0,0,0])'.format(np.rad2deg(angles[1])))
    print('cmd.rotate("z", {0}, "all", camera=0, origin=[0,0,0])'.format(np.rad2deg(angles[2])))

    print('\nNote, some axes may need to be rotate 180 deg, for example:\ncmd.rotate("z", 180, "all", camera=0, origin=[0,0,0])')

    if show_plot:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_set_axes_equal import set_axes_equal
        
        # move the centre of mass to the origin
        a += point
        R = Rz(angles[2]) * Ry(angles[1]) * Rx(angles[0])
        # rotate all the points, giving the a2 coordinate set
        a2 = np.matmul(R, a.T)
        a2 = np.array(a2.T)
        #print(planeFit(a2.T))
        #print(planeFit2(a2.T))
        
        # plotting
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(*a.T, color='C0')
        ax.scatter(*a2.T, color='C1')

        # label the axes
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        set_axes_equal(ax)
        plt.show()


if __name__ == '__main__':
    do_svd('a_coords.json')
    
