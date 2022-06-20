import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_set_axes_equal import set_axes_equal
from la_functions import aff2axangle, angle_between
import json


def get_rot_point(tr_mat, plane='xy'):
    from scipy.optimize import minimize
    def rotate(x):
        if plane == 'xy':
            point_i = np.array([x[0], x[1], 0., 1.])
        point_f = np.matmul(tr_mat, point_i)
        length = np.linalg.norm(point_f-point_i)
        #print(length)
        return length
    result = minimize(rotate, np.array((0, 0)))
    if result.success:
        point = np.array([result.x[0], result.x[1], 0])
        return point
    else:
        raise(Exception('Minimization has failed'))


def get_all_axangle(fn, show_plot=False):
    # load the transformation matrix(ces)
    f = open(fn)
    all_info = json.loads(f.read())
    f.close()
    print(all_info)
    print(all_info.keys())
    points = [[], []]
    rot_axes = []
    for k in all_info.keys():
        # load the transformation matrix
        print(k)
        tr_mat = all_info[k]['tr_matrix']
        tr_mat = np.array(tr_mat)
        #print(tr_mat)
        tr_mat = tr_mat.reshape([4, 4])
        print(tr_mat)
        # determine the axis of rotation, angle of rotation about this axis
        axis, angle, point1 = aff2axangle(tr_mat)
        # the point from this function is generally "rubbish"
        # it lies on the axis, so is mathematically correct,
        #print('Axis, angle (rad), angle (deg), point')
        #print(axis, angle, np.rad2deg(angle), point1)
        # we'll determine this point using scipy.minimize
        # of the point that passes through a given plane (e.g. the xy plane)
        point = get_rot_point(tr_mat)
        print('Axis, angle (rad), angle (deg), point')
        print(axis, angle, np.rad2deg(angle), point)
        points[0].append(point[0])
        points[1].append(point[1])
        rot_axes.append(axis)
    
    # plot the centre of rotation?
    if show_plot:
        plt.scatter(*points)
        plt.show()


if __name__ == '__main__':
   get_all_axangle('./out_transformation_matrices.json') 

   
