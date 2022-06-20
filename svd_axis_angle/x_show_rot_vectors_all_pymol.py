from pymol import cmd
import numpy as np
import json
import os
import scipy
from scipy.optimize import minimize

# ffmpeg -i image%04d.png -c:v libx264 -vf fps=10 -pix_fmt yuv420p out.mp4 -y

color_red = [1., 0., 0.]
color_blue = [0., 0., 1.]
color_green = [0., 1., 0.]


def make_cgo_list(p1, p2, color1=[1, 1, 1], color2=[1, 1, 1], radius=0.3):
    x1, y1, z1 = p1  # start point
    x2, y2, z2 = p2  # end point
    r1, g1, b1 = color1  # color start
    r2, g2, b2 = color2  # color end
    cgo_list = [9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2]
    return cgo_list


# this code taken fr transforms3d package
def aff2axangle(aff):
    """Return axis, angle and point fr affine

    Parameters
    ----------
    aff : array-like shape (4,4)

    Returns
    -------
    axis : array shape (3,)
       vector giving axis of rotation
    angle : scalar
       angle of rotation in radians.
    point : array shape (3,)
       point around which rotation is performed

    Examples
    --------
    >>> direc = np.random.random(3) - 0.5
    >>> angle = (np.random.random() - 0.5) * (2*math.pi)
    >>> point = np.random.random(3) - 0.5
    >>> R0 = axangle2aff(direc, angle, point)
    >>> direc, angle, point = aff2axangle(R0)
    >>> R1 = axangle2aff(direc, angle, point)
    >>> np.allclose(R0, R1)
    True

    Notes
    -----
    http://en.wikipedia.org/wiki/Rotation_matrix#Axis_of_a_rotation
    """
    R = np.asarray(aff, dtype=np.float)
    direction, angle = mat2axangle(R[:3, :3])
    # point: unit eigenvector of R33 corresponding to eigenvalue of 1
    L, Q = np.linalg.eig(R)
    i = np.where(abs(np.real(L) - 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError("no unit eigenvector corresponding to eigenvalue 1")
    point = np.real(Q[:, i[-1]]).squeeze()
    point /= point[3]
    return direction, angle, point
    
    
def mat2axangle(mat, unit_thresh=1e-5):
    """Return axis, angle and point fr (3, 3) matrix `mat`

    Parameters
    ----------
    mat : array-like shape (3, 3)
        Rotation matrix
    unit_thresh : float, optional
        Tolerable difference fr 1 when testing for unit eigenvalues to
        confirm `mat` is a rotation matrix.

    Returns
    -------
    axis : array shape (3,)
       vector giving axis of rotation
    angle : scalar
       angle of rotation in radians.

    Examples
    --------
    >>> direc = np.random.random(3) - 0.5
    >>> angle = (np.random.random() - 0.5) * (2*math.pi)
    >>> R0 = axangle2mat(direc, angle)
    >>> direc, angle = mat2axangle(R0)
    >>> R1 = axangle2mat(direc, angle)
    >>> np.allclose(R0, R1)
    True

    Notes
    -----
    http://en.wikipedia.org/wiki/Rotation_matrix#Axis_of_a_rotation
    """
    #print(mat.dtype)
    #M = np.asarray(mat, dtype=np.float)
    M = mat
    # direction: unit eigenvector of R33 corresponding to eigenvalue of 1
    L, W = np.linalg.eig(M.T)
    i = np.where(np.abs(L - 1.0) < unit_thresh)[0]
    if not len(i):
        raise ValueError("no unit eigenvector corresponding to eigenvalue 1")
    direction = np.real(W[:, i[-1]]).squeeze()
    # rotation angle depending on direction
    cosa = (np.trace(M) - 1.0) / 2.0
    if abs(direction[2]) > 1e-8:
        sina = (M[1, 0] + (cosa-1.0)*direction[0]*direction[1]) / direction[2]
    elif abs(direction[1]) > 1e-8:
        sina = (M[0, 2] + (cosa-1.0)*direction[0]*direction[2]) / direction[1]
    else:
        sina = (M[2, 1] + (cosa-1.0)*direction[1]*direction[2]) / direction[0]
    angle = math.atan2(sina, cosa)
    return direction, angle


def fit_line3d(data):
    # Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean = data.mean(axis=0)
    # Do an SVD on the mean-centered data.
    uu, dd, vv = np.linalg.svd(data - datamean)
    # Now vv[0] contains the first principal component, i.e. the direction
    # vector of the 'best fit' line in the least squares sense.
    #print(vv[0])
    spread = int(np.max(np.abs(data - datamean)))*3
    linepts = vv[0] * np.mgrid[-spread:spread:2j][:, np.newaxis]
    linepts += datamean
    return vv[0], datamean, linepts
    

def compare_rot(fr='a0gi_St', to='b0gi_St', skip=1, n=30, 
                compare_self=False, color_states=True, phosphates_only=True):
    print('numpy version', np.__version__)
    print('scipy version', scipy.__version__)
    
    # align (d0_St and (resi 1-894 or resi 1315-1454)), (c0_St and (resi 1-894 or resi 1315-1454))
    
    to_al = to
    fr_al = fr
    
    # align the two states
    # this returns a seven member list, the first member is the RMSD
    rmsd = cmd.align(to_al, fr_al)[0]
    print('RMSD of fit: ', rmsd)
    # get the transformation matrix of the rotated object
    rot = cmd.get_object_matrix(to)
    print('Transformation Matrix:')
    print(rot)
    # get the inverse of the transformation matrix
    rot = np.array(rot)
    rot = rot.reshape([4, 4])
    #print(rot)
    rot_inv = np.linalg.inv(rot)
    #print(rot_inv)
    
    # pymol wants this as a flattened list
    rot = rot.flatten().tolist()
    rot_inv = rot_inv.flatten().tolist()
    #print(rot)
    #print(rot_inv)
    
    # if we want to only look at the phosphates of the backbone
    if phosphates_only:
        sele = '({0}) and name P'
    else:
        sele = '{0}'
    
    # now we have a choice as to whether we want to compare with 
    # the rotate 'to' structure (essentially comparing to itself, but rotated)
    # or with the 'fr' structure (this makes the most sense)
    if compare_self:
        # get the rotated coordinates (selecting only P atoms here)
        cmd.select('sele', sele.format(to))
        xyz_rot = cmd.get_coords('sele', 1)
        # also show the translation of the centre of mass
        com2 = cmd.centerofmass(to)
    else:
        # get the coordinates for the fr structure
        cmd.select('sele', sele.format(fr))
        xyz_rot = cmd.get_coords('sele', 1)
        # also show the translation of the centre of mass
        com2 = cmd.centerofmass(to)
    
    # reverse the transformation
    cmd.transform_object(to, rot_inv)

    # get the coordinates for the non-rotated state
    cmd.select('sele', sele.format(to))
    xyz_notrot = cmd.get_coords('sele', 1)
    # also show the translation of the centre of mass
    com1 = cmd.centerofmass(to)
    
    # now both moving and fr should have null transformation matrices
    # if we just want to look at the change in position of the atoms chosen, 
    # we can stop here and display the relevant information
    
    # determine the angle of rotation and the axis of rotation
    # need to reshape the transformation matrix
    np_rot = np.array(rot)
    np_rot = np_rot.reshape([4, 4])
    rot_axis, angle, point = aff2axangle(np_rot)
    print('Axis Angle: (rot_axis, angle, point)', rot_axis, angle, point)
    
    # also, see if the minimum movers roughly aligns with the rotation axis 

    # find the movement of the COM
    a = np.array(com1)
    b = np.array(com2)
    #print(b-a)
    #print(np.linalg.norm(b-a))
    com_move = np.linalg.norm(b-a)

    norms = np.linalg.norm(xyz_notrot-xyz_rot, axis=1)
    #print(np.max(norms), np.min(norms))
    #print(np.max(norms), np.min(norms), np.mean(norms))

    argsort_norms = np.argsort(norms)
    argsort_min = argsort_norms[:n]
    # get the coordinates of the n minimum norms
    a = np.array(xyz_notrot)
    b = np.array(xyz_rot)
    a_min = a[argsort_min]
    b_min = b[argsort_min]
    ab_min = np.vstack((a_min, b_min))
    #print(a_min)
    #print(b_min)
    #print(norms[argsort_min])
    # use the ab_min points to fit a line in 3d
    line_axis, data_mean, line_points = fit_line3d(ab_min)
    
    # how does it compare to the transformation matrix results?
    print('Comparison of rotation axis and minimum movers axis.')
    print('Axis of rotation: {0}'.format(rot_axis))
    print('Minimum norms axis: {0}'.format(line_axis))
    print('Angle of rotation (deg): {0}'.format(np.rad2deg(angle)))
    print('COM of minimum norms {0}'.format(data_mean))
    #print(data_mean)  # the minimum movers centre of mass
    #print(line_points)  # the extent of the line_axis, to display
    
    # what is the direction of movement of the minimum movers?
    #print(norms)
    minmove_direction = xyz_notrot[argsort_min]-xyz_rot[argsort_min]
    #print(minmove_direction.shape)
    #print(minmove_direction)
    minmove_vector = np.mean(minmove_direction, axis=0)
    print('Minimum norms vector: {0}'.format(minmove_vector))
    #print(minmove_vector)
    # what is the length of the minmove_vector?
    print('Minimum norms vector length: {0}'.format(np.linalg.norm(minmove_vector)))
    # what is the norm of the movement
    minmove_norms = norms[argsort_min]
    #print(minmove_norms)
    print('Mean of minimum norms (Å): {0}'.format(np.mean(minmove_norms)))
    
    # visualise the results
    sorted_norms = np.sort(norms)
    #print("min norms")
    #print(sorted_norms[:n])
    #print(np.mean(sorted_norms[:n]))
    print("max norms")
    print(sorted_norms[-n:])
    #print(np.mean(sorted_norms[-n:]))
    
    # store all cgos
    all_cgos = []
    
    # colours
    if color_states:
        # what state are we looking at?
        # https://matplotlib.org/stable/tutorials/colors/colors.html
        # https://www.w3schools.com/colors/colors_picker.asp
        state = to[0]
        #print(state, state == 'd')
        if state == 'c':
            #color_red = [1., 0., 0.]
            this_color = [0.839, 0.153, 0.157] # tab red
        elif state == 'a':
            this_color = [0.122, 0.467, 0.706]  # tableau blue
        elif state == 'b':
            this_color = [1.0, 0.498, 0.055]  # tableau orange
        elif state == 'd':
            this_color = [0.173, 0.627, 0.173]  # tableau green
        else:
            this_color = color_blue
        #print(this_color)
        
    # show the rotation axis
    #this_cgo = make_cgo_list(line_points[0], line_points[1], this_color, this_color, 1.0)
    #all_cgos += this_cgo    
    
    # show the rotation vectors of selected coordinates
    # not including the minimum movers
    #n1 = len(xyz_notrot)
    index = np.where(norms >= sorted_norms[n])
    xyz_a = xyz_notrot[index]
    xyz_b = xyz_rot[index]
    for a, b in zip(xyz_a[::skip], xyz_b[::skip]):    
        #print(a, b)
        this_cgo = make_cgo_list(a, b)
        all_cgos += this_cgo
    #print(all_cgos)
    
    # create a ChimeraX bild file
    chimx_fn = '{0}{1}_rotation_vectors_axis.bild'.format(fr[:2],to[:2])
    chimx_f = open(chimx_fn, 'w')
    chimx_rad1 = 0.5
    chimx_f.write('.comment {0} {1}\n'.format(fr, to))
    chimx_f.write('.color {0} {1} {2}\n'.format(*this_color))
    # output this to a .bild file for chimeraX too
    for a, b in zip(xyz_a[::skip], xyz_b[::skip]):
        # print(a, b)
        astring = '.cylinder {0} {1} {2} {3} {4} {5} {6}\n'.format(*a, *b, chimx_rad1)
        chimx_f.write(astring)
    #chimx_f.close()
        
        
    # add the COM to the cgos
    all_cgos += make_cgo_list(com1, com2, color_red, color_red, 0.6)
    
    # show the movement of the n atoms that move the least
    #index = np.where(norms < sorted_norms[n])
    #xyz_a = xyz_notrot[index]
    #xyz_b = xyz_rot[index]
    #for a, b in zip(xyz_a, xyz_b):    
    #    #print(a, b)
    #    this_cgo = make_cgo_list(a, b, color_blue, color_blue, 0.6)
    #    all_cgos += this_cgo
    
    # showing the movement of the atoms which move the most
    #index = np.where(norms > sorted_norms[-n])
    #xyz_a = xyz_notrot[index]
    #xyz_b = xyz_rot[index]
    #for a, b in zip(xyz_a, xyz_b):    
    #    #print(a, b)
    #    this_cgo = make_cgo_list(a, b, color_green, color_green, 0.6)
    #    all_cgos += this_cgo
    
    # load all the graphics objects
    cmd.load_cgo(all_cgos, 'rv_{0}_{1}'.format(fr[:2], to[:2]))
    
    # print out some stuff to console
    print('Molecule COM movement = ' + '{0}'.format(com_move))
    #print('rot = ' + '{0}'.format(rot))
    #print('rot_inv = ' + '{0}'.format(rot_inv))
    
    #out_string = ''
    #out_string += fr + ','
    #out_string += to + ','
    #out_string += '{0},'.format(com_move)
    #out_string += '{0},'.format(rot)
    #out_string += '{0}\n'.format(rot_inv)
    
    # display the rotation vector compared to the minimum movers vector
    v = rot_axis
    
    def rotate(x):
        point_i = np.array([x[0], x[1], 0., 1.])
        point_f = np.matmul(np_rot, point_i)
        length = np.linalg.norm(point_f-point_i)
        #print(length)
        return length
    result = minimize(rotate, np.array((0, 0)))
    #print(result)
    #print(result.x)
    point = np.array([result.x[0], result.x[1], 0])
    print('Point of rotation axis passing through XY plane')
    print(point)
    #print(point.shape)
    #point = point.reshape((3, 1))
    #ax.scatter(*point, s=100, c='green')
    p0 = point
    
    # we want a new t to give a certain size 
    t = 100
    xyz = p0 + t*v
    p1 = p0 + -t*v
    
    # save these vectors for comparison later...
    fn = 'out_rotation_vectors.json'
    if not os.path.exists(fn):
        vecs = {}
    else:
        f = open(fn)
        s = f.read()
        f.close()
        if len(s):
            vecs = json.loads(s)
        else:
            vecs = {}
    f = open(fn, 'w')
    this_vec = '{0}_{1}'.format(fr[:2], to[:2])
    vecs[this_vec] = [p0.tolist(), v.tolist()]
    print('Output vectors to file: {0}'.format(vecs))
    f.write(json.dumps(vecs))
    f.close()
    
    # save some things to a file
    # the COM movers, and rotation/inverse matrices
    f = open('object_matrices.csv', 'w')  # append results to file
    #f.write(out_string)
    json_string = json.dumps([fr, to, com_move, rot, rot_inv])
    f.write(json_string + '\n')
    f.close()
    
    radius = 3.0
    #obj = [CYLINDER, p0[0], p0[1], p0[2], xyz[0], xyz[1], xyz[2], radius, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0]
    obj = [9.0, p1[0], p1[1], p1[2], xyz[0], xyz[1], xyz[2], radius, *this_color, *this_color]
    new_cgo = 'rv2_{0}_{1}'.format(fr[:2], to[:2])
    cmd.load_cgo(obj, new_cgo)
    
    astring = '.cylinder {0} {1} {2} {3} {4} {5} {6}\n'.format(p1[0], p1[1], p1[2], xyz[0], xyz[1], xyz[2], 1.5)
    chimx_f.write(astring)
    chimx_f.close()
    
    # save some things to a file
    # the COM movers, and rotation/inverse matrices
    f = open('object_matrices.csv', 'w')  # append results to file
    #f.write(out_string)
    json_string = json.dumps([fr, to, com_move, rot, rot_inv])
    f.write(json_string + '\n')
    f.close()
    # write initial and final coordinates to file
    f = open('init_final_points.txt', 'w')
    json_string = json.dumps([xyz_notrot.tolist(), xyz_rot.tolist()])
    f.write(json_string)
    f.close()


