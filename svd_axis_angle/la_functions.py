import numpy as np
import math as m


def Rx(theta):
    return np.matrix([[1, 0, 0],
                      [0, m.cos(theta), -m.sin(theta)],
                      [0, m.sin(theta), m.cos(theta)]])


def Ry(theta):
    return np.matrix([[m.cos(theta), 0, m.sin(theta)],
                      [0, 1, 0],
                      [-m.sin(theta), 0, m.cos(theta)]])


def Rz(theta):
    return np.matrix([[m.cos(theta), -m.sin(theta), 0],
                      [m.sin(theta), m.cos(theta), 0],
                      [0, 0, 1]])


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.

    https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
    https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def rotate_xy_plane(axis, translate=None, plane_width=2.0, step=0.2):
    if translate is None:
        translate = [0., 0., 0.]
    axis = np.array(axis)
    # make it a unit vector
    axis = axis / np.linalg.norm(axis)

    # start with points in the xy plane and rotate them to the above normal axis
    # create a meshgrid and translate these values with the transformation matrix
    x = np.arange(-plane_width, plane_width + step, step)
    y = np.arange(-plane_width, plane_width + step, step)
    xx, yy = np.meshgrid(x, y, sparse=False)
    zz = np.zeros(xx.shape)
    a = np.array([xx.flatten(), yy.flatten(), zz.flatten()])

    # we can do this various ways, but for now, will use Euler angles
    # we will rotate the axis vector to match the plane vector
    # make a rotation matrix then apply its inverse to the xy plane coords

    # what is the angle between yz and z=0
    theta = np.arctan(axis[1] / axis[2])
    # so I know I want to rotate around the x-axis by this amount
    phi1 = theta  # save this angle for later
    R = Rz(0) * Ry(0) * Rx(phi1)
    rot_Rzyz = np.matmul(R, axis)

    # now rotate the unit vector to align with the z-axis
    # we need to find the arctan of x/z
    theta = np.arctan(rot_Rzyz[0, 0] / rot_Rzyz[0, 2])
    # we need to move in the opposite direction
    theta = -theta
    theta1 = theta
    R = Rz(0) * Ry(theta1) * Rx(phi1)
    rot_Rzyz = np.matmul(R, axis)

    # The rotation matrix
    R = Rz(0) * Ry(theta1) * Rx(phi1)

    # now we have the rotation matrix, R
    # we want to apply the inverse of this to go from [0, 0, 1] to the given axis vector [a, b, c]
    # the inverse rotation matrix
    rinv = np.linalg.inv(R)
    a = np.matmul(rinv, a)

    # do any translation that is required
    for i, t in enumerate(translate):
        a[i] += t
    return a


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return np.array(vector) / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2':: """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def project_point(normal, plane_point, point):
    """
    project a point onto a plane given by
    :param normal: the normal vector to the plane you want to project onto
    :param plane_point: some point that lies on the plane
    :param point: the point you want projected onto the plane
    :return: the projected point on the plane
    """
    normal = np.array(normal, dtype=np.float)
    # the normal should be length 1
    length = np.linalg.norm(normal)
    normal /= length
    plane_point = np.array(plane_point, dtype=np.float)
    point = np.array(point, dtype=np.float)
    t = (np.dot(normal, plane_point) - np.dot(normal, point)) / np.dot(normal, normal)
    proj = point + t*normal
    return proj


def fit_line_3d(points, spread=3):
    """
    given a set of points, do SVD to
    :param points:
    :param spread:
    :return: a list with
    [0] the direction of the line in 3d
    [1] the centre of mass of the points
    [2] a bunch of points which lie on the line, with extent determined by 'spread'
    """
    # calculate the mean of the points, i.e. the 'center' of the cloud
    datamean = points.mean(axis=0)
    # do SVD on the mean-centered data to get the direction of the line
    uu, dd, vv = np.linalg.svd(points - datamean)
    # Now vv[0] contains the first principal component, i.e. the direction
    # vector of the 'best fit' line in the least squares sense.
    #print(vv[0])
    # Now generate some points along this best fit line, for plotting.
    spread = int(np.max(np.abs(points - datamean))) * spread
    linepts = vv[0] * np.mgrid[-spread:spread:2j][:, np.newaxis]
    linepts += datamean
    return vv[0], datamean, linepts


def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    #import numpy as np
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1))  # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1],
                                                                                                   points.shape[0])
    ctr = points.mean(axis=1)  # the centre of the point cloud
    x = points - ctr[:, np.newaxis]  # shift the points to the centre of mass
    a = np.dot(x, x.T)  # Could also use np.cov(x) here.
    #print('planeFit dot product')
    #print(a)
    svd_result = svd(a)
    #print('planeFit svd result')
    #for n, r in enumerate(svd_result):
    #    print(r)
    #print(svd_result)
    svd_0 = svd(a)[0][:, -1]
    #print(svd_0)
    return ctr, svd_0


def planeFit2(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    from numpy.linalg import svd
    #print(points.shape)
    #print('reshape')
    points = np.reshape(points, (np.shape(points)[0], -1))  # Collapse trailing dimensions
    #print(points.shape)
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1],
                                                                                                   points.shape[0])
    ctr = points.mean(axis=1)  # the centre of the point cloud
    #print(ctr, np.linalg.norm(ctr))
    if np.linalg.norm(ctr) > 0.01:
        x = points - ctr[:, np.newaxis]  # shift the points to the centre of mass
    else:
        x = points
    a = np.dot(x, x.T)  # Could also use np.cov(x) here.
    #print('planeFit dot product')
    #print(a)
    svd_result = svd(a)
    #print('planeFit svd result')
    #for n, r in enumerate(svd_result):
    #    print(r)
    #print(svd_result)
    svd_0 = svd(a)[0]
    #print(svd_0)
    #print()
    svd_u = [svd_0[:, 0], svd_0[:, 1], svd_0[:, 2]]
    return ctr, svd_u


# this code taken from transforms3d package
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

    Notes
    -----
    http://en.wikipedia.org/wiki/Rotation_matrix#Axis_of_a_rotation
    """
    print(mat.dtype)
    # M = np.asarray(mat, dtype=np.float)
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
        sina = (M[1, 0] + (cosa - 1.0) * direction[0] * direction[1]) / direction[2]
    elif abs(direction[1]) > 1e-8:
        sina = (M[0, 2] + (cosa - 1.0) * direction[0] * direction[2]) / direction[1]
    else:
        sina = (M[2, 1] + (cosa - 1.0) * direction[1] * direction[2]) / direction[0]
    angle = m.atan2(sina, cosa)
    return direction, angle
