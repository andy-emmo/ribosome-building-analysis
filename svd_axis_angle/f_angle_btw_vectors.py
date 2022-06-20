import numpy as np


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return np.array(vector) / np.linalg.norm(vector)


def angle_between(v1, v2, degrees=True):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    v1 = np.array(v1)
    v2 = np.array(v2)
    #print(v1.shape)
    if v1.shape[0] == 6:
        v1[3] -= v1[0]
        v1[4] -= v1[1]
        v1[5] -= v1[2]
        v1 = v1[3:]
        print(v1)
        print(v1.shape)
    if v2.shape[0] == 6:
        v2[3] -= v2[0]
        v2[4] -= v2[1]
        v2[5] -= v2[2]
        v2 = v2[3:]
        print(v2.shape)

    v1_u = unit_vector(v1)
    print(v1_u)
    v2_u = unit_vector(v2)
    print(v2_u)

    if degrees:
        return np.rad2deg(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    else:
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


if __name__ == '__main__':
    # test
    vector_1 = [0, 1]
    vector_2 = [1, 0]
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle = np.arccos(dot_product)
    print(angle)

    print(angle_between((1, 0, 0), (-1, 0, 0)))
    print(angle_between([0, 0, 1], [0.567, 1.87162, 13.9821]))

    # c0d0 compared to Z-axis
    print(angle_between([0, 0, 1], [4.803964331214903, 5.675695710353089, 98.5049652555547]))
    # a0c0 compared to Z-axis
    print(angle_between([0, 0, 1], [22.762659731605723, -14.564094511127209, 73.60151726272912]))
    # c0d0 compared to a0c0
    print(angle_between([4.803964331214903, 5.675695710353089, 98.5049652555547],
                        [22.762659731605723, -14.564094511127209, 73.60151726272912]))

    # these point/vector representations of the rotation axis need to be moved to the origin
    # taken from ChimeraX visualisation:
    # c0d0 .cylinder -2.7619520384163514 39.2888507169215 -98.5049652555547 4.803964331214903 5.675695710353089 98.5049652555547 1.5
    # a0c0 .cylinder -30.124166659441023 110.07237776236114 -73.60151726272912 22.762659731605723 -14.564094511127209 73.60151726272912 1.5
    # here the two-endpoint form can be interpreted as [offset_x, offset_Y, offset_Z, vector_x, vector_y, vector_Z]

    print('\n SSU head angle between (not aligned).')
    c0d0 = [-2.7619520384163514, 39.2888507169215, -98.5049652555547, 4.803964331214903, 5.675695710353089, 98.5049652555547]
    d0b0 = [-12.562503918952434, 48.52230303313735, -96.53155761531234, 4.22093170703418, -0.9241308730091191, 96.53155761531234]
    b0a0 = [-9.603649822765469, 88.61706871431915, -85.66680402795741, 20.481064590613197, -10.072873575701472, 85.66680402795741]
    a0c0 = [-30.124166659441023, 110.07237776236114, -73.60151726272912, 22.762659731605723, -14.564094511127209, 73.60151726272912]
    print('SSU body rot axis angle between c0d0 a0c0')
    print(angle_between(c0d0, a0c0))
    print('SSU body rot axis angle between c0d0 d0b0')
    print(angle_between(c0d0, d0b0))
    print('SSU body rot axis angle between c0d0 b0a0')
    print(angle_between(c0d0, b0a0))
    print('SSU body rot axis angle between d0b0 b0a0')
    print(angle_between(d0b0, b0a0))

    print('\n SSU head angle between (not aligned).')
    # for the head movement SSU body aligned
    c0d0 = [-17.211922358974633, 49.54004376725196, -111.32846446289183, 13.698760485134379, 5.63289830097904, 81.32846446289183]
    d0b0 = [-94.17250208255784, 86.09464817285912, 11.288889151065646, 84.62950924559719, 13.531448313576064, -41.288889151065646]
    b0a0 = [86.1619429179321, 5.0392817906824945, 30.56110243243934, -75.91727058865874, 78.70794765106052, -60.56110243243934]
    a0c0 = [-70.75138781963443, 108.37772469011429, -72.1120391804277, 63.09104605105845, 13.303738667741335, 42.11203918042771]
    print('SSU head rot axis angle between c0d0 a0c0')
    print(angle_between(c0d0, a0c0))
    print('SSU head rot axis angle between c0d0 d0b0')
    print(angle_between(c0d0, d0b0))
    print('SSU head rot axis angle between c0d0 b0a0')
    print(angle_between(c0d0, b0a0))
    print('SSU head rot axis angle between d0b0 b0a0')
    print(angle_between(d0b0, b0a0))

    print('\n SSU head angle between body aligned.')
    # for the head movement SSU body aligned
    c0d0 = [-89.18121083321103, 86.36756917503982, -78.51773098943693, 48.14996658314047, 15.64176688724185, 48.51773098943693]
    d0b0 = [-81.52776831077195, 77.70498045986835, -49.169018997952556, 90.5969144810657, 2.186730104874684, 19.169018997952556]
    b0a0 = [-83.67416333051558, 70.46623841032192, -41.7137919256817, 98.95213012847374, 8.877844850848028, 11.713791925681697]
    a0c0 = [2.9582202605652093, 112.88679292545626, -90.88258535175106, -31.5769913911754, -12.710389039372629, 60.88258535175106]
    print('SSU head rot axis angle between c0d0 a0c0')
    print(angle_between(c0d0, a0c0))
    print('SSU head rot axis angle between c0d0 d0b0')
    print(angle_between(c0d0, d0b0))
    print('SSU head rot axis angle between c0d0 b0a0')
    print(angle_between(c0d0, b0a0))
    print('SSU head rot axis angle between d0b0 b0a0')
    print(angle_between(d0b0, b0a0))

