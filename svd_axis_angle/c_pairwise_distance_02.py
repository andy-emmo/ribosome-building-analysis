from pdb_functions import read_pdb
import os
import numpy as np
import matplotlib.pyplot as plt


def pairwise_euclidean(matrices, is_protein=False, comparison='all', thresh=1.0, invert=True):
    #print(matrices)
    if len(matrices) == 1:
        # making a matrix against self
        a = np.array(matrices[0])
        b = np.array(matrices[0])
    else:
        a = np.array(matrices[0])
        b = np.array(matrices[1])
    #print(a.shape, b.shape)
    # create a matrix to store the results of the Euclidean distances
    cs = np.zeros((a.shape[0], b.shape[0], a.shape[1]))
    #print(cs.shape)
    #input()
    # matrix a will be the rows of the result matrix c
    # for every row of a, we compare against every column of residue from matrix b
    for b0 in range(b.shape[0]):
        #print(b0)
        this_a = a[:, :]
        this_b = b[b0, :]
        #print(this_a.shape)
        #print(this_b.shape)
        #res = np.linalg.norm(this_a - this_b, axis=2)
        #print(res.shape)
        #print(np.linalg.norm(this_a - this_b, axis=1))
        cs[:, b0, :] = np.linalg.norm(this_a - this_b, axis=2)
        #print(cs[:, b0, :].shape)
        #input()
    # we will compare the mean of all distances
    if comparison == 'all':
        c = np.mean(cs, axis=2)
    else:
        # compare just the ribose ring centres
        c = cs[:, :, 0]
    if invert:
        # for any values of zero, inverting will result in np.inf
        index = np.where(c < 0.001)
        c[index] = 0.001
        a = 1 / c
        index = np.where(a > thresh)
        a[index] = thresh
    else:
        a = c
    return a
    

def distance_between(pdb1, pdb2):
    fns = [pdb1, pdb2]
    all_coords = {}
    labels = []
    # get the coordinates
    for i, fn in enumerate(fns):
        print(i, fn)
        label = fn.split('_')[0]
        coords, this_seq, is_protein = read_pdb(pdb_fn=fn, gap_char='x')
        # coords returns
        # [ribose_ring, base_elements, O5', O3', C5']
        coords = np.array(coords)
        #print(coords.shape)
        #print(coords[66])
        coords = np.mean(coords, axis=1)
        #print(coords.shape)
        #print(coords[66])
        all_coords[label] = coords
        labels.append(label)

    # determine the euclidean norms
    all_norms = []
    all_norms_labels = []
    print(len(all_coords))
    print(labels)
    compare_to = 0
    for i in range(0, len(all_coords)):
        if i != compare_to:
            print(i, labels[i], labels[compare_to])
            this_label = '{0} - {1}'.format(labels[i].replace('gi', ''), labels[compare_to].replace('gi', ''))
            all_norms_labels.append(this_label)
            a = all_coords[labels[compare_to]]
            b = all_coords[labels[i]]
            ba = b-a
            #print(a[66], b[66], ba[66])
            this_norm = np.linalg.norm(ba, axis=1)
            #print(this_norm[66])
            #input()
            #print(this_norm.shape)
            all_norms.append(this_norm)

    all_norms = np.array(all_norms)
    #print(all_norms.shape)
    return all_norms
    
    
def plot_distances(all_norms, all_norms_labels):
    # plot the euclidean norms
    alpha = 1.0
    for n, an in enumerate(all_norms):
        this_label = all_norms_labels[n]
        print(this_label, this_label[0], this_label[0].lower())
        c = 'C8'
        if this_label[0].lower() == 'a':
            c = 'C0'
        elif this_label[0].lower() == 'b':
            c = 'C1'
        elif this_label[0].lower() == 'c':
            c = 'C3'
        elif this_label[0].lower() == 'd':
            c = 'C2'
        x = np.arange(1, len(an)+1)
        print(x[:10])
        plt.plot(x, an, label=all_norms_labels[n], color=c, alpha=alpha)

    print(plt.ylim())
    plt.ylim((-1.24838045618071, 27.4932030992973))
    print(plt.ylim())
    #plt.legend()
    #ax = plt.subplot(111)
    #plt.legend(bbox_to_anchor=(1.0, 1.0))  # bbox_to_anchor (lr, tb)
    plt.tight_layout(rect=(0.05, 0.05, 0.95, 0.95))  # (left, bottom, right, top)
    plt.xlabel('rRNA residue number')
    plt.ylabel('Euclidean norm (Å)')
    plt.axvline(895, color='k', linewidth=0.8)
    plt.axvline(1315, color='k', linewidth=0.8)
    plt.axvspan(895, 1315, alpha=0.1, color='k')
    plt.gcf().set_size_inches(6, 4)
    plt.savefig('c_pairwise_eucnorm.png', dpi=150)
    #plt.savefig('pairwise_eucnorm.svg')
    plt.show()
    

if __name__ == '__main__':
    distances = distance_between('./structures/c0_svd.pdb', './structures/d2_svd.pdb')
    plot_distances(distances, ['C0', 'D2'])
    
