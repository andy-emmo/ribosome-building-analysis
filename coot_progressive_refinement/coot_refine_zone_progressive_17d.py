from datetime import datetime
import sys
import gc  # garbage collector
import random
from subprocess import Popen, PIPE
import os
import pickle

# coot directory taken from ccpem-1.4.1/lib/python2.7/site-packages/coot/

# progressive fit of structure to electron density map
# multiple chains
# determines free memory (due to coot memory leak during refinement) and quits on low memory
# refines from geometric centre (TODO: refine from best fitted regions)
# autosave pdb every X number of refinements
# quits when a refinement exceeds 5.0 Angstrom
# checks for refined residues/nucleotides and marks as refined after movement less than X
# stops refinement from occurring in low density regions de novo, wait for chain to be fitted in this region
# progressive refinement from geometric centre


def get_euclidean_distance(a, b):
    print('calculating euclidean distance')
    print(a, b)
    dist = ((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)**0.5
    return dist


class ChainRefiner(object):
    def __init__(self, parent=None):
        """
        Parameters
        start_at (int): which residue number to start at, default 0.
        restart_every (int): number of refinements to do before going back to the COM, consider 10-30.
        stringency (float): if a residue moves less than this, consider it refined, consider 0.001-0.005.
        set_matrix(300): default is 60, use estimate refinement weight tool
        chains (list): ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] chains in this coot instance to fit
        density_cutoff (float): map density to consider a residue to fit, consider 0.035-0.05
          (0.05 for first fit, then 0.04 for outer stretches)
        save_pdb_every (int): number of refinements before saving the pdb, consider 5-15
        """

        # format the date for saving PDBs
        self.datestring = '%Y%m%d-%H%M%S'
        # a dict to store all the residues
        self.models = []
        self.chain_ids = []
        self.all_residues = {}
        self.which = 0
        self.total_residue_count = 0

        #set_matrix(1200.00)
        
        # a few coot settings
        # turn off smooth scrolling
        #set_smooth_scroll_flag(0)
        # turn on immediate replacement
        #self.replace_state = refinement_immediate_replacement_state()
        # set rigid body fraction of atoms in positive density
        #set_rigid_body_fit_acceptable_fit_fraction(0.2)
        #set_refine_with_torsion_restraints(0)
        #remove_planar_peptide_restraints()
        
        #self.refine_residue_sphere(self.all_residues[754])
        #self.do_refinement()

    def set_parameters(self):
        self.refinement_rounds = 50000
        if self.which == 1:
            # there are around ~30 residues in a 12 A radius sphere
            self.refinement_rounds = 50000
            self.radius = 12.0
            self.decimate = 12
            self.save_every = 15
        elif self.which == 2:
            # there are around ~1500 residues in a 50 A radius sphere
            self.radius = 50.0  # same as current key binding 'Y'
            self.decimate = 100  # TODO: might need a better way to decimate if getting rid of this many residues
            self.save_every = 2
        elif self.which == 3:
            self.refinement_rounds = 50000
            self.radius = 8.0
            self.decimate = 6
            self.save_every = 30
        elif self.which == 4:
            # there are around ~30 residues in a 12 A radius sphere
            self.refinement_rounds = 50000
            self.radius = 3.5
            self.decimate = 1
            self.save_every = 200
        elif self.which == 5:
            # there are around ~400 residues in a 30 A radius sphere
            self.refinement_rounds = 500
            self.radius = 30.0
            self.decimate = 50
            self.save_every = 5
        else:
            # there are around ~250 residues in a 25 A radius sphere
            self.radius = 25.0
            self.decimate = 100 
            self.save_every = 5

    def populate_residue_list(self):
        if os.path.exists('all_residues.pkl'):
            f = open('all_residues.pkl', 'rb')
            self.all_residues = pickle.load(f)
            f.close()
        else:
            self.all_residues = []
            # get all molecules
            self.models = model_molecule_list()  # returns a list of molecules
            print('models')
            print(self.models)
            self.chain_ids = []
            for imol in self.models:
                # get the centre of this molecule
                mc = molecule_centre(imol)
                print(mc)
                chains = run_scheme_command("(chain-ids {0})".format(imol))
                print(chains)
                self.chain_ids.append(chains)
                for chain in chains:
                    ric = residues_in_chain(imol, chain)
                    print('for chain in chains')
                    print(ric)
                    print(imol)
                    print(chain)
                    print(isinstance(ric, bool))
                    if not isinstance(ric, bool):  # debug: some structures return a False on residues_in_chain
                        for res in ric:
                            #res_info = residue_info(imol, *res)
                            #rc = residue_centre(imol, *res)  # this gives errors in COOT ccpem 0.9.2
                            rc = residue_centre(imol, *res[1:])
                            print('residue_centre')
                            print(rc)
                            try:
                                dist_from_centre = get_euclidean_distance(rc, mc)
                                refined = 0
                                this_res = [imol, res, rc, dist_from_centre, refined]
                                print(this_res)
                                self.all_residues.append(this_res)
                            except TypeError:
                                pass
                        #print(self.all_residues)
            # sort residues according to their distance from centre
            self.all_residues = sorted(self.all_residues, key=lambda x: x[3])
            # decimate all_residues list
            # i.e. only do every X number of residues
            self.total_residue_count = len(self.all_residues)
            self.all_residues = self.all_residues[::self.decimate]
            self.save_all_residues()

    def save_all_residues(self):
        # save all_residues to a file
        f = open('all_residues.pkl', 'wb')
        pickle.dump(self.all_residues, f)
        f.close()

    def refine_residue_sphere(self, this_res):
        imol = this_res[0]
        res = this_res[1]
        res_near = residues_near_residue(imol, res, self.radius)  # int:imol, list:residue, float:radius
        res_near = [res] + res_near
        print('There are {1} residues within {0} angstroms.'.format(self.radius, len(res_near)))
        # centre of the residue
        rc = residue_centre(imol, *res[1:])  # returns a list ['residue centre', [float, float, float]
        print('The centre of the current residue is {0}'.format(rc))
        # move the screen to the current residue
        # set_rotation_centre(*residue_centre(0, self.current_chain, self.res_diffs[self.current_residue][0], ""))
        try:
            set_rotation_centre(*rc)
        except TypeError as t:
            print(t)
        
        set_refinement_immediate_replacement(1)
        #refine_residues(imol, res_near)
        # SUCCESS!! But it misses the central residue
        #scm = "(refine-residues {0} (residues-near-residue {0} (list \"{2}\" {3} \"\") {1:0.1f}))".format(imol, self.radius, res[0], res[1])
        #print(scm)
        #run_scheme_command(scm)
        # Try it with a list constructor
        scm = "(refine-residues {0} (cons (list \"{2}\" {3} \"\") (residues-near-residue {0} (list \"{2}\" {3} \"\") {1:0.1f})))"
        print(res)
        res = res[1:]
        print(res)
        scm = scm.format(imol, self.radius, res[0], res[1])
        print('Executing Scheme command:')
        print(scm)
        run_scheme_command(scm)
        #morph_fit_residues(imol, res_near, 9.0)
        accept_moving_atoms()  # alias for accept_regularizement()
        set_refinement_immediate_replacement(0)
        print('RESIDUE CENTRE BEFORE AND AFTER:')
        print(rc, residue_centre(imol, *res))

    def residues_refined(self):
        # sort all_residues by their distance from centre
        self.all_residues = sorted(self.all_residues, key=lambda x: x[3])
        done = 0
        for res in self.all_residues:
            if res[4]:
                done += 1
        print('Completed {0} of {1} residues to be refined, that is {2:0.1f}%'.format(done, len(self.all_residues), 100.*done/len(self.all_residues)))
        self.datestring = '%Y%m%d-%H%M%S'
        fn = 'andy_log.txt'
        if not os.path.exists(fn):
        	f = open(fn, 'w')
        else:
        	f = open(fn, 'a')
        s = 'At {3} Completed {0} of {1} residues to be refined, that is {2:0.1f}%'
        s = s.format(done, len(self.all_residues), 100.*done/len(self.all_residues), datetime.now())
        if self.total_residue_count > 0:
        	f.write(s + ' Total residues {0}'.format(self.total_residue_count) + '\n')
        else:
        	f.write(s + '\n')
        f.close()

    def do_refinement(self):
    	"""this is the function that the coot menu item executes"""
    	self.set_parameters()
        self.populate_residue_list()
        self.residues_refined()
        
        i = 0  # how many refinements done this round
        j = 0  # the current residue to consider
        while i < self.refinement_rounds and j < len(self.all_residues):
            # if this residue hasn't been refined yet
            #print('doing refinement')
            #print(self.all_residues[j][4])
            if self.all_residues[j][4] < 1:
                #print('doing refinement')
                self.refine_residue_sphere(self.all_residues[j])
                self.all_residues[j][4] += 1
                i += 1
                j += 1
                # do a save every so often
                if not i % self.save_every and i != 0:
                    quick_save()
                    self.save_all_residues()
                    self.residues_refined()
                # check memory every round that finishes
                if self.check_memory():
                    i = self.refinement_rounds + 1
            else:
                #print('skipping, already done')
                j += 1
        quick_save()
        self.save_all_residues()
        self.residues_refined()

    def refine_chain(self):
        # file to save refined residues, this should save some time on rerunning refinements
        self.finished_fn = 'finished_{0}.txt'.format(self.current_chain)
        if os.path.exists(self.finished_fn):
            finished_file = open(self.finished_fn)
            # load file contents into finished list
            d = finished_file.read()
            d = d.split(',')
            self.finished_list = []
            for f in d:
                if len(f):
                    self.finished_list.append(int(f))
            print('loaded finished list: {0}'.format(self.finished_list))
            self.logfile.write('loaded finished list: {0}'.format(self.finished_list))
        else:
            self.finished_list = []

        # get the residues sorted by distance from COM
        # also gets residue ids
        self.do_get_coords()
        # set the offset for
        self.offset = self.restart_every

        while 1:
            if self.do_refinement_round():
                print('refinement round returned 1, exiting')
                return 1
            self.offset += self.restart_every
            if self.check_memory():
                return 1
            self.write_to_pdb()

            print('Refinements done this round: {0}'.format(self.refinements_done))
            self.logfile.write('Refinements done this round: {0}\n'.format(self.refinements_done))
            # length of the chain
            test_offset = self.offset > len(self.res_ids)
            astring = 'Offset length {1} greater than number of residues{2}: {0}'.format(test_offset, self.offset, len(self.res_ids))
            print(astring)
            astring = astring + '\n'
            self.logfile.write(astring)

            if self.refinements_done == 0 and test_offset:
                print('No refinements done this round, finishing')
                self.write_to_pdb()
                return 0
            #return

    def do_refinement_round(self):
        # do a refinement round
        self.refinements_done = 0
        print('doing refinement round')
        for i in range(self.start_at, self.start_at + self.offset):
            ## check this is right
            #self.current_residue = self.res_ids[i]
            self.current_residue = i
            print(i, len(self.res_diffs))
            if i < len(self.res_diffs):
                print(self.res_diffs[self.current_residue][0] not in self.finished_list)
                if self.res_diffs[self.current_residue][0] not in self.finished_list:
                    # if the return value is 1 some error has occurred
                    if self.do_residue_refinement():
                        print('do residue returned 1, exiting')
                        return 1
                    self.refinements_done += 1
                    # garbage collect?
                    #print(gc.get_count())
                    #print(gc.collect())
                else:
                    print('res {0} in finished list, skipping'.format(self.res_diffs[self.current_residue][0]))
                    #print(self.finished_list)
                    #print(self.res_diffs[self.current_residue][0])
                    #print(self.res_diffs[self.current_residue][0] not in self.finished_list)
                    #print(self.res_diffs[self.current_residue][0] in self.finished_list)
                    #input('...')
        # save finished list at end of refinement round
        self.save_finished_list()
        # on a successful round return 0
        print('finished refinement round')
        return 0

    def do_residue_refinement(self):
        # do the refinement of a single residue
        # determine whether or not to refine this residue
        #r = self.res_diffs[self.current_residue]
        do_this_res = 1

        try:
            if self.res_diffs[self.current_residue][0] in self.skip_list:
                do_this_res = 0
            #if self.res_diffs[self.current_residue][0] in self.finished_list:
            #    do_this_res = 0
        except Exception as e:
            print('Exception!')
            print(e)
            print(self.res_diffs)
            print(self.res_diffs[self.current_residue])
            print(self.res_diffs[self.current_residue][0])
            print(self.res_diffs[self.current_residue][0] in self.skip_list)
            print(self.res_diffs[self.current_residue][0] in self.finished_list)
            print(self.skip_list)
            print(self.finished_list)
            input('Waiting for input...')

        print('do this res {1}? {0}'.format(do_this_res, self.res_diffs[self.current_residue][0]))
        if not do_this_res:
            return 0

        # get a contiguous region around the current residue
        # do the upper stretch first
        upper = self.res_diffs[self.current_residue][0]
        for i in range(random.choice(self.length_choice)):
            this_res = self.res_diffs[self.current_residue][0] + i + 1
            if this_res in self.res_ids:
                upper = this_res
            else:  # not in list, keep last upper res
                break
        # now do the lower stretch
        lower = self.res_diffs[self.current_residue][0]
        for i in range(random.choice(self.length_choice)):
            this_res = self.res_diffs[self.current_residue][0] - i - 1
            if this_res in self.res_ids:
                lower = this_res
            else:  # not in list, keep last lower res
                break
                
        # before refinement, get the residue coords.
        centres_before = []
        for this_res in range(lower, upper + 1):
            print(this_res)
            rc = residue_centre(0, self.current_chain, this_res, '')
            print(rc)
            centres_before.append(rc)

        #this_density = []

        # before refinement, check density, if too low, don't fit yet...
        # not sure what the best way is to check this, mean density or density per res?
        check_density = 0
        density_count = 0
        density_sum = 0
        print('Checking density')
        for cb in centres_before:
            density_here = density_at_point(1, cb[0], cb[1], cb[2])
            print(density_here)
            density_count += 1
            density_sum += density_here
            # if the density is low for all centres, skip this refinement
            if density_here > self.density_cutoff:
                check_density += 1
        density_ratio = float(check_density) / float(density_count)
        print('Mean ratio:')
        print(density_ratio)
        print('Mean density:')
        print(density_sum / density_count)

        # logfile.write('Density in region too low skipping refinement {0}\n'.format(density_ratio))

        if density_ratio >= self.density_accept:
            # move the screen centre to the residue
            set_rotation_centre(*residue_centre(0, self.current_chain, self.res_diffs[self.current_residue][0], ""))
            # increment the number of refinements done
            self.refinements_done += 1
            # timestamp the log file
            #now = datetime.now()
            #now = now.strftime(self.datestring)
            #self.logfile.write(now + '\n')
            #self.logfile.write('Range of refinement {0}, {1} for res {2}\n'.format(lower, upper, r[0]))
            #print('Range of refinement {0}, {1} for res {2}'.format(lower, upper, r[0]))
            #self.logfile.write(
            #    '{0}, {1}, {2}, {3}, {4}, {5}\n\n'.format(n, r[0], start_at, offset, restart_every, finished))
            #print(n, r[0], start_at, offset, restart_every, finished)

            # do the refinement
            refine_zone(0, self.current_chain, lower, upper, "")
            accept_regularizement()

            # after refinement, get new residue coordinates
            centres_after = []
            these_reses = []
            for this_res in range(lower, upper + 1):
                # print(this_res)
                rc = residue_centre(0, self.current_chain, this_res, '')
                # print(rc)
                centres_after.append(rc)
                these_reses.append(this_res)

            # now check that the centres_before and centres_after
            # don't shift by an unreasonable amount
            for cb, ca, rid in zip(centres_before, centres_after, these_reses):
                x = (cb[0] - ca[0]) ** 2
                y = (cb[1] - ca[1]) ** 2
                z = (cb[2] - ca[2]) ** 2
                dist = (x + y + z) ** 0.5
                print('dDist: {0:0.5f} for resid {1}'.format(dist, rid))
                dapb = density_at_point(1, cb[0], cb[1], cb[2])
                dapa = density_at_point(1, ca[0], ca[1], ca[2])
                print(
                    'density at point before/after: {0:0.05f} {1:0.05f} {2:0.05f}'.format(dapb, dapa, dapa - dapb))

                if dist <= self.stringency:
                    if rid not in self.finished_list:
                        self.finished_list.append(rid)
                    print('Residue {0} finished refinement'.format(rid))

                if dist > 5.0:
                    print('Distance of refinement move exceeded 5.0 Angstrom, exiting refinement')
                    fn = 'error.txt'
                    if not os.path.exists(fn):
                        f2 = open(fn, 'w')
                    else:
                        f2 = open(fn, 'a')
                    f2.write('Distance of refinement move exceeded 5.0 Angstrom, exiting refinement\n')
                    f2.write(' This occurred at {0}'.format(datetime.now().strftime(self.datestring)))
                    f2.write(' for residue id:{0}\n'.format(rid))
                    f2.close()
                    return 1
        return 0

    def do_get_coords(self):
        # get all the residues in the chain we're working on

        self.res_ids = []
        self.res_info = []
        for r in ric:
            self.res_ids.append(r[1])
            self.res_info.append(r)

        # list to store the information
        self.res_coords = []
        # sum all the coords to get the naive COM
        x = 0
        y = 0
        z = 0
        res_count = 0
        for r in ric:
            print(r)
            rc = residue_centre(0, *r)
            # print(rc)
            resno = r[1]
            self.res_coords.append([resno, rc[0], rc[1], rc[2]])
            x += rc[0]
            y += rc[1]
            z += rc[2]
            res_count += 1
        # naive com of residues
        x /= float(res_count)
        y /= res_count
        z /= res_count
        #print('Residue COM: {0} {1} {2}'.format(x, y, z))

        # save the residue coordinates
        f = open('res_coords.txt', 'w')
        for r in self.res_coords:
            f.write('{0},{1},{2},{3}\n'.format(*r))
        f.close()

        # get the euclidean distance between the COM and all residues
        self.res_diffs = []
        for r in self.res_coords:
            dx = (x - r[1]) ** 2
            dy = (y - r[2]) ** 2
            dz = (z - r[3]) ** 2
            dist = (dx + dy + dz) ** 0.5
            self.res_diffs.append([r[0], r[1], r[2], r[3], dist])
        # sort the res_diffs according to diff
        self.res_diffs = sorted(self.res_diffs, key=lambda x: x[4])

    def do_get_density(self):
        # get the density of all residues in self.res_id
        self.do_get_coords()
        # now the residues and their x,y,z coords are stored in self.res_coords
        density_count = 0
        density_sum = 0
        print('Checking density')
        self.res_density = []
        for ri in self.res_info:
            # get all the atom information for this residues
            print(ri)
            imol = 0
            print(imol, ri[0], ri[1], ri[2])
            this_res_info = residue_info(imol, ri[0], ri[1], ri[2])
            density_for_residue = []
            # for each atom in this residue
            atom_posns = []
            atom_densities = []
            atom_density_sum = 0
            atom_count = len(this_res_info)
            print(atom_count)
            for this_info in this_res_info:
                print(this_info)
                map_imol = 1
                c = this_info[2]
                print(c[0], c[1], c[2])
                atom_density = density_at_point(map_imol, c[0], c[1], c[2])
                atom_densities.append(atom_densities)
                atom_density_sum += atom_density
            print('Total density for this residues atoms:')
            print(atom_density_sum)
            print('Average denisty for this residues atoms:')
            av_density = atom_density_sum/atom_count
            print(av_density)
            self.res_density.append([ri[1], av_density])
            # changing a thing
            # get density using the residue COM
            #density_here = density_at_point(1, ri[1], ri[2], ri[3])
            #print(density_here)
            #density_count += 1
            #density_sum += density_here
            #self.res_density.append([ri[0], density_here])

        #print('Mean density:')
        #print(density_sum / density_count)

        # save the residue density
        f = open('res_density.txt', 'w')
        for r in self.res_density:
            f.write('{0},{1}\n'.format(*r))
        f.close()
        return

    def do_get_relative_density(self):
        # get the density of all residues in self.res_id
        self.do_get_coords()
        self.create_coord_list()
        # now the residues and their x,y,z coords are stored in self.res_coords
        density_count = 0
        density_sum = 0
        print('Checking density')
        self.res_density = []
        for ri in self.res_info:
            # get all the atom information for this residues
            print(ri)
            imol = 0
            print(imol, ri[0], ri[1], ri[2])
            this_res_info = residue_info(imol, ri[0], ri[1], ri[2])
            density_for_residue = []
            # for each atom in this residue
            atom_posns = []
            atom_densities = []
            atom_density_sum = 0
            atom_count = len(this_res_info)
            print(atom_count)
            for this_info in this_res_info:
                print(this_info)
                map_imol = 1
                c = this_info[2]
                print(c[0], c[1], c[2])
                atom_density = density_at_point(map_imol, c[0], c[1], c[2])
                atom_densities.append(atom_densities)
                atom_density_sum += atom_density
            print('Total density for this residues atoms:')
            print(atom_density_sum)
            print('Average denisty for this residues atoms:')
            av_density = atom_density_sum/atom_count
            print(av_density)

            # changing a thing
            # get density using the residue COM
            #density_here = density_at_point(1, ri[1], ri[2], ri[3])
            # using the sphere coordinates, get the densities around this residue
            max_density = 0
            rc = residue_centre(0, *ri)
            for x in self.sphere_coord_list:
                print(rc[0], x[0], rc[0]+x[0])
                density_here = density_at_point(1, rc[0]+x[0], rc[1]+x[1], rc[2]+x[2])
                print(density_here)
                if density_here > max_density:
                    max_density = density_here
            print(max_density)
            self.res_density.append([ri[1], av_density, max_density])
            #print(density_here)
            #density_count += 1
            #density_sum += density_here
            #self.res_density.append([ri[0], density_here])

        #print('Mean density:')
        #print(density_sum / density_count)

        # save the residue density
        f = open('res_density.txt', 'w')
        print(self.res_density)
        for r in self.res_density:
            f.write('{0},{1},{2}\n'.format(*r))
        f.close()
        return

    def create_coord_list(self, radius=3.0, granular=4):
        def euc_distance(i, j, k):
            return (i ** 2 + j ** 2 + k ** 2) ** 0.5
        # create a number of coordinates within a sphere
        # in x, y, z
        step = radius / granular
        x = []
        # create x coords
        for r in range(-granular, 0):
            x.append(step * r)
        for r in range(0, granular + 1):
            x.append(step * r)
        self.sphere_coord_list = []
        for i in x:
            for j in x:
                for k in x:
                    dist = euc_distance(i, j, k)
                    if dist <= radius:
                        self.sphere_coord_list.append([i, j, k])

    def check_memory(self):
        p = Popen(['cat', '/proc/meminfo'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, error = p.communicate()
        # print(out)
        # print(mem)
        # print(type(mem))
        out = out.split('\n')
        check_this = 'MemAvailable'
        #check_this = 'MemFree'
        if self.radius > 20:
            mem_limit = 2048000
        elif self.radius > 10:
            mem_limit = 1024000
        else:
            mem_limit = 512000
        for o in out:
            print(o)
            if check_this in o:
                # print('found')
                mem = ''
                for i in o:
                    # print(i)
                    # print(i.isdigit())
                    if i.isdigit():
                        mem += i
        mem = int(mem)
        print('Free memory is {0} kb'.format(mem))
        if mem < mem_limit:
            print('Memory low, quitting refinement')
            return 1
        else:
            return 0

# functions to link to coot menu
def do_refine_from_centre(arg):
    print(arg)
    cr = ChainRefiner(arg)
    cr.do_refinement()
    return
    
def do_refine_from_centre_8(arg):
    #print('the_arg', arg)  # the arg is the menu item
    cr = ChainRefiner(arg)
    cr.which = 3
    cr.do_refinement()
    return
    
def do_refine_from_centre_25(arg):
    #print('the_arg', arg)  # the arg is the menu item
    cr = ChainRefiner(arg)
    cr.which = 0
    cr.do_refinement()
    return

def output_coords(arg):
    print(arg)
    cr = ChainRefiner(config='rna-outer')
    cr.do_get_coords()
    return

def output_density(arg):
    print(arg)
    cr = ChainRefiner(config='rna-outer')
    cr.do_get_density()
    return

def output_relative_density(arg):
    print(arg)
    cr = ChainRefiner(config='rna-outer')
    cr.do_get_relative_density()
    return

def reset_files(arg):
	os.remove('andy_log.txt')
	os.remove('all_residues.pkl')
	return

#now = datetime.now()
menu_object = coot_menubar_menu('andy')
add_simple_coot_menu_menuitem(menu_object, 'reset_files', reset_files)
add_simple_coot_menu_menuitem(menu_object, 'rfc_{0}'.format('8'), do_refine_from_centre_8)
add_simple_coot_menu_menuitem(menu_object, 'rfc_{0}'.format('25'), do_refine_from_centre_25)

#add_simple_coot_menu_menuitem(menu_object, 'output-coords', output_coords)
#add_simple_coot_menu_menuitem(menu_object, 'output-density', output_density)
#add_simple_coot_menu_menuitem(menu_object, 'output-relative-density', output_relative_density)

