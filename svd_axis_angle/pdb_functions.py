import numpy as np


aa_one_to_three = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly',
                   'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser',
                   'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}
aa_three_to_one = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                   'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                   'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def generate_alphanum():
    alphanum = []
    for i in range(128):
        a = chr(i)
        #print(a)
        if a.isalnum():
            alphanum.append(a)
    #print(alphanum)
    return alphanum


def convert_pdb_string_to_number(astring):
    atom_number = 0
    for n, i in enumerate(reversed(astring)):
        # print(n, i, i.isalpha(), ord(i), ord(i)-64)
        if i.isalpha():
            value = ord(i) - 55
        else:
            value = int(i)
        to_add = value * 10 ** n
        atom_number += to_add
        # print(to_add, atom_number)
    return atom_number


def desc_from_fn(file_path, get_all=False):
    head, tail = os.path.split(file_path)
    desc, extension = tail.split(os.extsep)
    if get_all:
        return [head, tail, desc, extension]
    else:
        return desc


def bio_cif_to_pdb(ciffile):
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB import PDBIO
    import os
    import shutil
    # get the filename info
    head, tail, desc, extension = desc_from_fn(ciffile, get_all=True)
    print([head, tail, desc, extension])
    # copy the cif file
    shutil.copy(ciffile, 'temp.cif')
    # BioPython objects
    parser = MMCIFParser()
    # read the structure
    structure = parser.get_structure(desc, 'temp.cif')
    print(structure)
    print(structure.id)
    # make a subdirectory to put the separate chains
    pdb_dir = os.path.join(head, desc)
    if not os.path.exists(pdb_dir):
        os.mkdir(pdb_dir)
    pdbfile = os.path.join(pdb_dir, '{0}.pdb'.format(desc))
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdbfile)
    os.remove('temp.cif')


def my_cif_to_pdb(ciffile, target_dir, pdb_id):
    # a quick and dirty cif parser
    # a cif file is similar to a star file
    # we're only interested in the ATOM info
    # so we'll search for the loop_ members with _atom_site
    # to get the id of the various elements
    atom_site = []
    which_loop = '_atom_site.'
    alphanum = generate_alphanum()
    fo = ''
    last_chain = ''
    last_model = ''
    atom_n = 1
    chain_n = 0
    chains = []
    with open(ciffile) as f:
        for line in f:
            line = line[:-1]  # remove newline char
            #print(line)
            #print(line[:4])
            if which_loop in line:
                #print(line)
                line = line.replace(which_loop, '')
                line = line.replace(' ', '')  # remove any whitespace
                #print(line)
                atom_site.append(line)
            if line[:4] == 'ATOM':
                #print(atom_site)
                atom_list = uncif_atom(line, atom_site)
                if atom_list[3] != last_chain:
                    #print(atom_list[3], last_chain)
                    #print(atom_list[12], last_model)
                    last_chain = atom_list[3]
                    atom_n = 1
                    if atom_list[12] != last_model:
                        chain_n = 0
                    else:
                        chain_n += 1
                        if chain_n > 61:
                            chain_n = 0
                    last_model = atom_list[12]
                    out_fn = "{0}_{1}.pdb".format(pdb_id, last_chain)
                    if not isinstance(fo, str):
                        chains.append(out_fn)
                        fo.write('TER')
                        fo.close()
                    fo = open(os.path.join(target_dir, out_fn), 'w')
                atom_list[0] = str(atom_n)
                atom_list[3] = alphanum[chain_n]
                pdb_atom = atom_to_pdb(*atom_list)
                fo.write(pdb_atom + '\n')
                atom_n += 1
    chains.append(out_fn)
    fo.write('TER')
    fo.close()
    #print(atom_site)
    # return the number of chains that we processed
    return chains


def parse_atom_site_loop():
    # from PDB 4v88
    loop1 = """_atom_site.group_PDB
    _atom_site.id
    _atom_site.type_symbol
    _atom_site.label_atom_id
    _atom_site.label_alt_id
    _atom_site.label_comp_id
    _atom_site.label_asym_id
    _atom_site.label_entity_id
    _atom_site.label_seq_id
    _atom_site.pdbx_PDB_ins_code
    _atom_site.Cartn_x
    _atom_site.Cartn_y
    _atom_site.Cartn_z
    _atom_site.occupancy
    _atom_site.B_iso_or_equiv
    _atom_site.pdbx_formal_charge
    _atom_site.auth_seq_id
    _atom_site.auth_comp_id
    _atom_site.auth_asym_id
    _atom_site.auth_atom_id
    _atom_site.pdbx_PDB_model_num"""
    # from PyMOL
    loop2 = """_atom_site.group_PDB
    _atom_site.id
    _atom_site.type_symbol
    _atom_site.label_atom_id
    _atom_site.label_alt_id
    _atom_site.label_comp_id
    _atom_site.label_asym_id
    _atom_site.label_entity_id
    _atom_site.label_seq_id
    _atom_site.pdbx_PDB_ins_code
    _atom_site.Cartn_x
    _atom_site.Cartn_y
    _atom_site.Cartn_z
    _atom_site.occupancy
    _atom_site.B_iso_or_equiv
    _atom_site.pdbx_formal_charge
    _atom_site.auth_asym_id
    _atom_site.pdbx_PDB_model_num"""
    loop = []
    for l in loop1.split('\n'):
        print(l)
        which = '_atom_site.'
        if which in l:
            l = l.replace(which, '')
            loop.append(l)
    print(loop)
    return loop


def uncif_atom(line, loop):
    # test_line = "ATOM   1      P  P     . U   A   1  1    ? -88.729  32.079   65.623  1.00 116.34 ? 1    U   A2 P     1 "
    # test_line = "ATOM   102504 O  OP2   . A   JA  36 1251 ? -13.248  -13.171  -13.189 1.00 220.99 ? 1251 A   A1 OP2   1 "
    line = line.split(' ')
    line = [x for x in line if len(x)]
    # print(line)
    atom_serial_number = line[loop.index('id')]
    atom_name = line[loop.index('label_atom_id')].replace('\"', '')
    residue_name = line[loop.index('label_comp_id')]
    # this can be either 'label_asym_id' or 'auth_asym_id'
    if 'auth_asym_id' in loop:
        chain = line[loop.index('auth_asym_id')]
    else:
        chain = line[loop.index('label_asym_id')]
    residue_number = line[loop.index('label_seq_id')]
    charge = ""
    x = line[loop.index('Cartn_x')]
    y = line[loop.index('Cartn_y')]
    z = line[loop.index('Cartn_z')]
    occupancy = line[loop.index('occupancy')]
    bfactor = line[loop.index('B_iso_or_equiv')]
    atom_symbol = line[loop.index('type_symbol')]
    model = line[loop.index('pdbx_PDB_model_num')]
    atom_list = [atom_serial_number, atom_name, residue_name, chain, residue_number,
                 x, y, z, occupancy, bfactor, atom_symbol, charge, model]
    return atom_list


def atom_to_pdb(atom_serial_number, atom_name, residue_name, chain,
                residue_number, x, y, z, occupancy, bfactor, atom_symbol, charge, model):
    """
    there are 15 entries to a PDB
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    -------------------------------------------------------------------------------------
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 20        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.
    """
    pdb_atom = "ATOM  "
    pdb_atom += atom_serial_number[:5].rjust(5)
    pdb_atom += atom_name.rjust(5)
    pdb_atom += " "  # alternate location indicator
    pdb_atom += residue_name.rjust(3)
    pdb_atom += chain[0].rjust(2)
    pdb_atom += residue_number.rjust(4)
    pdb_atom += "    "  # code for insertion of residues
    pdb_atom += x.rjust(8)
    pdb_atom += y.rjust(8)
    pdb_atom += z.rjust(8)
    pdb_atom += occupancy.rjust(6)
    pdb_atom += bfactor.rjust(6)
    pdb_atom += "".rjust(10)
    pdb_atom += atom_symbol.rjust(2)
    pdb_atom += charge.rjust(2)
    return pdb_atom


def pdb_line_to_(line):
    atom_name = line[13:16].replace(' ', '')
    residue_name = line[17:20].replace(' ', '')
    chain = line[20:22].replace(' ', '')
    resid = int(line[22:26])
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    print(atom_name, residue_name, chain, resid, x, y, z)
    return atom_name, residue_name, chain, resid, x, y, z


def split_pdb(pdb_fn, target_dir, pdb_id):
    f = open(pdb_fn)
    data = f.read()
    data = data.split('\n')
    # a dict to store the separate chains
    pdb_chains = {}
    chains = []
    last_chain = ''
    this_chain = []
    for d in data:
        #print(d[0:4])  # columns 1-4 in the description is 1-1:4
        if 'ATOM' in d:
            chain = d[20:22].replace(' ', '')
            # if there is a new chain, save the last one
            if chain != last_chain:
                # only save the chain if the last one is populated
                if len(this_chain):
                    pdb_chains[last_chain] = this_chain
                    this_chain = []  # reset the chain
            else:
                this_chain.append(d)
            last_chain = chain
    pdb_chains[last_chain] = this_chain
    # save the separated chains
    for k in pdb_chains.keys():
        this_chain = pdb_chains[k]
        this_fn = "{0}_{1}.pdb".format(pdb_id, k)
        chains.append(this_fn)
        f = open(os.path.join(target_dir, this_fn), 'w')
        for line in this_chain:
            f.write(line + '\n')
        f.close()
    # return the number of chains processed
    return chains


def model_to_chains(filename, target_dir, pdb_id, extension):
    # with a pdb or mmCIF file with filename, parse the model
    # and split into separate chains and save in target_dir
    print(filename, target_dir, extension)
    if extension == 'pdb':
        # parse pdb and split by chain
        chains = split_pdb(filename, target_dir, pdb_id)
    elif extension == 'cif':
        # parse mmCIF and split by chain, save as pdb
        chains = my_cif_to_pdb(filename, target_dir, pdb_id)
    else:
        raise(Exception("unknown file type"))
    return chains


def summarise_coords_protein(coords):
    # return only the CA for now
    these_coords = [coords['CA']]
    return these_coords


def summarise_coords_rna(coords):
    #for k in coords.keys():
    #    print(k)
    ribose_ring = []
    rr_elements = ['C1\'', 'C2\'', 'C3\'', 'C4\'', 'O4\'']
    for r in rr_elements:
        ribose_ring.append(coords[r])
    ribose_ring = np.array(ribose_ring)
    ribose_ring = np.mean(ribose_ring, axis=0)
    ribose_ring = ribose_ring.tolist()
    if 'N9' in coords.keys():
        #print('purine')
        base_elements = ['C4', 'C5', 'C8', 'N7', 'N9']
    else:
        #print('pyrimidine')
        base_elements = ['C2', 'C4', 'C5', 'C6', 'N3']
    base_ring = []
    for r in base_elements:
        try:
            base_ring.append(coords[r])
        except KeyError:
            count = 0
            for k in coords.keys():
                print(k)
                if k in base_elements:
                    count += 1
                else:
                    print('base element missing')
            # if the base is completely missing, just use the C1' of the ribose as a proxy
            # this shouldn't happen for many bases in the model
            if count == 0:
                base_ring.append(coords['C1\''])
            #raise Exception('base ring not complete')

    base_ring = np.array(base_ring)
    base_ring = np.mean(base_ring, axis=0)
    base_ring = base_ring.tolist()
    if 'C5\'' in coords.keys():
        these_coords = [ribose_ring, base_ring, coords['O5\''], coords['O3\''], coords['C5\'']]
    else:
        # a missing C5' is usually at the 5' end of the molecule
        # if this is the case, O5' might be missing too
        # c5_proxy = [['O5\''], ['C3\'']]  # proxy 1
        c5_proxy = ['C4\'']  # proxy 2
        try:
            c5 = [coords[x] for x in c5_proxy]
        except KeyError:
            count = 0
            for k in coords.keys():
                print(k)
                if k in c5_proxy:
                    count += 1
                else:
                    print('C5 proxy element missing')
            raise Exception('C5 incomplete')
        c5 = np.array(c5)
        c5 = np.mean(c5, axis=0)
        c5 = c5.tolist()
        #these_coords = [ribose_ring, base_ring, coords['O5\''], coords['O3\''], c5]
        these_coords = [ribose_ring, base_ring, c5, coords['O3\''], c5]
    #print(these_coords)
    return these_coords


def get_coords(input_pdb, ignore_chains=False):
    head, tail, pdb_id, extension = input_pdb
    print(extension)
    f = open(input_pdb)
    data = f.read()
    data = data.split('\n')
    # get the PDB id from the name
    # a dicts to store info
    this_pdb = []  # each residue info will be added to this
    this_res = {}  # each residue will be added to this
    # variables to store temporary info
    last_chain = ''  # to check if the chain changes, this isn't handled
    last_resid = 0  # resid will never be zero, always starts >= 1
    # variables for parsing mmCIFs
    cif_header = ''
    in_loop = False
    chain = 'z'
    for d in data:
        # print(d[0:4])  # columns 1-4 in the description is 1-1:4
        # cycle through all the 'ATOM' lines in the file
        if 'ATOM' in d[:5]:
            if '.pdb' in input_pdb:
                print(d)
                atom_name = d[13:16].replace(' ', '')
                residue_name = d[17:20].replace(' ', '')
                chain = d[20:22].replace(' ', '')
                resid = int(d[22:26])
                x = float(d[30:38])
                y = float(d[38:46])
                z = float(d[46:54])
                print(pdb_id, chain, resid, residue_name, atom_name, x, y, z)
            elif '.cif' in input_pdb:
                in_loop = False
                print(d)
                line = d.split(' ')
                line = [x for x in line if x != '']
                atom_name = line[3]
                residue_name = line[5]
                chain = line[16]
                resid = int(line[8])
                x = float(line[10])
                y = float(line[11])
                z = float(line[12])
                print(pdb_id, chain, resid, residue_name, atom_name, x, y, z)
            # add the atom info to the info for this residue
            this_res[atom_name] = [x, y, z]  # save the name and coords of atom
            # if there is a new chain, save the last one
            if chain != last_chain:
                # only save the chain if the last one is populated
                if len(last_chain):
                    raise (Exception("This script only handles PDBs with one chain"))
                else:
                    last_chain = chain
            # if we're on to a new residue, need to summarise and save info
            if last_resid != resid:
                if len(this_res.keys()) > 1:
                    # print(pdb_id, chain, resid, residue_name)
                    # print(this_res)
                    for r in this_res:
                        print(r, this_res[r])
                    coords = summarise_coords_rna(this_res)
                    res_info = [resid, residue_name, coords]
                    this_pdb.append(res_info)
                    # input('...')
            last_resid = resid
            last_chain = chain
        # for getting mmCIF header info
        else:
            if '.cif' in input_pdb:
                if 'loop' in d:
                    in_loop = True
                if in_loop:
                    cif_header += d + '\n'
    f.close()
    return this_pdb


def read_pdb(pdb_fn, gap_char='x'):
    # from a single PDB chain, we will retrieve
    # sequence, xyz coordinates
    f = open(pdb_fn)
    data = f.read()
    f.close()
    data = data.split('\n')
    this_res = {}
    this_pdb = []
    # a dict to store the separate chains
    chain = ''
    last_chain = ''
    this_seq = ''
    last_resid = 0
    is_protein = False
    these_coords = []
    # go through each line in the PDB file
    for d in data:
        # print(d[0:4])  # columns 1-4 in the description is 1-1:4
        if 'ATOM' in d:
            #print(d)
            # some PDBs have their atom numbers encoded with alphanumeric strings
            #print(d[6].isalpha())
            # we're ignoring atom number for now
            #if d[6].isalpha():
            #    atom_number = convert_pdb_string_to_number(d[6:11])
            #else:
            #    atom_number = int(d[6:11])
            atom_name = d[12:16].replace(' ', '')
            residue_name = d[17:20].replace(' ', '')
            chain = d[20:22].replace(' ', '')
            resid = int(d[22:26])
            x = float(d[30:38])
            y = float(d[38:46])
            z = float(d[46:54])
            #print(atom_number, atom_name, residue_name, chain, resid)
            #print(x, y, z)
            # for each atom, we store information on the atomic element and its coordinates
            this_res[atom_name] = [x, y, z]
            # get the sequence from the pdb as well
            # it is not likely that a residue will be less than one atom in length
            # so we don't need to capture the last
            if last_resid != resid:
                # if there's anything in this_res
                if len(this_res.keys()) > 1:
                    # print(pdb_id, chain, resid, residue_name)
                    # print(this_res)
                    #for r in this_res:
                    #    print(r, this_res[r])
                    # summarise coordinates or not?
                    if is_protein:
                        coords = summarise_coords_protein(this_res)
                    else:
                        coords = summarise_coords_rna(this_res)
                    these_coords.append(coords)
                    res_info = [resid, residue_name, coords]
                    #res_info = [resid, residue_name, this_res]
                    this_pdb.append(res_info)
                    this_res = {}
                # do we need to fill gaps in the sequence?
                # if so, the gap_char should be set to a unique char (not seen in AA or NA seqs, e.g. 'x')
                # if no gap_char is given, i.e., gap_char = '', then it will not be added to the sequence
                if resid - last_resid != 1:
                    #print(resid - last_resid)
                    gap = ''
                    for i in range(resid - last_resid - 1):
                        gap += gap_char
                    #print(len(gap), gap)
                    this_seq += gap
                # need to work out if it's a protein or nucleic acid
                if len(residue_name) == 1:
                    this_seq = this_seq + residue_name
                    is_protein = False
                elif len(residue_name) == 3:
                    this_seq = this_seq + aa_three_to_one[residue_name]
                    is_protein = True
                    # if we're on to a new residue, need to summarise and save info
                # for debugging
                #print(this_pdb, this_seq, is_protein)
                #for res in this_pdb:
                #    print(res)
                #input()

            #print(chain, last_chain)
            if chain != last_chain and last_chain != '':
                raise(Exception('We are only processing one chain at a time, why am I seeing another chain?'))

            last_resid = resid
            last_chain = chain
    # we could return this_pdb as well
    # but it contains redundant information
    return these_coords, this_seq, is_protein


def calc_vossvolvox(input_pdb, vvv_path='./bin/vossvolvox', clean_generated_files=True):
    import subprocess
    this_fn = input_pdb[:-4]
    #print(this_fn)
    stderr = ''
    # paths for vossvolvox
    xyzr_path = os.path.join(vvv_path, 'pdb_to_xyzr')
    vol_path = os.path.join(vvv_path, 'vol')
    vdw_path = os.path.join(vvv_path, 'vdw')

    # create output filenames and description
    # if it's a cif, we'll be changing the input_pdb name
    output_xyzr = this_fn + '.xyzr'
    output_vol = this_fn + '.vol'
    output_vdw = this_fn + '.vdw'
    #desc = desc_from_fn(input_pdb)
    #print(output_xyzr, output_vol, output_vdw)

    # check that paths exist
    for path in [xyzr_path, vol_path, vdw_path]:
        if not os.path.exists(path):
            print('{0} does not exist, cannot do VOSSVOLVOX'.format(path))
            return 0.0, 0.0

    # convert PDB to xyzr
    # e.g. ./pdb_to_xyzr 5xyi_a.pdb > 5xyi_a.xyzr
    cmd = [xyzr_path, input_pdb]
    #print(cmd)
    # note atmtypenumbers needs to be in working directory
    f = open(output_xyzr, 'w')
    result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    stderr += result.stderr
    f.close()

    # calculate volume
    # after doing make vol in src directory
    # e.g. Volume.exe -i 1a01-noions.xyzr -p 1.5 -g 0.5
    cmd = [vol_path, '-i', output_xyzr, '-p', '1.5', '-g', '0.5']
    f = open(output_vol, 'w')
    result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    stderr += result.stderr
    f.close()

    # calculate vdw
    cmd = [vdw_path, '-i', output_xyzr, '-p', '1.5', '-g', '0.5']
    f = open(output_vdw, 'w')
    result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    stderr += result.stderr
    f.close()

    # read the output files and return the vol and vdw metrics
    f = open(output_vol)
    data = f.read()
    f.close()
    data = data.split('\t')
    #print(data)
    vol = float(data[2])

    f = open(output_vdw)
    data = f.read()
    f.close()
    data = data.split('\t')
    #print(data)
    vdw = float(data[2])

    # delete files from analysis?
    if clean_generated_files:
        os.remove(output_xyzr)
        os.remove(output_vol)
        os.remove(output_vdw)

    return vol, vdw, stderr

