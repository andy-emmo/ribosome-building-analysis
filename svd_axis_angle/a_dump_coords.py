from pymol import cmd
import json  

def dump_coords():
    # here we simply want to dump the points of the current selection
    # just the phosphates are enough to do this satisfactorily
    # for our structure we select the ribosome small subunit body with
    # select c0_St and (resi 1-894 or resi 1316-1454) and name P
    # and then 
    # dump_coords()
    
    coords = cmd.get_coords('sele', 1)
    f = open('a_coords.json', 'w')
    f.write(json.dumps(coords.tolist()))
    f.close()
    
    # next, use the python script with prefix b_to do SVD and return the translation and 
    # Euler angles to move the structure to the orgin and align the structure so that the 
    # least variation is aligned with the (0,0,1) axis 
    
