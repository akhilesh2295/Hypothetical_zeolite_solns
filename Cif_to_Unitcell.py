#!/usr/bin/env python3
#
# Python script that converts a CIF file (Crystallographic Information File)
# into a configuration file for Gromacs or LAMMPS.


import sys

import logging
# logging.basicConfig(level=logging.DEBUG)

from math import *

from CifFile import CifFile, CifBlock


def float_with_error(x):
    """
    some value in cif accompanies error like "1.234(5)
    """
    if "?" in x:
        return 0
    pos = x.find("(")
    if pos >= 0:
        x = x[:pos]
    return float(x)


# =============================================================================

# Tell the user how to use this script, and exits the script.
def print_usage():
    # logging.getLogger().setLevel(logging.INFO)
    # logging.root.setLevel(logging.NOTSET)
    logger = logging.getLogger()
    logger.info('*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *')
    logger.info('                             --- HOW TO USE  ---')
    logger.info('Examples:')
    logger.info('  %s  -i crystal.cif  -o unitcell.xyz  -r' % sys.argv[0])
    logger.info('  %s  -i crystal.cif  -o conf.gro  -b 5 5 5  -r' % sys.argv[0])
    logger.info('  %s  -i crystal.cif  -o unitcell.cif' % sys.argv[0])
    logger.info('')
    logger.info('List of optional arguments:')
    logger.info('  -i  filename     CIF file with description of crystal')
    logger.info('  -o  filename     Output configuration file.  Extension must be one of the following:')
    logger.info('                     .xyz         XYZ file')
    logger.info('                     .lammpstrj   LAMMPS trajectory file')
    logger.info('                     .gro         Gromacs file')
    logger.info('                     .cif         copy of original CIF file with additional atoms')
    logger.info('                                  (looks nice with Jmol)')
    logger.info('  -b int int int   Box size in terms of the unit cell. For example, if the unit cell')
    logger.info('                   consists of 2 atoms and has size 1x2x3 A, then -b 4 5 6 creates a box')
    logger.info('                   of size 4x10x18 A consisting of 2x4x5x6=240 atoms. Default is -b 1 1 1')
    logger.info('  -r               Make the box rectangular (default is shape of the unit cell)')
    logger.info('*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *')
    # print('fail')
    exit(1)


# =============================================================================

# Shows an error message and the usage.
def print_error(msg):
    logger = logging.getLogger()
    logger.error('{0}'.format(msg))
    print_usage()


# =============================================================================

# Converts an "_atom_type_label" into an element name.
def extract_element(label):
    logger = logging.getLogger()

    elem2 = ['He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti',
             'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
             'Rb','Sr','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
             'Sb','Te','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
             'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt',
             'Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',
             'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']

    if (label[0:2] in elem2):
        return label[0:2]

    elem1 = ['H','B','C','N','O','F','P','S','K','V','Y','I','W','U']

    if (label[0] in elem1):
        return label[0]

    # logger.warn('could not convert "%s" into element name!' % label)
    return label

# =============================================================================

# Basic XYZ format.
def write_xyz(atoms, box, f):

    # Write the number of atoms.
    N = len(atoms)
    f.write('%d\n' % N)

    # Write a comment line with the box size.
    if (len(box) == 6):
        # box = (ax,bx,by,cx,cy,cz)
        ax = box[0]
        bx = box[1]
        by = box[2]
        cx = box[3]
        cy = box[4]
        cz = box[5]
        f.write('Crystal created from CIF file. Box vectors:')
        f.write(' a= %10.5f %10.5f %10.5f' % (ax, 0.0, 0.0))
        f.write(' b= %10.5f %10.5f %10.5f' % (bx, by, 0.0))
        f.write(' c= %10.5f %10.5f %10.5f\n' % (cx, cy, cz))
    else:
        # box = (ax,by,cz)
        f.write('Crystal created from CIF file. Box size:') 
        f.write(' %10.5f %10.5f %10.5f\n' % box)

    # Write atom data (units are Angstroms).
    # The argument "atoms" has format ('Si', x, y, z) for example
    for i in range(N):
        f.write('%-10s %10.6f %10.6f %10.6f\n' % atoms[i])


# =============================================================================


# Read CIF file, and extract the necessary info in the form of a dictionary.
# E.g., the value of "_cell_volume" can be found with data['_cell_volume'].
def read_cif_for_obabel(fNameIn):
    logger = logging.getLogger()

    data = {}

    # print('inside read_cif')
    # print(fNameIn)
    # Open the CIF file and read all the lines into a list of strings.
    try:
        # print(fNameIn)
        # print(os.getcwd)
        f = open(fNameIn, 'r')
        # print('opened file successfully')
        lines = []
        oko_fix = False
        for line in f:
            stripped = line.strip()
            if (len(stripped) > 0):
                #fix for OKO
                #http://stackoverflow.com/questions/79968/split-a-string-by-spaces-preserving-quoted-substrings-in-python
                import shlex
                if stripped[0] == '_':
                    cols = shlex.split(stripped) #shlex protects the quotation
                    if len(cols) > 2:
                        stripped = "{0} '{1}'".format(cols[0]," ".join(cols[1:]))
                        logger.debug("Line (Fixed): {0}".format(stripped))
                        oko_fix = True
                lines.append(stripped)
        if oko_fix:
            fixedfilename = fNameIn + ".fix"
            fixedfile = open(fixedfilename,"w")
            fixedfile.write("\n".join(lines)+"\n")
            fixedfile.close()
            fNameIn = fixedfilename
        f.close()
    except:
        logger.error("Failed to open CIF file '{0}'".format(fNameIn))
        sys.exit()

    # Use the CifFile parser to extract the data.  Although there might be
    # multiple data blocks, we'll only use the first one.
    cif_file = CifFile(fNameIn)
    
    for db in cif_file:
        data_block = db
        break


    try:

        # Extract some parameters, and convert them to floats.
        data['_cell_length_a']    = float_with_error(data_block['_cell_length_a'])
        data['_cell_length_b']    = float_with_error(data_block['_cell_length_b'])
        data['_cell_length_c']    = float_with_error(data_block['_cell_length_c'])
        data['_cell_angle_alpha'] = float_with_error(data_block['_cell_angle_alpha'])
        data['_cell_angle_beta']  = float_with_error(data_block['_cell_angle_beta'])
        data['_cell_angle_gamma'] = float_with_error(data_block['_cell_angle_gamma'])
        if data_block.has_key('_cell_volume'):
            data['_cell_volume']      = float_with_error(data_block['_cell_volume'])

        # Get the symbolic operations that define the space group.  In a CIF file
        # that's the part that looks like:
        #
        # loop_
        # _symmetry_equiv_pos_as_xyz
        #   'x,y,z'
        #   'y,x,2/3-z'
        #   '-y,x-y,2/3+z'
        #   '-x,-x+y,1/3-z'
        #   '-x+y,-x,1/3+z'
        #   'x-y,-y,-z'
        #
        # In some cases it's called "_space_group_symop_operation_xyz" apparently?!?!
        data['_symmetry_equiv_pos_as_xyz'] = []

        try:
            xyz = data_block["_symmetry_equiv_pos_as_xyz"]

        except KeyError:
            try:
                xyz = data_block["_space_group_symop_operation_xyz"]
            except KeyError:
                # logger.error("Missing item in CIF file: need either '_symmetry_equiv_pos_as_xyz' or '_space_group_symop_operation_xyz'. Sys exit commented out so that obabel output can be processed further. If this error shows without obabel use, investigate")
                logger.error("Critical error but this function is designed only for obabel. Investigate if not using obabel approach")
                # sys.exit()

        # for op_xyz in xyz:

        #     if (op_xyz[0] == '\''):
        #         data['_symmetry_equiv_pos_as_xyz'].append(op_xyz[1:-1])
        #     else:
        #         data['_symmetry_equiv_pos_as_xyz'].append(op_xyz)


        # Add x,y,z of the atoms to "data", but make sure to convert
        # e.g. "0.1549(8)" to "0.1549".
        data['_atom_site_label'] = data_block['_atom_site_label']

        data['_atom_site_fract_x'] = []
        for str_x in data_block['_atom_site_fract_x']:
            data['_atom_site_fract_x'].append( float_with_error(str_x))

        data['_atom_site_fract_y'] = []
        for str_y in data_block['_atom_site_fract_y']:
            data['_atom_site_fract_y'].append( float_with_error(str_y))

        data['_atom_site_fract_z'] = []
        for str_z in data_block['_atom_site_fract_z']:
            data['_atom_site_fract_z'].append( float_with_error(str_z))

    
    except KeyError as e:
        logger.error("Error!  Missing item in file.")
        logger.error(e)
        sys.exit()


    #print "ALL DATA:"
    #print data
    #print

    # Return the extracted data.
    return data



# Read CIF file, and extract the necessary info in the form of a dictionary.
# E.g., the value of "_cell_volume" can be found with data['_cell_volume'].
def read_cif(fNameIn):
    logger = logging.getLogger()

    data = {}

    # print('inside read_cif')
    # print(fNameIn)
    # Open the CIF file and read all the lines into a list of strings.
    try:
        # print(fNameIn)
        # print(os.getcwd)
        f = open(fNameIn, 'r')
        # print('opened file successfully')
        lines = []
        oko_fix = False
        for line in f:
            stripped = line.strip()
            if (len(stripped) > 0):
                #fix for OKO
                #http://stackoverflow.com/questions/79968/split-a-string-by-spaces-preserving-quoted-substrings-in-python
                import shlex
                if stripped[0] == '_':
                    cols = shlex.split(stripped) #shlex protects the quotation
                    if len(cols) > 2:
                        stripped = "{0} '{1}'".format(cols[0]," ".join(cols[1:]))
                        logger.debug("Line (Fixed): {0}".format(stripped))
                        oko_fix = True
                lines.append(stripped)
        if oko_fix:
            fixedfilename = fNameIn + ".fix"
            fixedfile = open(fixedfilename,"w")
            fixedfile.write("\n".join(lines)+"\n")
            fixedfile.close()
            fNameIn = fixedfilename
        f.close()
    except:
        logger.error("Failed to open CIF file '{0}'".format(fNameIn))
        sys.exit()

    # Use the CifFile parser to extract the data.  Although there might be
    # multiple data blocks, we'll only use the first one.
    cif_file = CifFile(fNameIn)
    
    for db in cif_file:
        data_block = db
        break


    try:

        # Extract some parameters, and convert them to floats.
        data['_cell_length_a']    = float_with_error(data_block['_cell_length_a'])
        data['_cell_length_b']    = float_with_error(data_block['_cell_length_b'])
        data['_cell_length_c']    = float_with_error(data_block['_cell_length_c'])
        data['_cell_angle_alpha'] = float_with_error(data_block['_cell_angle_alpha'])
        data['_cell_angle_beta']  = float_with_error(data_block['_cell_angle_beta'])
        data['_cell_angle_gamma'] = float_with_error(data_block['_cell_angle_gamma'])
        if data_block.has_key('_cell_volume'):
            data['_cell_volume']      = float_with_error(data_block['_cell_volume'])

        # Get the symbolic operations that define the space group.  In a CIF file
        # that's the part that looks like:
        #
        # loop_
        # _symmetry_equiv_pos_as_xyz
        #   'x,y,z'
        #   'y,x,2/3-z'
        #   '-y,x-y,2/3+z'
        #   '-x,-x+y,1/3-z'
        #   '-x+y,-x,1/3+z'
        #   'x-y,-y,-z'
        #
        # In some cases it's called "_space_group_symop_operation_xyz" apparently?!?!
        data['_symmetry_equiv_pos_as_xyz'] = []

        try:
            xyz = data_block["_symmetry_equiv_pos_as_xyz"]

        except KeyError:
            try:
                xyz = data_block["_space_group_symop_operation_xyz"]
            except KeyError:
                logger.error("Missing item in CIF file: need either '_symmetry_equiv_pos_as_xyz' or '_space_group_symop_operation_xyz'.")
                sys.exit()

        # Copy the x,y,z symmetry group operations.  Remove the quotes if there
        # are any.
        for op_xyz in xyz:

            if (op_xyz[0] == '\''):
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz[1:-1])
            else:
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz)


        # Add x,y,z of the atoms to "data", but make sure to convert
        # e.g. "0.1549(8)" to "0.1549".
        data['_atom_site_label'] = data_block['_atom_site_label']

        data['_atom_site_fract_x'] = []
        for str_x in data_block['_atom_site_fract_x']:
            data['_atom_site_fract_x'].append( float_with_error(str_x))

        data['_atom_site_fract_y'] = []
        for str_y in data_block['_atom_site_fract_y']:
            data['_atom_site_fract_y'].append( float_with_error(str_y))

        data['_atom_site_fract_z'] = []
        for str_z in data_block['_atom_site_fract_z']:
            data['_atom_site_fract_z'].append( float_with_error(str_z))

    
    except KeyError as e:
        logger.error("Error!  Missing item in file.")
        logger.error(e)
        sys.exit()


    #print "ALL DATA:"
    #print data
    #print

    # Return the extracted data.
    return data


# =============================================================================
# =============================================================================


def parse_commandline(filetypes=['.xyz', '.cif']):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Read arguments given.

    # Default settings.
    fNameIn = ''
    fNameOut = ''
    Nx = 1
    Ny = 1
    Nz = 1
    make_rect_box = False


    # Read the arguments.  We expect at least 4.
    if (len(sys.argv) <= 4):
        print_usage()

    i = 1
    while (i < len(sys.argv)):


        # Check if the name of the input file was given.
        if (sys.argv[i] == '-i'):

            # Make sure a file name is given.
            if (i+1 == len(sys.argv)):
                print_error('no input file name given')

            fNameIn = sys.argv[i+1]
            i = i + 2


        # Check if the name of the output file was given.
        elif (sys.argv[i] == '-o'):

            # Make sure a file name is given.
            if (i+1 == len(sys.argv)):
                print_error('no output file name given')

            # Check we have a valid file extension.
            fNameOut = sys.argv[i+1]
            unknown = True

            for ext in filetypes:
                if (fNameOut.endswith(ext)):
                    unknown = False

            if (unknown):
                print_error('unknown file extension of output file')

            i = i + 2


        # Check if the box size was given.
        elif (sys.argv[i] == '-b'):

            # Make sure 3 integers are given.
            if (i+3 >= len(sys.argv)):
                print_error('need 3 integers to indicate box size')

            Nx = int(sys.argv[i+1])
            Ny = int(sys.argv[i+2])
            Nz = int(sys.argv[i+3])

            if (Nx == 0  or  Ny == 0  or  Nz == 0):
                print_error('box size integers need to be larger than zero')

            i = i + 4


        # Check if the final configuration should be in a rectangular shape, or in
        # the same shape as the unit cell.
        elif (sys.argv[i] == '-r'):

            make_rect_box = True
            i = i + 1


        # Anything else is wrong.
        else:
            print_error('invalid argument "%s"' % sys.argv[i])

    return fNameIn, fNameOut, (Nx,Ny,Nz), make_rect_box

def read_and_process(fNameIn, Nbox, make_rect_box=False):
    logger = logging.getLogger()
    Nx,Ny,Nz = Nbox
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Read input file.


    # Make sure an input file was given.
    if (fNameIn == ''):
        print_error('no input file given.  Use:  -i filename')


    # Open the CIF file and read the data.
    # print('inside read_and_process')
    data = read_cif(fNameIn)


    # Extract lengths and angles from the CIF file.
    La = float(data['_cell_length_a'])
    Lb = float(data['_cell_length_b'])
    Lc = float(data['_cell_length_c'])
    alpha = radians(float(data['_cell_angle_alpha']))
    beta = radians(float(data['_cell_angle_beta']))
    gamma = radians(float(data['_cell_angle_gamma']))
    volume = 0
    if '_cell_volume' in data:
        volume = float(data['_cell_volume'])


    # Extract the symmetry operations.  This will be a list of strings such as:
    #    ['x,y,z', 'y,x,2/3-z', '-y,x-y,2/3+z', '-x,-x+y,1/3-z', ... ]
    ops = data['_symmetry_equiv_pos_as_xyz']

    # For proper evaluation, we need to convert "2/3" into "2./3", etc. to prevent
    # integer division which would turn e.g. 2/3 into 0.
    for i in range(len(ops)):
        ops[i] = ops[i].replace("0/", "0./") # also for e.g. 10/9
        ops[i] = ops[i].replace("1/", "1./")
        ops[i] = ops[i].replace("2/", "2./")
        ops[i] = ops[i].replace("3/", "3./")
        ops[i] = ops[i].replace("4/", "4./")
        ops[i] = ops[i].replace("5/", "5./")
        ops[i] = ops[i].replace("6/", "6./")
        ops[i] = ops[i].replace("7/", "7./")
        ops[i] = ops[i].replace("8/", "8./")
        ops[i] = ops[i].replace("9/", "9./")
    #    ops[i] = ops[i].replace("/", "./")


    # Get the atom labels and coordinates.
    labels = data['_atom_site_label']
    fX = [ float(s) for s in data['_atom_site_fract_x'] ]
    fY = [ float(s) for s in data['_atom_site_fract_y'] ]
    fZ = [ float(s) for s in data['_atom_site_fract_z'] ]

    # Create a list of 4-tuples, where each tuple is an atom:
    #   [ ('Si', 0.4697, 0.0, 0.0),  ('O', 0.4135, 0.2669, 0.1191),  ... ]
    atoms = [ (labels[i], fX[i], fY[i], fZ[i]) for i in range(len(labels)) ]

    # Make sure that all atoms lie within the unit cell.  Also convert names such
    # as 'Oa1' into 'O'.
    for i in range(len(atoms)):
        (name,xn,yn,zn) = atoms[i]
        # print(xn,yn,zn)
        xn = (xn + 10.0) % 1.0
        yn = (yn + 10.0) % 1.0
        zn = (zn + 10.0) % 1.0
        name = extract_element(name)
        atoms[i] = (name,xn,yn,zn)

    # Update the user.
    # logger.info('Loaded a CIF file with {0} atom coordinates and {1} symmetry operations.'.format(len(atoms), len(ops)))

    # Two atoms are on top of each other if they are less than "eps" away.
    eps = 0.01  # in Angstrom

    # print(atoms)
    # For each atom, apply each symmetry operation to create a new atom.
    imax = len(atoms)
    i=0
    while (i < imax):

        label,x,y,z = atoms[i]

        for op in ops:

            # Python is awesome: calling e.g. eval('x,y,1./2+z') will convert the
            # string into a 3-tuple using the current values for x,y,z!
            # logger.debug("OP: {0}".format(op))
            # print("OP: {0}".format(op))
            op = op.lower()
            xn,yn,zn = eval(op)
            # print(xn,yn,zn)
            # Make sure that the new atom lies within the unit cell.
            xn = (xn + 10.0) % 1.0
            yn = (yn + 10.0) % 1.0
            zn = (zn + 10.0) % 1.0
            # print(xn,yn,zn)
            # Check if the new position is actually new, or the same as a previous
            # atom.
            new_atom = True
            for at in atoms:
                if (abs(at[1]-xn) < eps  and  abs(at[2]-yn) < eps  and  abs(at[3]-zn) < eps):
                    # if (abs(11.589-xn) < eps  and  abs(2.733-yn) < eps  and  abs(11.289-zn) < eps):
                    #     print(xn)
                    new_atom = False

                    # Check that this is the same atom type.
                    if (at[0] != label):
                        print_error('invalid CIF file: atom of type %s overlaps with atom of type %s' % (at[0],label))

            # If the atom is new, add it to the list!
            if (new_atom):
                # print(label, xn, yn, zn)
                atoms.append( (label,xn,yn,zn) )  # add a 4-tuple
        # print(atoms[6])

        # Update the loop iterator.
        i = i + 1
        imax = len(atoms)


    # Sort the atoms according to type alphabetically.
    atoms = sorted(atoms, key=lambda at: at[0])
    atoms.reverse()


    # Done with creating the unit cell.  Update the user.
    # logger.info('Created a unit cell consisting of %d atoms.' % len(atoms))

    # logger.info('Fractional coordinates:')
    # for atom in atoms:
    #     logger.info('%10s  %.3f  %.3f  %.3f' % atom)



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Create a larger box made of several unit cells: the super cell.


    atomlist = []
    for atom in atoms:

        # Get label and fractional coordinates.
        label,xf,yf,zf = atom

        for i in range(Nx):
                x = i+xf

                for j in range(Ny):
                    y = j+yf

                    for k in range(Nz):
                        z = k+zf
                        atomlist.append( (label,x,y,z) ) # add 4-tuple

    atoms = atomlist
    # print(atoms)

    # If the user wants us to create a copy of the current CIF file, with
    # additional atoms, then do that.  Note that the atoms here have *fractional*
    # coordinates!
    ##if (fNameOut.endswith('.cif')):

    ##    write_cif(fNameIn, atoms, fNameOut)

    ##    logger.info('Done writing extended CIF file (%d atoms in total).' % len(atoms))
    ##    exit(0)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Convert the fractional coordinates into real coordinates.


    # The primitive vectors a,b,c are such that 
    #
    #   cos(alpha) = b.c / |b||c|
    #   cos(beta)  = a.c / |a||c|
    #   cos(gamma) = a.b / |a||b|
    #
    # with the convention
    #
    #   a = La*xhat
    #   b = bx*xhat + by*yhat
    #   c = cx*xhat + cy*yhat + cz*zhat
    #
    cosa = cos(alpha)
    sina = sin(alpha)
    cosb = cos(beta)
    sinb = sin(beta)
    cosg = cos(gamma)
    sing = sin(gamma)

    cosa2 = cosa * cosa
    cosb2 = cosb * cosb
    sing2 = sing * sing

    ax = La

    bx = Lb * cosg
    by = Lb * sing

    cx = Lc * cosb
    cy = Lc * (cosa - cosg*cosb) / sing
    cz = Lc * sqrt( 1 - (cosa2 + cosb2 - 2*cosg*cosb*cosa) / sing2 )


    # Use the volume to check if we did the vectors right.
    V = ax*by*cz
    if volume != 0 and ( abs(V - volume) > 0.1):
        print_error('volume does not match that calculated from primitive vectors')


    # Check if we have a rectangular box.
    if (abs(bx) < eps  and  abs(cx) < eps  and abs(cy) < eps):
        make_rect_box = True


    # Update the user.
    # logger.info('The primitive unit cell vectors are:')
    # logger.info('   a = [%6.4f, %6.4f, %6.4f]' % (ax,0,0))
    # logger.info('   b = [%6.4f, %6.4f, %6.4f]' % (bx,by,0))
    # logger.info('   c = [%6.4f, %6.4f, %6.4f]' % (cx,cy,cz))
    # logger.info('This gives a volume of %f A^3 (CIF file indicates it is %f A^3)' % (V,volume))

    # print(Nx, Ny, Nz)
    # Determine the box size.
    Lx = Nx * La
    Ly = Ny * Lb
    Lz = Nz * Lc
    # print(atoms)
    # print(len(atoms))
    for i in range(len(atoms)):

        # Get label and fractional coordinates.
        label,xf,yf,zf = atoms[i]

        xa = xf * ax  # contribution of a-vector to the x-coordinate of this atom
        #ya = 0       # a-vector has no y-component, so does not affect y of atom
        #za = 0       # a-vector has no z-component, so does not affect z of atom

        xb = yf * bx  # contribution of b-vector to the x-coordinate of this atom
        yb = yf * by  # contribution of b-vector to the y-coordinate of this atom
        #zb = 0       # b-vector has no z-component, so does not affect z of atom

        xc = zf * cx  # contribution of c-vector to the x-coordinate of this atom
        yc = zf * cy  # contribution of c-vector to the y-coordinate of this atom
        zc = zf * cz  # contribution of c-vector to the z-coordinate of this atom

        # Add all contributions.
        xn = xa + xb + xc
        yn = yb + yc
        zn = zc
        # print('xn, yn, zn',xn,yn,zn)
        # print('xn/yn', xn/yn)
        # following if condition makes the atoms go from negative to the other side of the lattice requested box
        if (make_rect_box):
            xn = (xn + Lx) % Lx
            yn = (yn + Ly) % Ly
            zn = (zn + Lz) % Lz

        atoms[i] = (label, xn, yn, zn)

    # print(atoms)
    # Determine the box-vector.
    if make_rect_box:
        box = (Lx, Ly, Lz)
    else:
        #box = (Lx, Ly, Lz, Nx*cx, Ny*cy, Nz*cz)
        box = (ax*Nx, bx*Ny, by*Ny, cx*Nz, cy*Nz, cz*Nz)
    return atoms, box, La, Lb, Lc, cosa, sina, cosb, sinb, cosg, sing
    # return atoms, box


def read_and_report(fNameIn):
    # Make sure an input file was given.
    if (fNameIn == ''):
        print_error('no input file given.  Use:  -i filename')


    # Open the CIF file and read the data.
    # print('inside read_and_process')
    data = read_cif_for_obabel(fNameIn)


    # Extract lengths and angles from the CIF file.
    La = float(data['_cell_length_a'])
    Lb = float(data['_cell_length_b'])
    Lc = float(data['_cell_length_c'])
    alpha = radians(float(data['_cell_angle_alpha']))
    beta = radians(float(data['_cell_angle_beta']))
    gamma = radians(float(data['_cell_angle_gamma']))
    volume = 0
    if '_cell_volume' in data:
        volume = float(data['_cell_volume'])

    cosa = cos(alpha)
    sina = sin(alpha)
    cosb = cos(beta)
    sinb = sin(beta)
    cosg = cos(gamma)
    sing = sin(gamma)

    cosa2 = cosa * cosa
    cosb2 = cosb * cosb
    sing2 = sing * sing

    ax = La

    bx = Lb * cosg
    by = Lb * sing

    cx = Lc * cosb
    cy = Lc * (cosa - cosg*cosb) / sing
    cz = Lc * sqrt( 1 - (cosa2 + cosb2 - 2*cosg*cosb*cosa) / sing2 )


    # Use the volume to check if we did the vectors right.
    V = ax*by*cz
    if volume != 0 and ( abs(V - volume) > 0.1):
        print_error('volume does not match that calculated from primitive vectors')


    return La, Lb, Lc, ax, bx, by, cx, cy, cz


def write(fNameOut, atoms, box, fNameIn=""):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Create the output file.

    try:
        fOut = open(fNameOut, 'w')

        if (fNameOut.endswith('.xyz')):
            write_xyz(atoms, box, fOut)

        elif (fNameOut.endswith('.cif')):
            write_cif(fNameIn, atoms, fNameOut)

    except:
        print_error('Failed to write to output file')


    # logger.info('Created output file %s (%d atoms in total).' % (fNameOut, len(atoms)))
    fOut.close()


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def main():
    fNameIn, fNameOut, Nbox, make_rect_box = parse_commandline()
    atoms, box = read_and_process(fNameIn, Nbox, make_rect_box)
    write(fNameOut, atoms, box, fNameIn=fNameIn)

if __name__ == "__main__":
    main()