# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 10:09:16 2020

@author: Adi
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from  matplotlib.collections import LineCollection
import collections
import scipy.io
import os
import shutil


#-------------Definition--------------------

def write_new_projwfc(filproj, nspin) :
    '''
    Writing To a new file without unessacary things

    Parameters
    ----------
    nspin : Integer
        DESCRIPTION :   output of read_scf_out()
                        1 is  spin unpolarized
                        2 is  spin polarized 
                        3 is  non-collinear and spin-orbit 
                        
    filproj :  String
        filproj as in projwfc.x input
        This File is made after running projwfc.x 
        File format is 
        
         State#  Atom# Atom Orbital   wfc#  l  j  m
         Where, 
         l = [0, 1, 2, ...] --->  [s, p, d, ....]
         j = Total angular momentum
         m = One Electronic State in each Orbital with maximum of 2
        

    Returns
    -------
    Writes a new file in name of projwfc_up.new
    Retruns the name of the new file

    '''       
    if nspin == 3 :
        split      =   '    T    T\n'
        footer     =   ['.projwfc_up'] 
        print('NON-COLLINEAR WITH SPIN-ORBIT CASE')
        
    elif nspin == 1 :
        split      =   '    F    F\n'
        footer     =   ['.projwfc_up']
        print('UNPOLARIZED SPIN CASE')
        
    elif nspin == 2 :
        split      =   '    F    F\n' 
        footer    =   ['.projwfc_up', '.projwfc_down']
        print('COLINEAR POLARIZED SPIN CASE')
        
    else  :  
        print('Problem with nspin in definition write_new_projwfc()')
        sys.exit(0)
        
    filename    =   [] 
    for foot in footer :
        file       =     filproj + foot
        projwfc    =     open(file).read().split(split)[1]
        if len(projwfc) > 5 :  
            file      +=     '.new' 
            filename.append(file)
            towrite     =    open(file, 'w', newline='\n')
            towrite.write(projwfc)
            towrite.close()
    
        else :  
            print('Split string not found in %s')
            sys.exit(0)
    
    return(filename)

def read_states_projwfc_out(filename, nspin) :
    '''
    To read the state informations from Output of Projwfc.x

    Parameters
    ----------
    nspin : Integer
        DESCRIPTION :   output of read_scf_out()
                        1 is  spin unpolarized
                        2 is  spin polarized 
                        3 is  non-collinear and spin-orbit 
    filename : String
        DESCRIPTION : Contains Terminal Output of Projwfc.x

    Returns
    -------
    State information as lists
    if nspin is 1 or  2 :
    states format -->  [state, atom#, atom, wfc#, l, m]
    if nspin is 3 :
    states format -->  [state, atom#, atom, wfc#, l, j, m_j]

    '''
    data     =   open(filename).read()\
                             .split('Atomic states used for projection')[1]\
                             .split(' k =   ')[0]\
                             .split('state #')[1:]    
    states   =   []
    if nspin == 3 :
        for state in range(len(data)) :
            arr  = []
            arr.append(int(data[state].split(':')[0]))
            arr.append(int(data[state].split('atom')[1].split('(')[0]))
            arr.append(str(data[state].split('(')[1].split(')')[0]))
            arr.append(int(data[state].split('wfc')[1].split('(')[0]))
            foo     =   data[state].split('(')[2].split(')')[0].split()
            foo     =   [i.split('=') for i in foo]
            index   =   0
            while True :
                if foo[index][0] == 'l'   : 
                    l  =  (int(foo[index][1]))
                    index  += 1
                if foo[index][0] == 'j'   : 
                    j   = (float(foo[index][1]))
                    index  += 1
                if foo[index][0] == 'm_j' : 
                    if foo[index][1]  ==  '':
                        m_j   =   (float(foo[index+1][0]))
                    else :
                        m_j   =   (float(foo[index][1]))
                    index  += 1
                if index  ==  3 :
                    break 
            for i in [l, j, m_j] :
                arr.append(i)
            states.append(arr) 
        
    elif nspin == 1 or nspin == 2 :
        for state in range(len(data)) :
            arr  = []
            arr.append(int(data[state].split(':')[0]))
            arr.append(int(data[state].split('atom')[1].split('(')[0]))
            arr.append(str(data[state].split('(')[1].split(')')[0]))
            arr.append(int(data[state].split('wfc')[1].split('(')[0]))
            foo     =   data[state].split('(')[2].split(')')[0].split()
            foo     =   [i.split('=') for i in foo]
            index   =   0
            while True :
                if foo[index][0] == 'l'   : 
                    l  =  (int(foo[index][1]))
                    index  += 1
                elif foo[index][0] == 'm' : 
                    if foo[index][1]  ==  '':
                        m   =   (int(foo[index+1][0]))
                    else :
                        m   =   (int(foo[index][1]))
                    index  += 1
                elif index  ==  2 :
                    break 
            for i in [l, m] :
                arr.append(i)
            states.append(arr)
            
    else : 
        print('Problem in definition read_states_projwfc_out()\
              for nspin value other than 1, 2 and 3')
        sys.exit(0)
           
    return(states)

def quantum_state(nspin, state, states):
    '''
    Finds the quantum state of an electron from the states list
    for nspin = 3 , i,e, non-collinear, spin-orbit case

    Parameters
    ----------
    nspin : Integer
        DESCRIPTION :   output of read_scf_out()
                        1 is  spin unpolarized
                        2 is  spin polarized 
                        3 is  non-collinear and spin-orbit 
                        
    state  :  Integer
        DESCRIPTION : State Index
        
    states :  List
        DESCRIPTION : Output of  read_states_projwfc_out()
        
    Returns
    -------
    Atom     Index as integer
    Orbital  Index as integer as in format of projwfc.x
             
                 For l=0 :
                  0 s             
                 For l=1:
                  1 p_z  ;  2 p_x  ;  3 p_y 
                  
                 For l=2:
                  4 d_z2  ;  5 d_zx  ;  6 d_zy  ;  7 d_(x2-y2)  ;  8 d_xy  
                  
                 For l=3:
                  9 f_z3  ;  10 f_zx2 ;  11 f_yz2 ;  12 f_z(x2-y2) ;  
                  13 f_xyz ;  14 f_x(x2-3y2)  ;  15  f_y(3x2-y2)
                  
                Further info at  
                https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html#idm89

    Spin     Index as integer for nspin = 3
             In the output of projwfc.x, first line of the 
             two orbital lines in each shell is 0 and the 
             second one is 1.
             
             For nspin = 1 or 2, spin value is provieded
             according to the file projwfc_up or projwfc_down
             in the later definitions
    
    '''
    atom     =   states[state-1][1] - 1  #  -1 because python index
    l        =   states[state-1][4]
    
    if nspin == 1 or nspin == 2 :
        m = states[state-1][5]
        
    elif nspin == 3 :
        j        =   states[state-1][5]
        m_j      =   states[state-1][6]
        
    else : 
        print('Problem in definition read_states_projwfc_out()\
              for nspin value other than 1, 2 and 3')
        sys.exit(0)
    
    
    
    norbitals   =   0                   #  For Counting number of
                                        #  orbitals for a given shell
    # Finding orbital and spin index for given state
    # Loop over all the shells of an atom, and count the number
    # of orbitals for that shell, e.g,
    # 4p and 5p -->
    # [0, 1, 2] and [3, 4, 5]
    for i in range(len(states)) :
        
        # Assuming that bands are always even numbered
        # Assigining spin index alternatively for nspin 3
        # Doesn't matters for nspin 1 or 2
        # because it will depend of file projwfc_up
        # or projwfc_down
        if  i  !=  2  :  
            spin  =  0         # Odd Spins are Up
        else :
            spin  =  1         # Even Spins are Down
        
        if states[i][1]  ==   (atom + 1)  :      
            if states[i][4]  == l :    
                
                if nspin ==  3  :              
                    if  states[i][5]  ==  j  :
                        if states[i][6]  ==  m_j  :
                            orbital_position  =  norbitals
                            
                elif nspin == 1 or nspin == 2 :
                    if states[i][5]  ==  m  :
                            orbital_position  =  norbitals
                
                else : 
                    print('Problem in definition quantum_states() \
                          for nspin value other than 1, 2 and 3')
                    sys.exit(0)
                
                norbitals   +=   1
                
    # Example of  below code :            
    # For collected form of 4p and 5p -->
    # [0, 1, 2] and [3, 4, 5]            
    # Now merge the orbital index like 3rd index 
    # of 2nd part is numbered as index 1 
    # thus summing all the orbitals of p shell
    # Although spin indices are still separately taken care
    # of for nspin = 3 in the above loop
                  
    if     l  ==  0 :
        orbital  =  0
        
    elif   l  ==  1 :
        orbital  =  orbital_position  %  3  +  1
        
    elif   l  ==  2 :
        orbital  =  orbital_position  %  5  +  4
        
    elif   l  ==  3 :
        orbital  =  orbital_position  %  7  +  9
        
    else  :  
        print('There are more than f orbitals')
        print('Please modify the definition quntum \
              states of this code')
        sys.exit(0)

    return(atom, orbital, spin)


def read_projwfc_new(filenames, states, nkpts, nbnd, nspin) :
    '''
    Creating a 5 Dimensional Array containng  Band Weights

    Parameters
    ----------
    filenames : List of String
        DESCRIPTION : Output of write_new_projwfc()    
        
    states : Integer
        DESCRIPTION : Output of read_states_projwfc_out()
        
    nkpts : Integer
        DESCRIPTION : Number of K-Points (output of 
                        read_nbnd_nkpts_projwfc_out()), manually specify
                        for previous versions of QE
                        
    nbnd : Integer
        DESCRIPTION : Number of Bands as in nscf output (output of
                        read_nbnd_nkpts_projwfc_out()), manually specify
                        for previous versions of QE
            
    nspin : Integer
        DESCRIPTION :   output of read_scf_out()
                        1 is  spin unpolarized
                        2 is  spin polarized 
                        3 is  non-collinear and spin-orbit 

    Returns
    -------
    Band weights in array format
    [nkpts, nbnd, spins, natoms, norbital]
    for all calculations,  spins = 2
    

    '''            
    natoms     =   max([states[i][1] for i in range(len(states))])
    nshell     =   max([states[i][4] for i in range(len(states))]) + 1
    norbitals  =   0
    for  i in range(nshell) :
        norbitals += 2 * i + 1
        
    band_weights = np.zeros([nkpts, nbnd, 2 , natoms, norbitals])
    scale   =   1
    
    # Reading Projwfc_up for nspin = 1, 2 or 3
    for file_index in range(len(filenames)) :
        start   =   0
        stop    =   (nkpts * nbnd)
        filename  =  filenames[file_index]
            
        projwfc = open(filename, 'r')       
        
        for line_index, lines in enumerate(projwfc) :
                        
            if line_index == start :
                #  Checking if we have reached the last line which is empty
                if not lines :  
                    # print('We have reached the last line of %s' %filename)
                    break
            
                state = int(lines.split()[0])
                atom, orbital, spin = quantum_state(nspin   =   nspin,
                                                    state   =   state,
                                                    states  =   states)
                # print(atom, orbital, spin)
            
            elif line_index > start and line_index <= stop :
    
                kpt     =   int(lines.split()[0]) - 1
                bnd     =   int(lines.split()[1]) - 1
                weight  =   float(lines.split()[2])
                
                if nspin == 1  or  nspin ==  2 :
                    spin    =     file_index
                    
                if file_index == 1 :
                    kpt  -=  nkpts
                
                band_weights[kpt, bnd, spin, atom, orbital] += weight * scale
                
                if line_index == stop :
                    
                    start +=  (nkpts * nbnd) + 1
                    stop  =  start + (nkpts * nbnd)
        
        projwfc.close()
        
    return(band_weights)


def read_scf_out(filename = 'scf.out'):
    '''
    Reads the Terminal Output of Scf
    Parameters
    ----------
    filename : String
        DESCRIPTION : Terminal Output of Scf

    Returns
    -------
    Fermi Energy as float in eV
    '''
    fermi_level  =  float(open(filename, 'r').read().split('the Fermi energy is')\
                    [1].split()[0])
        
    
    return(fermi_level)

def read_band_energies(filename):
    '''
    Read the total Band Structure from bands.x filband.gnu

    Parameters
    ----------
    filename : String
        DESCRIPTION :  Bands.x filband output in gnu format

    Returns
    -------    
    Total Band structure in Numpy Array format with 
    subtracted fermi energy, and 1st coloumn contains k-path
              band_1, band_2, band_3, ....
    kpoint_1   0.1     0.2     0.3
    Kpoint_2   0.4      0.5      0.6
    Kpoint_3   0.7      0.8      0.9
    ...
    
    '''
    gnu = open(filename + '.gnu', 'r').readlines()  
    
    efermi   =   read_scf_out()      #in eV
    print('The Fermi Level is at %10.4f  eV is subtracted in plots'\
          %efermi)
    
    grand_array = []
    array = []
    for lines in gnu :
        to_append  =   lines.split()
        if len(to_append) == 0 :
            grand_array.append(array)
            array = []
            continue
        
        else :
            array.append(to_append)
    
    # Format of below code is
    # [ [kpoint, band_1], [kpoint, band_2], ... ]
    bands   =    np.array(grand_array, dtype=float)
    
    num_bands    =    len(bands[:,0,0])
    num_k        =    len(bands[0,:,0])
    
    band_array   =    np.zeros([num_k, num_bands+1], dtype=float)

    for band_index in range(num_bands) :
        for k_index in range(num_k) :
            band_array[k_index, band_index + 1] = bands[band_index, 
                                                    k_index, 1] - efermi
            band_array[k_index, 0]   =   bands[0, k_index, 0]  

  
    return(band_array)


def read_kpoints(filename) :
    '''
    Read bands.in and extract kpoints informations

    Parameters
    ----------
    filename : String
        DESCRIPTION :  Terminal output of pw.x for bands.in

    Returns
    -------
    knames  :   Names of the High Symmetry K Points
    '''
    kpoints_   =   open(filename).read().split('K_POINTS')[1].split('\n')[2:]
    # Leaving all the blank strings
    kpoints    =   []
    for points in kpoints_  :
        if len(points.split()) > 2 :
            kpoints.append(points)
            
    knames    =   []
    for i in range(len(kpoints)) :
        if kpoints[i] :  # Checking  if the string is not empty
            kname   =   r'$\rm' + kpoints[i].split()[-1].split('!')[1] + '$'
            knames.append(kname)
            
    return(knames)

def read_bandsx_out(filename) :
    '''
    Read the k-ticks in reciprocal distance
    
    Parameters
    ----------
    filename : String
        DESCRIPTION  :   Terminal output of bands.x

    Returns
    -------
    K-ticks in reciprocal distance in list of floats
    '''
    bandsx_out   =   open(filename).read().split('x coordinate')[1:]                  
         
    kticks   =   []
    for lines in bandsx_out :
        data   =   float(lines.split()[0])
        kticks.append(data)
                
    return(kticks)


def normalize(array) :
    
    vmin       =    np.min(array)
    array     +=    vmin              # Shift up -ve to 0  if any
    vmin       =    np.min(array)   
    array     -=    vmin              # Shift down +ve min to 0
    vmax       =    np.max(array) 
    array     /=    vmax              # Scale maximum value to 1
    
    return(array)

def get_common_from_list(my_list) :
    '''
    Finds the reduced number of elements from a list
    
    Parameters
    ----------
    my_list : List
        DESCRIPTION   :  Input List

    Returns
    -------
    Input List with Reduced Items
    '''
    
    dict_list = collections.Counter()
    for i  in my_list :
        dict_list[i] += 1
        
    my_list   =    []
    for key in dict_list :
        my_list.append(key)
    
    return(my_list)

def adding_weights(band_weights, states, atoms_name, orbitals_name,
                   atoms=[-1], orbitals=[-1],
                   spins=[-1]
                   ) :
    '''
    Adds the weights for the given number of atoms, orbitals,
    spins, with the information of projections.

    Parameters
    ----------
    band_weights : 5-D array
        DESCRIPTION  :   output of read_projwfc_up_new()
        
    states   :   List 
        DESCRIPTION  :  Output of read_states_projwfc_out()
        
    atoms : List, optional
        DESCRIPTION  :  List of atoms to Project,
                        The default is [-1].
    orbitals : List, optional
        DESCRIPTION  :  List of Orbitals to Project.
                        The default is [-1].
    spins : List, optional
        DESCRIPTION  :  List of Spin states to Project.
                        The default is [-1].

    Returns
    -------
    Added Band Weights in 2D array of format
    [nkpoints,  nbands]
    List of Projections in format
    [atom_name,  orbital_state,  spin_state]
    Printing informations as Strings
    '''    
    spins_name   =   ['Spin Up', 'Spin down'] 
    
    orbitals_all    =   ['0 : s', 
            
                '1 : p_z', '2 : p_x', '3 : p_y',
            
                '4 : d_z2', '5 : d_zx', '6 : d_zy',
                '7 : d_(x2 - y2)', '8 : d_xy',
            
                '9 : f_z3', '10 : f_zx2', 
                '11 : f_yz2', '12 : f_z (x2 - y2)', 
                '13 : f_xyz', '14 : f_x(x2 - 3y2)', 
                '15 : f_y(3x2 - y2)'   ]
    
    data         =    band_weights[:,:,0,0,0] * 0
    projected    =    []
    for atom in atoms :
        for orbital in orbitals :            
            for spin in spins :
                
                projected.append([atoms_name[atom], 
                                  orbitals_all[orbital],
                                  spins_name[spin]])
                
                data   +=    band_weights[:,:, spin,
                                          atom, orbital]
                
    # To be in safe side for 1 atom, 1 orbital and 1 spin
    projected.append([' ', ' ', ' '])   

    towrite     =   '-' * 20 + '\n'
    towrite     +=   'Projections are done on  :\n'
    towrite     +=   'Atom       Orbital      Spin\n'
    towrite     +=   '-' * 20 + '\n'
    
    
    for i in range(len(projected)) :
        for j in range(len(projected[0])) :
            towrite   +=   '   %s   '  %(projected[i][j])
        towrite  +=  '\n'
                
    return(data, projected, towrite)

def max_atom_orbital(states) :
    '''
    Reads the Maximum number of Atoms and Orbitals in the 
    System
    Parameters
    ----------
    states : List
        DESCRIPTION  :   output of read_states_projwfc_out()

    Returns
    -------
    Names of Atoms and Orbitals with index as list

    '''
    
    # Naming the atoms and orbitals 
    atoms_index     =   [(states[i][1] - 1) for i in range(len(states))]    
    atoms_index     =    get_common_from_list(atoms_index)
    atoms_index.sort()
    
    atoms_name      =   []
    for i in atoms_index :
        for j in range(len(states)) :
            
            if i   ==  states[j][1] - 1 :
                atoms_name.append(' %s : %s '   %(i, states[j][2]))
        
    atoms_name   =   get_common_from_list(atoms_name)
                
    _          =   None
    shell      =   [states[i][4] for i in range(len(states))]
    shell      =   get_common_from_list(shell)
    shell.sort()
    
    
    orbitals_index   =   []
    for i in shell :
        
        if   i  ==  0   :             #  s shell
            orbitals_index.append(0)
            
        elif i  ==  1   :             #  p shell
             for j  in range(1,4)  :
               orbitals_index.append(j)
            
        elif i  ==  2   :             #  d shell
            for j  in range(4,9) :
                orbitals_index.append(j)  
            
        elif i  ==  3   :             #  f shell
            for j  in range(9,16) : 
                orbitals_index.append(j)  
            
        else  :  
            print('More than f orbitals Present')
            print('Please change Line number 477 in this code')
            sys.exit(0)
            
        
    orbitals_all    =   ['0 : s', 
                    
                        '1 : p_z', '2 : p_x', '3 : p_y',
                    
                        '4 : d_z2', '5 : d_zx', '6 : d_zy',
                        '7 : d_(x2 - y2)', '8 : d_xy',
                    
                        '9 : f_z3', '10 : f_zx2', 
                        '11 : f_yz2', '12 : f_z (x2 - y2)', 
                        '13 : f_xyz', '14 : f_x(x2 - 3y2)', 
                        '15 : f_y(3x2 - y2)'   ]
    
    orbitals_name   =   [orbitals_all[i] for i in orbitals_index] 
    
    towrite      =   'Maximum Atoms and Orbitals :\n['
    for atom in atoms_name :
        towrite +=   ' %s ,'  %atom
    towrite     +=    '] Total : %d\n['  %len(atoms_name)
    
    for orbital in orbitals_name :
        towrite +=   ' %s ,'  %orbital
        
    towrite     +=   '] Total : %d\n'  %(len(orbitals_name))
    
    return(atoms_name, orbitals_name, towrite)


def Total_Plot(x, ys, figsize, title, linewidth, 
               dpi, color, window_on, 
               window_lim,  win_color,
                xticks            =     None, 
                yticks            =     None, 
                xnames            =     None
                ) :
    '''
    Gives only the Total Band Structure as Plot
    
    Parameters
    ----------
    x : Array
        DESCRIPTION   :  Kpath
        
    ys : Array
        DESCRIPTION   :  Bands
        
    window_on  :  True  or False 
    DESCRIPTION  :  Wheather to plot wannier window
                    Default : False
        
    window_lim  :  List of lists
        DESCRIPTION  :  Containts number  of lists like :
                        [ [min1, max1], [min2, max2], ...]
        
    window_color  :  List of strings
        DESCRIPTION : Contains the list of colour of  window
                      ['g', 'b']  default : 'g'for  Green
        
    figsize :  List with two integers
        DESCRIPTION   :  Decides the Size of the Figure for same dpi
        
    title : String
        DESCRIPTION   :  Title of the Graph
        
    linewidth : Integer or Float
        DESCRIPTION   :  Decides the width of each y
        
    dpi : Integer
        DESCRIPTION   :  Number of pixels per inch
        
    xticks : List, optional
        DESCRIPTION   :  Ticks at K path 
                         The default is None.
                         
    yticks : List, optional
        DESCRIPTION   :  Ticks in the y-axis
                         The default is None.
                         
    xnames : List of Strings, optional
        DESCRIPTION   :  Names at xticks
                         The default is None.
                         
    color : String, optional
        DESCRIPTION   :  Colour of the ys
                         The default is 'afmhot_r'.
                        
    Returns
    -------
    Plots as matplotlib figure and axes
    '''
    # Plotting
    fig, axs   =    plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    factor     =    0.05 * abs( np.max(ys) - np.min(ys) )
    yrange     =    [int(np.min(ys) - factor), 
                     int(np.max(ys) + factor)]
    xrange     =     [min(x), max(x)]
    # if xticks == None   :   xticks     =    np.linspace(np.min(x), 
    #                                                     np.max(x), 3) 
    if yticks == None   :   yticks     =    np.linspace(yrange[0], 
                                                        yrange[1], 3) 
            
    axs.plot(x , ys, color=color, lw=linewidth)
    
    labelsize  =   max(figsize) * 3
    axs.set_xlim(xrange[0], xrange[1])
    axs.set_ylim(yrange[0], yrange[1])
    axs.set_title(title, fontsize=labelsize)
    axs.set_xticks(xticks)
    axs.set_yticks(yticks)
    if xnames !=  None   :   axs.set_xticklabels(xnames, 
                                                 rotation='horizontal')
    
    for xtick in xticks :
        axs.axvline(x=xtick, ymin=min(yrange),
                ymax=max(yrange), lw=linewidth*0.8, 
                color='k', alpha=1)
    axs.axhline(y=0, xmin=xrange[0], xmax=xrange[1], 
              linestyle='--',  lw=linewidth*1.2,
              color='r', alpha=1)
    
    # Plotting Window
    fig, axs    =      window(fig          =     fig,
                              axs          =     axs,
                              window_on    =     window_on,
                              window       =     window_lim,
                              color        =     win_color,
                              labelsize    =     labelsize*1.2
                              )   
    
    plt.xticks(fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    plt.ylabel(ylabel=r'$E$ - $E_{f}$  (eV)',  fontsize=labelsize*1.2)
    
    return(fig, axs)

def Projected_Plot(x, ys, weights, figsize, title, 
                   linewidth, dpi, window_on, 
                   window_lim, win_color,
                    xticks            =     None, 
                    yticks            =     None, 
                    xnames            =     None, 
                    cmap              =     'afmhot_r',
                    clr_tick_space    =      3
                    ) :
    '''
    Gives Projected Band Structure as Plot
    
    Parameters
    ----------
    x : Array
        DESCRIPTION   :  Kpath
        
    ys : Array
        DESCRIPTION   :  Bands
        
    weights  :  2-D array
        DESCRIPTION : Contains the indivisual or added 
                      added band weights from adding_weights()
                      
    window_on  :  True  or False 
        DESCRIPTION  :  Wheather to plot wannier window
                        Default : False
        
    window_lim  :  List of lists
        DESCRIPTION  :  Containts number  of lists like :
                        [ [min1, max1], [min2, max2], ...]
        
    window_color  :  List of strings
        DESCRIPTION : Contains the list of colour of  window
                      ['g', 'b']  default : 'g'for  Green
        
    figsize :  List with two integers
        DESCRIPTION   :  Decides the Size of the Figure for same dpi
        
    title : String
        DESCRIPTION   :  Title of the Graph
        
    linewidth : Integer or Float
        DESCRIPTION   :  Decides the width of each y
        
    dpi : Integer
        DESCRIPTION   :  Number of pixels per inch
        
    xticks : List, optional
        DESCRIPTION   :  Ticks at K path 
                         The default is None.
                         
    yticks : List, optional
        DESCRIPTION   :  Ticks in the y-axis
                         The default is None.
                         
    xnames : List of Strings, optional
        DESCRIPTION   :  Names at xticks
                         The default is None.
                         
    color : String, optional
        DESCRIPTION   :  Colour of the ys
                         The default is 'afmhot_r'.

    Returns
    -------
    Plots as matplotlib figure and axes
    '''
    
    
    
    ys         =    np.transpose(ys)
    weights    =    np.transpose(weights)

    # Plotting
    fig, axs   =    plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    factor     =    0.05 * abs( np.max(ys) - np.min(ys) )
    yrange     =    [int(np.min(ys) - factor), 
                     int(np.max(ys) + factor)]
    # yrange     =    [int(round(int(np.min(ys) - factor)/10))*10, 
    #                  int(round(int(np.max(ys) + factor)/10))*10]
    xrange     =     [min(x), max(x)]
    # if xticks == None   :   xticks     =    np.linspace(np.min(x), 
    #                                                     np.max(x), 3) 
    if yticks == None   :   yticks     =    np.linspace(yrange[0], 
                                                        yrange[1], 3) 

    # Normalize the weights from 0 to 1 for relative scaling
    # weights    =    normalize(weights)
    # vmin       =    np.min(weights)
    # vmax       =    np.max(weights)
    clr_ticks  =    []#np.linspace(vmin, vmax, clr_tick_space)
    # norm       =    plt.Normalize(vmin, vmax)

    labelsize  =   max(figsize) * 3
    axs.set_xlim(xrange[0], xrange[1])
    axs.set_ylim(yrange[0], yrange[1])
    axs.set_title(title, fontsize=labelsize)

    for i in range(len(ys)) :
        y          =   ys[i]
        points     =   np.array([x, y]).T.reshape(-1, 1, 2)
        segments   =   np.concatenate([points[:-1], points[1:]], axis=1)


        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, cmap=cmap)#, norm=norm)
        # Set the values used for colormapping
        lc.set_array(weights[i])
        lc.set_linewidth(linewidth)
        line = axs.add_collection(lc)

    cb         =   fig.colorbar(line, 
                                ax=axs, ticks=clr_ticks, 
                                )
    cb.ax.tick_params(labelsize=labelsize)
    axs.set_xticks(xticks)
    axs.set_yticks(yticks)
    if xnames !=  None   :   axs.set_xticklabels(xnames, 
                                                 rotation='horizontal')
    
    for xtick in xticks :
        axs.axvline(x=xtick, ymin=min(yrange),
                ymax=max(yrange), lw=linewidth*0.8, 
                color='k', alpha=1)
    axs.axhline(y=0, xmin=xrange[0], xmax=xrange[1], 
              linestyle='--',  lw=linewidth*1.2,
              color='r', alpha=1)
    
    # Plotting Window
    fig, axs    =      window(fig          =     fig,
                              axs          =     axs,
                              window_on    =     window_on,
                              window       =     window_lim,
                              color        =     win_color,
                              labelsize    =     labelsize*1.2
                              )

    plt.xticks(fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    plt.ylabel(ylabel=r'$E$ - $E_{f}$  (eV)',
               fontsize=labelsize*1.2) 

    return(fig, axs)


def project_bands(filband        =     'bands.dat',
                  projwfc_out    =     'projwfc.out',
                  filproj        =     'bands',
                  bands_in       =     'bands.in',
                  bandsx_out     =     'bandsx.out',
                  atoms          =      [-1],
                  orbitals       =      [-1],
                  spins          =      [-1],
                  cmap           =      'afmhot_r',
                  figsize        =      [20,16],
                  linewidth      =       5,
                  dpi            =       100,
                  title          =      'Projected Band Structure',
                  save_fig       =      'Bands.png',
                  band_weights   =       np.array([0]),
                  band_array     =       np.array([0]),
                  atoms_name     =       None, 
                  orbitals_name  =       None,
                  projection     =       True,
                  color          =       'b',
                  window_on      =       False,
                  window_lim     =       None,
                  win_color      =       'lime',
                  matlab         =       False,
                  nspin          =       1
                  ) :
    
    
    #----------Reading and Aranging Datas from Ouput--------
    nspin      =   nspin
    states     =   read_states_projwfc_out(projwfc_out, nspin)
    
    if len(band_array)  ==  1 : 
        band_array      =   read_band_energies(filband)
        atoms_name, orbitals_name    =   max_atom_orbital(states = states)[:2]
        
    bands      =   band_array[:,1:]
    kpath      =   band_array[:,0]
    nkpts      =   len(bands[:,0])
    nbnd       =   len(bands[0,:])
    

    towrite = max_atom_orbital(states)[2]
    with open('states.txt', 'w') as file :
        file.write(towrite)
        
    if len(band_weights)  ==  1 :
        filproj         =   write_new_projwfc(filproj, nspin)
        band_weights    =   read_projwfc_new(filproj, 
                                            states, nkpts, 
                                            nbnd, nspin)
            
    band_weights_continue  =  band_weights
        
        
    knames          =   read_kpoints(bands_in) 
    kticks          =   read_bandsx_out(bandsx_out)

    plt.ioff()    
    #------------Plotting--------------------    
    if projection  ==  True  :
        
        band_weights, projected, towrite \
                       = adding_weights(band_weights   =   band_weights,
                                        states         =   states,
                                        atoms_name     =   atoms_name, 
                                        orbitals_name  =   orbitals_name,
                                        atoms          =   atoms,
                                        orbitals       =   orbitals,
                                        spins          =   spins
                                        )
        
        fig, ax       =  Projected_Plot(  x              =   kpath, 
                                          ys             =   bands, 
                                          weights        =   band_weights,
                                          title          =   title,
                                          linewidth      =   linewidth,
                                          dpi            =   dpi,
                                          cmap           =   cmap,
                                          xticks         =   kticks,
                                          xnames         =   knames,
                                          figsize        =   figsize,
                                          window_on      =   window_on,
                                          window_lim     =   window_lim,
                                          win_color      =   win_color 
                                          )
        
    elif projection  ==  False :
        fig, ax       =  Total_Plot(x              =   kpath, 
                                    ys             =   bands, 
                                    figsize        =   figsize,
                                    title          =   title,
                                    linewidth      =   linewidth,
                                    dpi            =   dpi,
                                    color          =   color,
                                    xticks         =   kticks,
                                    xnames         =   knames,
                                    window_on      =   window_on,
                                    window_lim     =   window_lim,
                                    win_color      =   win_color 
                                    )
                
    print(towrite)   
    fig.savefig(save_fig)
    
    if matlab  ==  True :
        MATLAB_output(band_energies   =    bands, 
                      band_weights    =    band_weights_continue, 
                      kpath           =    kpath,
                      kticks          =    kticks,
                      knames          =    knames      
                      )
        
    return(fig, ax, band_weights_continue, band_array, 
           atoms_name, orbitals_name)


def window(fig, axs, window_on, window, labelsize, color,
           interval         =    0.05, 
           transparency    =    0.06
           ) :
    '''
    Shows a  transparent Window on the Plot
    Very useful to study disentanglement windows in Wannier90

    Parameters
    ----------
    fig : Matplotlib Figure object
        DESCRIPTION  :  Figure Contianing the Band Structure
        
    axs : Matplotlib Axes object
        DESCRIPTION  :  Figure Contianing the Band Structure
        
    window_on :  True or False
        DESCRIPTION  :  Says weather wannier window should be on
        
    window : List with two Contents of Float type
        DESCRIPTION  :  Decide the minimum and Maximum of the window
        
    Color  :  String
        DESCRIPTION  :   Decides the colour of the Window
    labelsize    :   Integer or Float
        DESCRIPTION  :  Font size for x and y Labels
        
    ticks   :   Numpy Array of Integer and Float type
        DESCRIPTION  :   ticks on y-axis
        
    Returns
    -------
    Matplotlib Figure and Axes Object.
    '''
    
    
    # Plotting Window
    if window_on  == True   :
        for i  in range(len(window)) :
            print('Wannier Window is On with Range : ')
            print(window[i])
            xlim         =         axs.get_xlim()
            yticks       =         axs.get_yticks()
            for j in range(len(window[i]))  :
                np.append(yticks,  window[i][j])
            win_array    =         np.arange(window[i][0], 
                                             window[i][1], 
                                             interval)   
            for hline in win_array  :
                axs.axhline(y      =     hline, 
                           xmin    =     xlim[0], 
                           xmax    =     xlim[1],  
                           color   =     color[i],
                           alpha   =     transparency,
                           lw      =     20
                           )
            for hline in [win_array[0], win_array[-1]]  :
                axs.axhline(y      =     hline, 
                           xmin    =     xlim[0], 
                           xmax    =     xlim[1],  
                           color   =     color[i],
                           lw      =     5
                           )
        axs.set_yticks(ticks  =   yticks)  
          
    else   :
        print('Wannier Window is Off')
        
    return(fig, axs)

def MATLAB_output(band_energies, band_weights, kpath, kticks, knames) :
    '''
    Generates the Matlab Ouputs 

    Parameters
    ----------
    band_energies : 2-D Numpy Array 
        DESCRIPTION   :  Contains energies in Dimension 
                        [kpoints, nbands]
    band_weights : 5-D Numpy Array
        DESCRIPTION   :   Contains all the projection inforamation
                          as in output of read_projwfc_new()
    kpath : 1-D Numpy Array
        DESCRIPTION   :   Contains the Distance of the k-points
    kticks : 1-D Numpy Array
        DESCRIPTION   :   Contains the position of High Symmetry points
    knames : List of Strings
        DESCRIPTION   :   Contains the names of High Symmetry points

    Returns
    -------
    Writes the arrays into .mat files in MATLAB_OUTPUT Folder.

    '''
    
    cwd  =  os.getcwd()
    folder = 'MATLAB_OUTPUT'
    os.chdir(cwd)
    try         :     
        os.mkdir(folder)
    except      :
        shutil.rmtree(folder)
        os.mkdir(folder)
        
    print('.MAT FILES IN MATLAB_OUTPUT FOLDER')    
    os.chdir(folder)
    
    scipy.io.savemat("BAND_ENERGIES.mat", 
                     {'kpoint_band_spin' : band_energies })
    scipy.io.savemat("BAND_WEIGHTS.mat", 
                     {'kpoint_band_spin_atom_orbital' : band_weights} )
    scipy.io.savemat("KPOINTS.mat", {'kpoints' : kpath})
    scipy.io.savemat("KTICKS.mat", {"kticks" : kticks})
    scipy.io.savemat("KNAMES.mat", {'knames' : knames} )
    
    os.chdir(cwd)
    
    return(1)

#-------------------------------------------#
#                Initialization             #
#-------------------------------------------#

# Plotting Total Band Structure

atoms       =        [0,1,2,3,4,5]
orbitals    =        [0, 1, 2, 3, 4, 5, 6, 7, 8]
spins       =        [0,1]
color       =        'b'   
savefile    =        'Total_Band_Structure.png'
yticks      =        [-6, 0, 6]
window_on   =        False
window_frz  =        [-1.6, 3]
frz_color   =        'lime'
window_dis  =        [3,    5]
dis_color   =        'aqua'
window_lim  =        [window_frz, window_dis]
win_color   =        [frz_color, dis_color] 
nspin       =        2


fig, ax, band_weights, band_array, atoms_name, orbitals_name =\
                  project_bands( 
                                atoms          =     atoms,
                                orbitals       =     orbitals,
                                spins          =     spins,
                                title          =     'Total Band Structure',
                                save_fig       =     savefile,
                                projection     =     False,
                                color          =     color,
                                window_on      =     window_on,
                                window_lim     =     window_lim,
                                win_color      =     win_color,
                                matlab         =     True,
                                nspin          =     nspin
                                )
                 
ax.set_yticks(ticks=yticks)
ax.set_ylim(yticks[0], yticks[2])
fig.savefig(savefile)


# Plotting Projected Band Structure

atoms_list       =    [[0]]
atom_names       =    ['Fe']
orbitals_list    =    [[0], [1,2,3], [4,5,6,7,8]]
orbital_names    =    ['s', 'p', 'd']
spins_list       =    [[0], [1]]
spin_names       =    ['Up', 'Down']
cmap             =    'viridis_r'


for i in range(len(atoms_list)) :
    atom        =     atom_names[i]
    atoms       =     atoms_list[i]
    
    for j in range(len(orbitals_list)) :
        orbital    =  orbital_names[j]
        orbitals   =   orbitals_list[j]
        
        for k in range(len(spins_list)) :
            spin     =   spin_names[k]
            spins    =   spins_list[k]
        
            title         =   'Atom : %s   Orbital : %s  Spin : %s\n'  \
                                %(atom, orbital, spin)
            savefile      =   'Projection_%s_%s_%s.png'    \
                                %(atom, orbital, spin)
            
            fig, ax = project_bands(  
                                    atoms          =   atoms,
                                    orbitals       =   orbitals,
                                    spins          =   spins,
                                    title          =   title,
                                    save_fig       =   savefile,
                                    cmap           =   cmap,
                                    band_weights   =   band_weights,
                                    band_array     =   band_array,
                                    atoms_name     =   atoms_name, 
                                    orbitals_name  =   orbitals_name,
                                    window_on      =   window_on,
                                    window_lim     =   window_lim,
                                    win_color      =   win_color,
                                    nspin          =   nspin
                                    )[:2]
        
            ax.set_yticks(ticks=yticks)
            ax.set_ylim(yticks[0], yticks[2])
            fig.savefig(savefile)


band_weights  =  None   # Clear Memory

#-----------Rough work----------------------
