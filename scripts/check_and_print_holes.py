#!/usr/bin/python

import os, subprocess, sys, getopt, glob, time, re
from time import sleep
from multiprocessing import Process

# A few useful functions
def command_output(command):
    try:
        # Note that subprocess.check_output retuns a byte string
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        out, err = pipe.communicate()
        return out
    except:
        print('Error in {0}\n'.format(command)) 
        sys.exit(1)

def add_holes(holes, treeholes, directory, boxfile):
    command = '{0} -r {1} \'{2}\''.format(treeholes, directory, boxfile)
    byte_string = command_output(command)
    if not byte_string: return
    byte_string = byte_string.replace('root','')
    if len(byte_string) == 0 : byte_string = 'root'
    new_holes = set(byte_string.rstrip().split('\n'))
    holes |= set(h for h in new_holes if len(h) < 121)

def add_holes_from_file(holes,fp) :
    try:
        for line in open(fp):
            hole = line.rstrip()
            if hole[0] == '1' or\
               hole[0] == '0' :
                if len(hole) < 81 :  
                    holes.add(hole)
    except:
        print('Error loading holes file {0}\n'.format(fp))

def run_refine(command, dest_dir) :
    pid = os.getpid()
    pid_file = dest_dir + '/' + str(pid) + '.pid'
    with open(pid_file,'a') as fp : 
        fp.write(command + '\n')
        fp.close()
    returnCode = subprocess.call(command, shell=True)
    if returnCode == 0:
        with open(pid_file,'a') as fp :
            fp.write('\ncompleted')
        sys.exit(0)
    else:
        with open(pid_file,'a') as fp :
            fp.write('\nfailed')
        sys.exit(1)

if __name__ == '__main__' :

    args = sys.argv[1:]
    if len(args) != 1:
        print('Usage: check_and_print_holes tree_dir')
        sys.exit(2)

    # Executables
    treecat = '/u/yarmola/margulis-search/bin/treecat'
    treeholes = '/u/yarmola/margulis-search/bin/treecat --open_holes'
    treecheck = '/u/yarmola/margulis-search/bin/treecat --mark -s'

    # Set up the rest of the arguments
    tree_dir = args[0]

    # Check for incomplete trees
    subprocess.call('{0} -r {1} \'{2}\''.format(treecheck, tree_dir, ''), shell=True)

    # Get holes. Note, treecat will check that all files are complete trees
    holes = set();
    add_holes(holes, treeholes, tree_dir, '')

    for h in holes:
      print("{}\n".format(h))
