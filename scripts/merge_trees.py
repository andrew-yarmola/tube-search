#!/usr/bin/env python3

import os, subprocess, sys, getopt, glob, time, re
from time import sleep
from multiprocessing import Process

# A few useful functions
def command_output(command):
    try:
        # Note that subprocess.check_output retuns a byte string
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        out, err = pipe.communicate()
        return out.decode('utf8')
    except:
        print('Error in {0}\n'.format(command)) 
        sys.exit(1)

def add_holes(holes, treeholes, directory, boxfile):
    command = '{0} -r {1} \'{2}\''.format(treeholes, directory, boxfile)
    out = command_output(command)
    if not out: return
    out = out.replace('root', '')
    if len(out) == 0 : out = 'root'
    new_holes = set(out.rstrip().split('\n'))
    holes |= set(h for h in new_holes)

def run_merge(command, dest_dir) :
    pid = os.getpid()
    pid_file = dest_dir + '/' + str(pid) + '.pid'
    with open(pid_file,'a') as fp : 
        fp.write(command + '\n')
        fp.close()
    return_code = subprocess.call(command, shell=True)
    if return_code == 0:
        with open(pid_file,'a') as fp :
            fp.write('\ncompleted')
        sys.exit(0)
    else:
        with open(pid_file,'a') as fp :
            fp.write('\nfailed')
        sys.exit(1)

if __name__ == '__main__' :

    try:
        opts, args = getopt.getopt(sys.argv[1:],'s:d:c:',['size=','depth=','child_limit='])
    except getopt.GetoptError as err:
        print(str(err))
        print('Usage: merge_tree [-s,--size <max_subtree_size>] ' +
              '[-d,--depth <min_depth>] [-c,--child_limit <limit>] dest_dir src_dir(s)')
        sys.exit(2)

    if len(args) < 2:
        print('Usage: merge_tree [-s,--size <max_subtree_size>] ' +
              '[-d,--depth <min_depth>] [-c,--child_limit <limit>] dest_dir src_dir(s)')
        sys.exit(2)

    # Executables
    treemerge = './treemerge'
    treecheck = './treecat --mark -s'
    treeholes = './treecat --open_holes'

    # Set up the rest of the arguments
    dest_dir = args[0]
    src_dirs = ' '.join(args[1:])
    size = '100000'
    depth = '6'
    child_limit = 8
    depth_limit = 200

    # Get config
    for opt, val in opts:
        if opt in ('-s', '--size'):
            size = int(val)
        if opt in ('-d', '--depth'):
            depth = int(val)
        if opt in ('-c', '--child_limit'):
            child_limit = int(val)
            
    # Get holes. Note, treecat will check that all files are complete trees
    holes = set();
    add_holes(holes, treeholes, dest_dir, '')

    # Get done words
    done = set()
    try:
        done = set([os.path.basename(boxfile).replace('.out','')
            for boxfile in glob.glob(dest_dir + '/*.out')])
    except:
        print('Error reading {0}\n'.format(dest_dir)) 
        sys.exit(1)

    print("Launching Merge")

    # Launch the refine runs
    active_pid_to_hole = {};
    run_count = 0
    child_count = 0
    wait_for_holes = False
    while True:
        sleep(0.01) # We don't need to to run the main loop to death
        open_holes = holes - done
        best_hole = '1'*1000
        for hole in open_holes:
            if len(hole) < len(best_hole) or best_hole == '1'*1000 :
                best_hole = hole 
        if ((len(open_holes) == 0 and run_count == 0 and len(done) == 0)
            or len(best_hole) == 0):
            best_hole = 'root'

        if len(best_hole) > depth_limit:
            if child_count > 0 :
                wait_for_holes = True
            else :
                # We only break if we don't have any more refine processes running
                break
        else :
            wait_for_holes = False

        # We now check for completed refine processes.
        if (child_count >= child_limit or (child_count > 0
            and len(open_holes) == 0) or wait_for_holes):
            iter_dict = dict(active_pid_to_hole)
            for done_pid, done_hole in iter_dict.items() :
                pid_file = dest_dir + '/' + str(done_pid) + '.pid'
                status = command_output('tail -1 {0}'.format(pid_file))
                if 'completed' in status :
                    # We should check the output either way to make sure it is clean 
                    subprocess.call('{0} {1} \'{2}\''.format(
                        treecheck, dest_dir, done_hole), shell=True)

                    print('Completed {0} {1}\n'.format(done_hole,done_pid))
                    add_holes(holes, treeholes, dest_dir, done_hole)

                    child_count -= 1
                    del active_pid_to_hole[done_pid]
                    os.remove(pid_file)
                    continue

                elif 'failed' in status :
                    # We should check the output either way to make sure it is clean 
                    subprocess.call('{0} {1} \'{2}\''.format(treecheck, dest_dir, done_hole), shell=True)
                    # If there was an error refining
                    print('Error with pid {0} and hole {1}\n'.format(done_pid, done_hole))
                    done.remove(done_hole)
                    child_count -= 1
                    del active_pid_to_hole[done_pid]
                    os.remove(pid_file)
                    continue
                else :
                    continue
            sleep(0.01) # We don't need to to run the main loop to death
            continue        

        # If we make it here. We are running refine
        print('Open hole count: {0}\n'.format(len(open_holes)))
        print('Best hole: {0}\n'.format(best_hole))

        out = dest_dir + '/' + best_hole + '.out'
        err = dest_dir + '/' + best_hole + '.err'

        if best_hole == 'root':
          best_hole = '""'

        command = '{0} {1} {2} {3} {4} > {5} 2> {6}'.format(
            treemerge, size, depth, best_hole, src_dirs, out, err)

        print('Running with count {1}: {0}\n'.format(command, run_count))
        merge = Process(target=run_merge, args=(command, dest_dir,))
        merge.start()
        pid = merge.pid

        child_count += 1   
        run_count += 1
        done.add(best_hole)
        active_pid_to_hole[pid] = best_hole
