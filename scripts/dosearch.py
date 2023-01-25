#!/usr/bin/env python3

import os, subprocess, sys, getopt, glob, time, re
from time import sleep
from multiprocessing import Process

# A few useful functions
def command_output(command):
  try:
    # Note that subprocess.check_output retuns a byte string
    pipe = subprocess.Popen(command,
        stdout=subprocess.PIPE, shell=True)
    out, err = pipe.communicate()
    return out.decode('utf8')
  except:
    print('Error in {0}\n'.format(command)) 
    sys.exit(1)

def add_holes(holes, treeholes, directory, boxfile):
  command = '{0} -r {1} \'{2}\''.format(
      treeholes, directory, boxfile)
  out_string = command_output(command)
  if not out_string: return
  out_string = out_string.replace('root','')
  if len(out_string) == 0 : out_string = 'root'
  new_holes = set(out_string.rstrip().split('\n'))
  holes |= set(h for h in new_holes if len(h) < 145)

def add_holes_from_file(holes,fp):
  try:
    for line in open(fp):
      hole = line.rstrip()
      if (hole[0] == '1' or hole[0] == '0') and len(hole) < 200:  
        holes.add(hole)
  except:
    print('Error loading holes file {0}\n'.format(fp))

def add_words(words, fp):
  try:
    for line in open(fp):
      word = line.rstrip()
      if word[0].isdigit() or word[0] == 'X' or word[0] == 'H':
        continue
      else:
        if ',' in word :
          word = re.findall('\((.*?),.*\)', word)[-1]
        elif '(' in word :
          word = re.findall('\((.*?)\)', word)[-1]
        words.add(word)
  except:
    print('Error loading words file {0}\n'.format(fp))

def run_refine(command, dest_dir) :
  pid = os.getpid()
  pid_file = dest_dir + '/' + str(pid) + '.pid'
  with open(pid_file,'a') as fp : 
    fp.write(command + '\n')
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
    opts, args = getopt.getopt(sys.argv[1:],
        'w:p:b:c:d:h:r:i:t:s:IM:R:n:',
      ['words=','impossible=', 'child_limit=','depth_limit=',
        'holes=','refine=','invent_depth=',
        'truncate_depth','word_search_depth=', 'improve',
        'cosh_marg=', 'sinh_rad=', 'name='])
  except getopt.GetoptError as err:
    print(str(err))
    args = []

  if len(args) != 2:
    print('Usage: dosearch [-w,--words <words_file>] ' + 
      '[-p,--impossible <impossible_file>] ' +
      '[-b,--bad_relators <bad_relator_file>] ' +
      '[-c,--child_limit <limit>] [-d,--depth_limit <limit>]' +
      '[-h,--holes <holes_file>] [-I,--improve]' + 
      '[-M,--cosh_marg <limit>] [-R,--sinh_rad <limit>]' +
      '[-n,--name] <name> src_dir dest_dir')
    sys.exit(2)

  # Executables
  treecat = './treecat'
  treeholes = './treecat --open_holes'
  treecheck = './treecat --mark -s'
  refine = './refine_marg'

  # Set up the rest of the arguments
  src_dir = args[0]
  dest_dir = args[1]
  child_limit = 8
  depth_limit = 330

  cosh_mu_upper = '1.32'
  sinh_tube_upper = '1.48'

  max_size = '500000'
  max_depth = '150'
  truncate_depth = '6'
  invent_depth = '42'
  word_search_depth = '60'
  fill_holes = ''
  improve_tree = ''
  impossible_file = 'none'
  words_file = 'none'
  bad_relator_file = 'none'
  name = 'none'

  # Get config
  holes_file = None
  seen_words = set()
  for opt, val in opts:
    if opt in ('-w', '--words'):
      words_file = val
    if opt in ('-p', '--impossible'):
      impossible_file = val
    if opt in ('-b', '--bad_relators'):
      bad_relator_file = val
    if opt in ('-c', '--child_limit'):
      child_limit = int(val)
    if opt in ('-d', '--depth_limit'):
      depth_limit = int(val)
    if opt in ('-h', '--holes'):
      holes_file = val
    if opt in ('-r', '--refine'):
      refine = val
    if opt in ('-i', '--invent_depth'):
      invent_depth = val
    if opt in ('-t', '--truncate_depth'):
      truncate_depth = val
    if opt in ('-s', '--word_search_depth'):
      word_search_depth = val
    if opt in ('-I', '--improve'):
      improve_tree = ' --improve_tree'
    if opt in ('-M', '--cosh_marg'):
      cosh_mu_upper = val
    if opt in ('-R', '--sinh_rad'):
      sinh_tube_upper = val
    if opt in ('-n', '--name'):
      name = val

  add_words(seen_words, words_file)

  # Check for incomplete trees
  subprocess.call('{0} -r {1} \'{2}\''.format(
    treecheck, dest_dir, ''), shell=True)

  # Get holes. Note, treecat will check that all 
  # files are complete trees
  holes = set();
  if holes_file :
    add_holes_from_file(holes, holes_file)
  else :
    add_holes(holes, treeholes, dest_dir, '')

  # Get done holes
  done = set()
  try:
    done = set([os.path.basename(boxfile).replace('.out','')
      for boxfile in glob.glob(dest_dir + '/*.out')])
  except:
    print('Error reading {0}\n'.format(dest_dir)) 
    sys.exit(1)

  print(f"Launching Refine: {len(holes)} holes and {len(done)} done files")

  # Launch the refine runs
  active_pid_to_hole = {};
  refine_run_count = 0
  child_count = 0
  wait_for_holes = False
  failed_holes = set()
  while True:
    # We don't need to to run the main loop to death
    # since we aren't using os.wait
    sleep(0.001)
    open_holes = holes - done
    best_hole = '1'*400
    if len(open_holes) == 0 and refine_run_count == 0 and len(done) == 0:
      best_hole = 'root'
    elif len(open_holes) > 0:
      def len_and_lex(x):
        return f"{len(x):10b}" + x
      best_hole = min(open_holes, key=len_and_lex)

    if len(best_hole) > depth_limit:
      if child_count > 0 :
        wait_for_holes = True
      else :
        # We only break if we don't have any more refine processes running
        break
    else :
      wait_for_holes = False

    # We now check for completed refine processes.
    if (child_count >= child_limit or (child_count > 0 and 
        len(open_holes) == 0) or wait_for_holes):
      iter_dict = dict(active_pid_to_hole)
      for done_pid, done_hole in iter_dict.items():
        pid_file = dest_dir + '/' + str(done_pid) + '.pid'
        status = command_output('tail -1 {0}'.format(pid_file))
        if 'completed' in status :
          # We should check the output either way to make sure it is clean 
          subprocess.call('{0} {1} \'{2}\''.format(
            treecheck, dest_dir, done_hole), shell=True)

          print('Completed {0} {1}\n'.format(done_hole,done_pid))
          add_holes(holes, treeholes, dest_dir, done_hole)

          num_patched = command_output(
              'grep -c Patched {0}/{1}.err; exit 0'.format(
                dest_dir, done_hole)).rstrip()
          num_unpatched = command_output(
              'grep -c Unpatched {0}/{1}.err; exit 0'.format(
                dest_dir, done_hole)).rstrip()
          num_holes = command_output(
              'if [ -f {0}/{1}.out ]; then grep -c HOLE {0}/{1}.out; else echo 0; fi; exit 0'.format(
                dest_dir, done_hole)).rstrip()

          print('Holes: {0} patched, {1} unpatched, {2} open holes\n'.format(
            num_patched, num_unpatched, int(num_holes)))

          box_words = set()
          add_words(box_words, '{0}/{1}.out'.format(dest_dir, done_hole))        
          new_words = box_words - seen_words
          seen_words |= new_words

          bad_holes = command_output(
              'grep HOLE {0}/{1}.err | cut -d " " -f 2; exit 0'.format(
                dest_dir, done_hole)).rstrip().split('\n')
          failed_holes.update(bad_holes)

          if len(new_words) > 0: 
            with open(words_file, 'a') as f:
              for word in new_words:
                print('Adding word {0}'.format(word))
                f.write('(' + word + ',)' + '\n')

          child_count -= 1
          del active_pid_to_hole[done_pid]
          os.remove(pid_file)
          continue

        elif 'failed' in status :
          # We should check the output either way to make sure it is clean 
          subprocess.call('{0} {1} \'{2}\''.format(
            treecheck, dest_dir, done_hole), shell=True)
          # If there was an error refining
          print('Error with pid {0}\n'.format(done_pid))
          print('Error refining hole {0}\n'.format(done_hole))
          done.remove(done_hole)
          child_count -= 1
          del active_pid_to_hole[done_pid]
          os.remove(pid_file)
          continue
        else :
          continue
      # We don't need to to run the main loop to death
      # since we aren't using os.wait
      sleep(0.001)
      continue

    # If we make it here. We are running refine
    print('Open hole count: {0} {1}\n'.format(len(open_holes), time.time()))
    print('Best hole: {0}\n'.format(best_hole))
    if len(failed_holes) > 0:
      print('Deepest failed hole: {}\n'.format(sorted(failed_holes, key=len)[-1]))
      if len(open_holes) % 100 == 0:
        with open('deep_holes_' + name, 'w') as fp:
          num = min(len(failed_holes), 1000000)
          fp.write('\n'.join(sorted(failed_holes, key=len_and_lex)[:num]))
        with open('open_holes_' + name, 'w') as fp:
          num = min(len(open_holes), 1000000)
          fp.write('\n'.join(sorted(open_holes, key=len_and_lex)[:num]))
    else:
      print('Deepest failed hole: None\n')

    out = dest_dir + '/' + best_hole + '.out'
    err = dest_dir + '/' + best_hole + '.err'

    if best_hole == 'root':
      pid_word_search_depth = '-1'
    else: 
      pid_word_search_depth = word_search_depth

    # subdivide more at the beginning
    if len(open_holes) < child_limit - child_count and len(best_hole) < 40:
      invent_depth_local = str(max(int(invent_depth) // 2, 4))
    else:
      invent_depth_local = invent_depth

    treecat_command = '{0} {1} {2}'.format(treecat, src_dir, best_hole)
    refine_command = refine + \
        fill_holes + \
        improve_tree + \
        ' --box ' + best_hole + \
        ' --max_depth ' + max_depth + \
        ' --truncate_depth ' + truncate_depth + \
        ' --invent_depth ' + invent_depth_local + \
        ' --max_size ' + max_size + \
        ' --words ' + words_file + \
        ' --word_search_depth ' + pid_word_search_depth + \
        ' --impossible ' + impossible_file + \
        ' --bad_relators ' + bad_relator_file + \
        ' -m ' + cosh_mu_upper + \
        ' -r ' + sinh_tube_upper + \
        ' > ' + out  + ' 2> ' + err

    first_command = '{0} {1} {2} | head -1'.format(treecat, src_dir, best_hole)
    first = command_output(first_command).rstrip()

    print('|' + first + '|')

    if first[:1] == 'H' or len(first) == 0: # HOLE
      treecat_command = 'echo 1'

    command = treecat_command + ' | ' + refine_command
    print('Running with run count {1}: {0}\n'.format(
      command, refine_run_count))
    sys.stdout.flush()
    refine_run = Process(target=run_refine, args=(command, dest_dir,))
    refine_run.start()
    pid = refine_run.pid

    child_count += 1   
    refine_run_count += 1
    done.add(best_hole)
    active_pid_to_hole[pid] = best_hole
    done_failed = set()
    for h in failed_holes:
      if h.startswith(best_hole):
        done_failed.add(h)
    failed_holes.difference_update(done_failed)
