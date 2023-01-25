#!/usr/bin/env python3

import os, subprocess, sys, getopt, glob, time, re, csv
from numpy import sinh, cosh, exp, isnan, isinf

test_data_file = 'margulis_testing/testing_data.csv' 
test_bin = 'test'

# A few useful functions
def command_output(command):
    try :
        # using subprocess.check_output returns empty for some reason.
        # using getoutput instead.
        return subprocess.getoutput(command)
    except:
        print('Error in {0}\n'.format(command)) 
        sys.exit(1)

def test_bounds(box, four_cosh_margulis, exp_2t, debug = False, debug_file = None) :
  print(box)
  out = command_output('./{} {}'.format(test_bin, box))
  (margulis_out, exp_2t_out) = out.splitlines()[-2:]
  margulis_match = re.fullmatch('.*margulis.*between\s(\S*)\s\((\S*)\)\sand\s(\S*)\s\((\S*)\)', margulis_out)
  exp_2t_match = re.fullmatch('.*2t.*between\s(\S*)\s\((\S*)\)\sand\s(\S*)\s\((\S*)\)', exp_2t_out)
  four_cosh_margulis_lb_lb = float(margulis_match.group(1))
  four_cosh_margulis_up_up = float(margulis_match.group(3))
  exp_2t_lb_lb = float(exp_2t_match.group(1))
  exp_2t_up_up = float(exp_2t_match.group(3))
  if debug :
    print('Asserting that {} < {} < {} and {} < {} < {}'.format(four_cosh_margulis_lb_lb, four_cosh_margulis, four_cosh_margulis_up_up,
                                                                exp_2t_lb_lb, exp_2t, exp_2t_up_up)) 
    if not debug_file : debug_file = box
    with open(debug_file, 'w') as fp :
      fp.write(out)

  assert not (isnan(four_cosh_margulis_lb_lb) or isnan(four_cosh_margulis_up_up) or 
              isnan(exp_2t_lb_lb) or isnan(exp_2t_up_up))
  assert not (isinf(four_cosh_margulis_lb_lb) or isinf(four_cosh_margulis_up_up) or 
              isinf(exp_2t_lb_lb) or isinf(exp_2t_up_up))
  assert four_cosh_margulis_lb_lb < four_cosh_margulis 
  assert four_cosh_margulis_up_up > four_cosh_margulis 
  assert exp_2t_lb_lb < exp_2t 
  assert exp_2t_up_up > exp_2t 

if __name__ == '__main__' :
  skipped = set()
  failed = set()
  with open(test_data_file) as test_data :
    test_data_reader = csv.reader(test_data)
    next(test_data_reader)
    for row in test_data_reader :
      name = row[0]
      sinh_ortho = sinh(eval(row[5]))
      # ortho is out of bounds 
      # TODO make ortho.real be the largest box coordinate
      if abs(sinh_ortho.real) > 5.65 :
        skipped.add(name)
        continue
      four_cosh_margulis = 4*cosh(float(row[2]))
      exp_2t = exp(2*float(row[6]))
      boxes = eval(row[11])
      for box in boxes :
        print("Working on " + name)
        try :
          test_bounds(box, four_cosh_margulis, exp_2t)
        except :
          print('Debug testing {} with 4cosh(margulis) = {} and exp(2t) = {} at box {}'.format(name, four_cosh_margulis, exp_2t, box))
          try :
            test_bounds(box, four_cosh_margulis, exp_2t, debug = True, debug_file = name + '_debug')
          except :
            failed.add(name)
  print('Failed for manifolds: {}'.format(sorted(failed)))
