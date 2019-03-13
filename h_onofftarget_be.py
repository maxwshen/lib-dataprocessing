# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, math, pickle, imp
sys.path.append('/home/unix/maxwshen/')
import fnmatch
import numpy as np
from collections import defaultdict
from mylib import util
from mylib import compbio
import pandas as pd

# Default params
inp_dir = _config.OUT_PLACE + 'g_poswise_be/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
lib_design = pd.read_csv(_config.DATA_DIR + 'library_12kChar.csv')


##
# Functions
##
def get_match_count(row):
  expected_guide_start_12kchar = 22

  grna = row['gRNA (20nt)']
  seq = row['Sequence context (56nt)']
  seq_g = seq[expected_guide_start_12kchar : expected_guide_start_12kchar + 20]

  mm = 0
  for nt1, nt2 in zip(grna, seq_g):
    if nt1 == nt2:
      mm += 1
  return mm

def on_off_target(nm):
  df = pd.read_csv(inp_dir + '%s.csv' % (nm), index_col = 0)

  df['Edited ULMI fraction'] = 1 - (df['Unedited ULMI'] / df['Total ULMI'])


  dd = defaultdict(list)
  timer = util.Timer(total = len(df))
  for idx, row in df.iterrows():
    name = row['Name (unique)']
    design_cat = row['Design category']

    match_count = get_match_count(row)

    if design_cat == 'guideseq':
      category = 'Off-target series'
      subcategory = name.split('_')[2] # gene name
    elif design_cat == 'mismatch':
      category = 'Mismatch series'
      subcategory = name.split('_')[1] # series number
    elif design_cat == 'chipseq':
      category = 'Chip series'
    elif design_cat == 'vivo':
      category = 'vivo'
      subcategory = 'vivo'
    else:
      assert match_count == 20, 'fail'
      category = 'On-target'
      subcategory = 'On-target'

    dd['Match count'].append(int(match_count))
    dd['Category'].append(category)
    dd['Subcategory'].append(subcategory)

    timer.update()

  for col in dd:
    df[col] = dd[col]

  df.to_csv(out_dir + '%s.csv' % (nm))

  return


##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
  print('Generating qsub scripts...')
  qsubs_dir = _config.QSUBS_DIR + NAME + '/'
  util.ensure_dir_exists(qsubs_dir)
  qsub_commands = []

  num_scripts = 0
  for bc in exp_design['Name']:
    if '12kChar' not in bc:
      continue
    if 'Cas9' in bc:
      continue

    command = 'python %s.py %s' % (NAME, bc)
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + 'q_%s_%s.sh' % (script_id, bc)
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -V -l h_rt=4:00:00 -wd %s %s &' % (_config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))

  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return


@util.time_dec
def main(nm = ''):
  print(NAME)

  # Function calls

  for bc in exp_design['Name']:
    if '12kChar' not in bc:
      continue
    if 'Cas9' in bc:
      continue
    if 'UT' in bc:
      continue
  
    print(bc)
    on_off_target(bc)

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1])
  else:
    gen_qsubs()