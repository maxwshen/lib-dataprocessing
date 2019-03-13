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
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')

##
# Functions
##
def poswise_baseediting(nm, inp_place, lib_design, edit_type, lib_nm):
  dd = defaultdict(list)

  all_positions = ['pos%s' % (s) for s in range(-9, 29 + 1)]

  timer = util.Timer(total = len(lib_design))
  for idx, row in lib_design.iterrows():
    exp = row['Name (unique)']
    inp_dir = inp_place + exp + '/'
    inp_fn = inp_dir + 'wildtype.txt'
    timer.update()
    
    if os.path.exists(inp_fn):
      designed_seq = row['Sequence context (56nt)']
      d, total, unedited_ct = process_wt(inp_fn, designed_seq, edit_type, lib_nm)

      dd['Total ULMI'].append(total)
      dd['Unedited ULMI'].append(unedited_ct)
      edited_ct = total - unedited_ct
      dd['Edited ULMI'].append(edited_ct)
      dd['Edited fraction'].append(edited_ct / total)

      for pos_nm in d:
        dd[pos_nm].append(d[pos_nm])

    else:
      dd['Total ULMI'].append(np.nan)
      dd['Unedited ULMI'].append(np.nan)
      dd['Edited ULMI'].append(np.nan)
      dd['Edited fraction'].append(np.nan)

      for pos_nm in all_positions:
        dd[pos_nm].append(np.nan)

  for col in dd:
    lib_design[col] = dd[col]

  lib_design.to_csv(out_dir + nm + '.csv')

  return

def parse_header(header):
  count = int(header.split('_')[0].replace('>', ''))
  ulmi_count = int(header.split('_')[2])
  return count, ulmi_count

def process_wt(inp_fn, designed_seq, edit_type, lib_nm):
  if edit_type == 'CtoT':
    target_base = 'T'
    ref_base = 'C'
  elif edit_type == 'AtoG':
    target_base = 'G'
    ref_base = 'A'


  prefix_len = len('GATGGGTGCGACGCGTCAT')
  # expected_cutsite = prefix_len + 28

  if lib_nm == '12kChar':
    expected_cutsite_from_zero = 39
  elif lib_nm in ['AtoG', 'CtoT']:
    expected_cutsite_from_zero = 28

  offset = expected_cutsite_from_zero - 18

  start_pos = -9 + offset
  end_pos = 29 + offset + 1

  total = 0
  unedited_ct = 0
  d = dict()
  for idx in range(start_pos, end_pos):
    pos_nm = 'pos%s' % (idx - offset)
    if designed_seq[idx] == ref_base:
      d[pos_nm] = 0
    else:
      d[pos_nm] = np.nan

  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        count, ulmi_count = parse_header(line.strip())
        total += ulmi_count
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 2:
        ref = line.strip()

        num_edits = 0
        for idx in range(start_pos, end_pos):
          pos_nm = 'pos%s' % (idx - offset)
          if read[prefix_len + idx] == target_base and ref[prefix_len + idx] == ref_base:
            d[pos_nm] += ulmi_count
            num_edits += 1

        if num_edits == 0:
          unedited_ct += ulmi_count

  for key in d:
    d[key] /= total

  return d, total, unedited_ct

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
    if 'Cas9' in bc or 'UT' in bc:
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


def determine_base_editor_type(nm):
  edit_type = None

  ct_editors = [
    'AID',
    'BE4',
    'BE4-CP1028',
    'CDA',
    'evoAPOBEC',
    'eA3a',
    'eA3A',
  ]
  for ct_editor in ct_editors:
    if ct_editor in nm:
      edit_type = 'CtoT'

  ag_editors = [
    'ABE',
    'ABE-CP1040',
  ]
  for ag_editor in ag_editors:
    if ag_editor in nm:
      edit_type = 'AtoG'

  return edit_type

@util.time_dec
def main(nm = ''):
  print(NAME)

  # Function calls
  if '12kChar' in nm:
    lib_nm = '12kChar'
  elif 'AtoG' in nm:
    lib_nm = 'AtoG'
  elif 'CtoT' in nm:
    lib_nm = 'CtoT'

  edit_type = determine_base_editor_type(nm)

  inp_dir = _config.OUT_PLACE + 'c6_polish_%s/%s/' % (lib_nm, nm)
  lib_design = pd.read_csv(_config.DATA_DIR + 'library_%s.csv' % (lib_nm), index_col = 0)
  poswise_baseediting(nm, inp_dir, lib_design, edit_type, lib_nm)

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1])
  else:
    gen_qsubs()