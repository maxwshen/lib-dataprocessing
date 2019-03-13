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
def get_all_mutations(nm, inp_place, lib_design, lib_nm):
  dd = defaultdict(list)

  timer = util.Timer(total = len(lib_design))
  for idx, row in lib_design.iterrows():
    exp = row['Name (unique)']
    inp_dir = inp_place + exp + '/'
    inp_fn = inp_dir + 'wildtype.txt'
    timer.update()
    
    if not os.path.exists(inp_fn):
      continue

    if lib_nm != 'PAMvar':
      designed_seq = row['Sequence context (56nt)']
    else:
      designed_seq = row['Sequence context (61nt)']

    d, total, unedited_ct = process_wt(inp_fn, designed_seq, lib_nm)

    for pos_nm in d:
      for ref_nt in d[pos_nm]:
        ref_nt_count = sum(d[pos_nm][ref_nt].values())
        for obs_nt in d[pos_nm][ref_nt]:
          ct = d[pos_nm][ref_nt][obs_nt]

          dd['Position'].append(pos_nm)
          dd['Ref nt'].append(ref_nt)
          dd['Obs nt'].append(obs_nt)
          dd['Count'].append(ct)
          dd['Ref nt count'].append(ref_nt_count)

          dd['Name (unique)'].append(exp)
          dd['Design category'].append(row['Design category'])
          dd['Total count'].append(total)
          edited_ct = total - unedited_ct
          dd['Edited count'].append(edited_ct)
          # dd['Edited fraction'].append(edited_ct / total)

  df = pd.DataFrame(dd)
  df.to_csv(out_dir + nm + '.csv')

  return

def parse_header(header):
  count = int(header.split('_')[0].replace('>', ''))
  ulmi_count = int(header.split('_')[2])
  return count, ulmi_count

def process_wt(inp_fn, designed_seq, lib_nm):
  prefix_len = len('GATGGGTGCGACGCGTCAT')
  # expected_cutsite = prefix_len + 28

  if lib_nm == '12kChar':
    start_pos, end_pos = 0, 56
    # expected_cutsite_from_zero = 39
  elif lib_nm in ['AtoG', 'CtoT']:
    start_pos, end_pos = 0, 56
    # expected_cutsite_from_zero = 28
  elif lib_nm == 'PAMvar':
    start_pos, end_pos = 0, 61
    # expected_cutsite_from_zero = 30

  # offset = expected_cutsite_from_zero - 18

  # start_pos = -9 + offset
  # end_pos = 29 + offset + 1

  total = 0
  unedited_ct = 0
  d = dict()
  reference = designed_seq[start_pos : end_pos]
  for idx, ref_nt in zip(range(start_pos, end_pos), reference):
    pos_nm = 'pos%s' % (idx)
    # pos_nm = 'pos%s' % (idx - offset)
    d[pos_nm] = dict()
    d[pos_nm][ref_nt] = dict()
    for nt2 in list('ACGT'):
      d[pos_nm][ref_nt][nt2] = 0

  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        count, ulmi_count = parse_header(line.strip())
        total += count
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 2:
        ref = line.strip()

      if i % 4 == 3:
        qs = [ord(s)-33 for s in line.strip()]

        num_edits = 0
        for idx in range(start_pos, end_pos):
          pos_nm = 'pos%s' % (idx)
          # pos_nm = 'pos%s' % (idx - offset)
          obs_nt = read[prefix_len + idx]
          ref_nt = ref[prefix_len + idx]
          q = qs[prefix_len + idx]

          if obs_nt == '-' or ref_nt == '-':
            # Should be very rare in wt
            continue

          if q < 30:
            continue

          d[pos_nm][ref_nt][obs_nt] += count
          if ref_nt != obs_nt:
            num_edits += 1

        if num_edits == 0:
          unedited_ct += count

  # Note: ref_nt is redundant since ref is always the same
  # for pos_nm in d:
  #   for ref_nt in d[pos_nm]:
  #     total = sum(d[pos_nm][ref_nt].values())

  #     for obs_nt in d[pos_nm][ref_nt]:
  #       try:
  #         d[pos_nm][ref_nt][obs_nt] /= total
  #       except ZeroDivisionError:
  #         pass

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
    if 'UT' not in bc:
      continue

    command = 'python %s.py %s' % (NAME, bc)
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + 'q_%s_%s.sh' % (script_id, bc)
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -V -l h_rt=2:00:00,h_vmem=1G -wd %s %s &' % (_config.SRC_DIR, sh_fn))

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
  if '12kChar' in nm:
    lib_nm = '12kChar'
  elif 'AtoG' in nm:
    lib_nm = 'AtoG'
  elif 'CtoT' in nm:
    lib_nm = 'CtoT'
  elif 'PAMvar' in nm:
    lib_nm = 'PAMvar'


  inp_dir = _config.OUT_PLACE + 'c6_polish_%s/%s/' % (lib_nm, nm)
  lib_design = pd.read_csv(_config.DATA_DIR + 'library_%s.csv' % (lib_nm))
  get_all_mutations(nm, inp_dir, lib_design, lib_nm)

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1])
  else:
    gen_qsubs()