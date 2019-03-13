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
# Alignment manipulations
##
def detect_wildtype(read, genome):
  ts1, ts2 = trim_start_end_dashes(read), trim_start_end_dashes(genome)
  num_dels = len(ts1.replace('-', ' ').split()) - 1
  num_ins = len(ts2.replace('-', ' ').split()) - 1
  if num_dels == 0 and num_ins == 0:
    if '-' not in ts1 and '-' not in ts2:
      return True
  return False

def trim_start_end_dashes(seq):
  alphabet = set(list(seq))
  if '-' in alphabet:
    alphabet.remove('-')
  alphabet = list(alphabet)
  start = min([seq.index(s) for s in alphabet])
  end = min([seq[::-1].index(s) for s in alphabet])
  if end > 0:
    return seq[start:-end]
  else:
    return seq[start:]

##
# Iterator
##
def get_poswise_baseedits(nm, inp_place, lib_design, lib_nm, start_idx, end_idx):
  data = dict()
  stats = defaultdict(list)

  folds = os.listdir(inp_place)

  rows = lib_design.iloc[start_idx : end_idx + 1].iterrows()
  timer = util.Timer(total = end_idx + 1 - start_idx)
  for idx, row in rows:
    timer.update()
    exp = row['Name (unique)']

    if lib_nm != 'PAMvar':
      designed_seq = row['Sequence context (56nt)']
    else:
      designed_seq = row['Sequence context (61nt)']

    fold_dd = defaultdict(lambda: 0)
    for fold in folds:
      inp_fn = inp_place + '%s/%s.txt' % (fold, exp)
      
      if not os.path.exists(inp_fn):
        continue

      d, total, unedited_ct = process_aligns(inp_fn, designed_seq, lib_nm)

      if exp not in data:
        data[exp] = d
      else:
        data[exp] += d

      fold_dd['Total count'] += total
      fold_dd['Unedited count'] += unedited_ct
      fold_dd['Edited count'] += total - unedited_ct

    stats['Name (unique)'].append(exp)
    for col in fold_dd:
      stats[col].append(fold_dd[col])

  stats_df = pd.DataFrame(stats)
  stats_df.to_csv(out_dir + '%s_stats_%s_%s.csv' % (nm, start_idx, end_idx))

  import pickle
  pkl_fn = out_dir + '%s_%s_%s.pkl' % (nm, start_idx, end_idx)
  with open(pkl_fn, 'wb') as f:
    pickle.dump(data, f)

  return

def parse_header(header):
  count = int(header.split('_')[0].replace('>', ''))
  return count

def process_aligns(inp_fn, designed_seq, lib_nm):
  prefix_len = len('GATGGGTGCGACGCGTCAT')
  # expected_cutsite = prefix_len + 28

  if lib_nm == '12kChar':
    start_pos, end_pos = 0, 56
    target_site_len = 56
    # expected_cutsite_from_zero = 39
  elif lib_nm in ['AtoG', 'CtoT']:
    start_pos, end_pos = 0, 56
    target_site_len = 56
    # expected_cutsite_from_zero = 28
  elif lib_nm == 'PAMvar':
    start_pos, end_pos = 0, 61
    target_site_len = 61
    # expected_cutsite_from_zero = 30

  # offset = expected_cutsite_from_zero - 18

  # start_pos = -9 + offset
  # end_pos = 29 + offset + 1

  total = 0
  unedited_ct = 0

  nts = list('ACGT')
  mapper = {nts[s]: s for s in range(len(nts))}
  d = np.zeros( (target_site_len, len(nts)) )

  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        count = parse_header(line.strip())
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 2:
        ref = line.strip()

      if i % 4 == 3:
        if not detect_wildtype(read, ref):
          continue

        total += count

        num_edits = 0
        qs = [ord(s)-33 for s in line.strip()]
        for idx in range(start_pos, end_pos):
          # pos_nm = 'pos%s' % (idx)
          # pos_nm = 'pos%s' % (idx - offset)
          obs_nt = read[prefix_len + idx]
          ref_nt = ref[prefix_len + idx]
          q = qs[prefix_len + idx]

          if obs_nt == '-' or ref_nt == '-':
            # Should be very rare in wt
            continue

          if q < 30:
            continue

          d[idx][mapper[obs_nt]] += count
          if ref_nt != obs_nt:
            num_edits += 1

        if num_edits == 0:
          unedited_ct += count

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
    if 'Cas9' in bc:
      continue

    for start_idx in range(0, 12000, 2000):
      command = 'python %s.py %s %s %s' % (NAME, bc, start_idx, start_idx + 1999)
      script_id = NAME.split('_')[0]

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, bc, start_idx)
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
def main(nm = '', start_idx = '', end_idx = ''):
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

  inp_dir = _config.OUT_PLACE + 'c_alignment_%s/%s/' % (lib_nm, nm)
  lib_design = pd.read_csv(_config.DATA_DIR + 'library_%s.csv' % (lib_nm))
  get_poswise_baseedits(nm, inp_dir, lib_design, lib_nm, int(start_idx), int(end_idx))

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1], start_idx = sys.argv[2], end_idx = sys.argv[3])
  else:
    gen_qsubs()