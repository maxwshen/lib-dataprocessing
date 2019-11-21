# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, math, pickle
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
def get_combinatorial_baseedits(nm, inp_place, lib_design, lib_nm, edit_type, start_idx, end_idx):

  data = dict()
  minq = dict()

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

    inp_fns = []
    for fold in folds:
      inp_fn = inp_place + '%s/%s.txt' % (fold, exp)
      if os.path.exists(inp_fn):
        inp_fns.append(inp_fn)

    df, mq = process_aligns(inp_fns, designed_seq, lib_nm, edit_type)
    if df is not None:
      data[exp] = df

    if mq is not None:
      minq[exp] = mq

  import pickle
  pkl_fn = out_dir + '%s_%s_%s.pkl' % (nm, start_idx, end_idx)
  with open(pkl_fn, 'wb') as f:
    pickle.dump(data, f)

  minq_fn = out_dir + '%s_minq_%s_%s.pkl' % (nm, start_idx, end_idx)
  with open(minq_fn, 'wb') as f:
    pickle.dump(minq, f)

  return

def parse_header(header):
  count = int(header.split('_')[0].replace('>', ''))
  return count

def process_aligns(inp_fns, designed_seq, lib_nm, edit_type):
  prefix_len = len('GATGGGTGCGACGCGTCAT')

  if lib_nm == '12kChar':
    start_pos, end_pos = 0, 56
    target_site_len = 56
    grna_end_pos = 42
  elif lib_nm in ['AtoG', 'CtoT']:
    start_pos, end_pos = 0, 56
    target_site_len = 56
    grna_end_pos = 31
  elif lib_nm == 'PAMvar':
    start_pos, end_pos = 0, 61
    target_site_len = 61
    grna_end_pos = 33
  elif lib_nm == 'LibA':
    start_pos, end_pos = 0, 55
    target_site_len = 55
    grna_end_pos = 30

  pos0_idx = grna_end_pos - 21

  if edit_type == 'CtoT':
    target_bases = ['C', 'G']
  elif edit_type == 'AtoG':
    target_bases = ['A', 'C']
  elif edit_type == 'UT':
    target_bases = ['C', 'G', 'A']

  dd = defaultdict(list)
  mq = dict()

  for inp_fn in inp_fns:
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
          if '-' in ref:
            continue

          qs = [ord(s)-33 for s in line.strip()]

          for idx in range(grna_end_pos):          
            obs_nt = read[prefix_len + idx]
            ref_nt = ref[prefix_len + idx]
            q = qs[prefix_len + idx]

            if ref_nt not in target_bases:
              continue

            if obs_nt == '-' or ref_nt == '-' or q < 30:
              obs_nt = '.'
            if q < 30:
              obs_nt = '.'

            pos = idx - pos0_idx
            col_nm = '%s%s' % (ref_nt, pos)
            dd[col_nm].append(obs_nt)
            # if len(dd[col_nm]) == 1 and max([len(dd[s]) for s in dd]) != 1:
              # import code; code.interact(local=dict(globals(), **locals()))

            if q >= 30:
              if col_nm not in mq:
                mq[col_nm] = q
              else:
                mq[col_nm] = min(q, mq[col_nm])

          dd['Count'].append(1)

  if len(dd) == 0:
    return None, None

  df = pd.DataFrame(dd)
  nt_cols = [s for s in df.columns if s != 'Count']

  if len(nt_cols) == 0:
    return None, None

  df['Grouped count'] = df.groupby(nt_cols)['Count'].transform('sum')
  df = df.drop_duplicates(subset = nt_cols)
  df['Count'] = df['Grouped count']
  df = df.drop(columns = ['Grouped count'])
  df = df.sort_values(by = 'Count', ascending = False)
  df = df.reset_index()
  return df, mq

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


  if 'UT' in nm:
    edit_type = 'UT'
  else:
    edit_type = determine_base_editor_type(nm)

  inp_dir = _config.OUT_PLACE + 'c_alignment_%s/%s/' % (lib_nm, nm)
  lib_design = pd.read_csv(_config.DATA_DIR + 'library_%s.csv' % (lib_nm))

  get_combinatorial_baseedits(nm, inp_dir, lib_design, lib_nm, edit_type, int(start_idx), int(end_idx))

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1], start_idx = sys.argv[2], end_idx = sys.argv[3])
  else:
    gen_qsubs()