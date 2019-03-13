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
def count_reads(exp, inp_dir, lib_design):

  dd = defaultdict(list)
  timer = util.Timer(total = len(lib_design['Name (unique)']))
  for nm in lib_design['Name (unique)']:
    ctd = get_counts_subfold(inp_dir + nm + '/')

    dd['Name'].append(nm)

    try:
      dd['Total count'].append(ctd['Total count'])
    except:
      import code; code.interact(local=dict(globals(), **locals()))
    dd['Total ULMI count'].append(ctd['Total ULMI count'])

    dd['WT count'].append(ctd['WT count'])
    dd['WT ULMI count'].append(ctd['WT ULMI count'])

    dd['Indel count'].append(ctd['Indel count'])
    dd['Indel ULMI count'].append(ctd['Indel ULMI count'])

    timer.update()

  df = pd.DataFrame(dd)
  df.to_csv(out_dir + '%s.csv' % (exp))
  return

def get_counts_subfold(inp_dir):
  if not os.path.isdir(inp_dir):
    return None, None

  ctd = defaultdict(lambda: 0)
  for fn in os.listdir(inp_dir):
    count, ulmi_count = get_counts_file(inp_dir + fn)
    ctd['Total count'] += count
    ctd['Total ULMI count'] += ulmi_count
    if 'wildtype' in fn:
      ctd['WT count'] += count
      ctd['WT ULMI count'] += ulmi_count
    if 'del' in fn or 'ins' in fn:
      ctd['Indel count'] += count
      ctd['Indel ULMI count'] += ulmi_count
  return ctd

def get_counts_file(inp_fn):
  counts, ulmi_counts = 0, 0
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        header = line.strip()
        count = int(header.split('_')[0].replace('>', ''))
        ulmi_count = int(header.split('_')[2])

        counts += count
        ulmi_counts += ulmi_count

  return counts, ulmi_counts

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
  for lib_type in ['12kChar', 'AtoG', 'CtoT', 'PAMvar']:
    if lib_type in nm:
      suffix = lib_type

  inp_dir = _config.OUT_PLACE + 'c6_polish_%s/%s/' % (suffix, nm)
  lib_design = pd.read_csv(_config.DATA_DIR + 'library_%s.csv' % (suffix))
  count_reads(nm, inp_dir, lib_design)

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1])
  else:
    gen_qsubs()