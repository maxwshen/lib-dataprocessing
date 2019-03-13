# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, math, pickle, imp
sys.path.append('/home/unix/maxwshen/')
import fnmatch
import numpy as np
from collections import defaultdict
from mylib import util
import pandas as pd

# Default params
inp_dir = _config.OUT_PLACE + 'e_newgenotype_Cas9/'
NAME = util.get_fn(__file__)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')

@util.time_dec
def main():
  print(NAME)
  dtypes = {'Category': str, 'Count': float, 'Genotype Position': float, 'Indel with Mismatches': str, 'Ins Fivehomopolymer': str, 'Ins Template Length': float, 'Ins mh2': str, 'Ins p2': str, 'Inserted Bases': str, 'Length': float, 'Microhomology-Based': str, '_ExpDir': str, '_Experiment': str, '_Sequence Context': str, '_Cutsite': int}

  for nm in exp_design['Name']:
    if 'Cas9' not in nm and 'UT' not in nm:
      continue
    if '12kChar' not in nm:
      continue

    mdf = pd.DataFrame()
    for start in range(0, 12000, 1000):
      fn = inp_dir + '%s_genotypes_%s.csv' % (nm, start)
      df = pd.read_csv(fn, index_col = 0, dtype = dtypes)
      mdf = pd.concat([mdf, df], ignore_index = True)

    mdf.to_csv(inp_dir + '%s.csv' % (nm))
    print('Wrote to %s.csv' % (inp_dir + '%s' % (nm)))
  print('Done')
  return


if __name__ == '__main__':
  main()