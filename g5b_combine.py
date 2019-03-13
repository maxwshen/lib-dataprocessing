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
inp_dir = _config.OUT_PLACE + 'g5_combin_be/'
NAME = util.get_fn(__file__)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')


@util.time_dec
def main():
  print(NAME)

  for nm in exp_design['Name']:
    if 'Cas9' in nm:
      continue

    mdf = pd.DataFrame()
    data = None
    for start_idx in range(0, 12000, 2000):
      data_fn = inp_dir + '%s_%s_%s.pkl' % (nm, start_idx, start_idx + 1999)
      with open(data_fn, 'rb') as f:
        temp_d = pickle.load(f)
        if data is None:
          data = temp_d
        else:
          for key in temp_d:
            data[key] = temp_d[key]

    # Data
    with open(inp_dir + '%s.pkl' % (nm), 'wb') as f:
      pickle.dump(data, f)
    print('Wrote to %s.pkl' % (nm))


  print('Done')
  return


if __name__ == '__main__':
  main()