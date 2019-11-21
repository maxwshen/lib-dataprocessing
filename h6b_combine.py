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
inp_dir = _config.OUT_PLACE + 'h6_anyindel/'
NAME = util.get_fn(__file__)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')


@util.time_dec
def main():
  print(NAME)

  combined_d = dict()

  for nm in exp_design['Name']:
    print(nm)
    mdf = pd.DataFrame()
    for start_idx in range(0, 12000, 2000):
      d = pickle.load(open(inp_dir + '%s_%s_%s.pkl' % (nm, start_idx, start_idx + 1999), 'rb'))
      for key in d:
        combined_d[key] = d[key]

    # Data
    with open(inp_dir + '%s.pkl' % (nm), 'wb') as f:
      pickle.dump(combined_d, f)

  print('Done')
  return


if __name__ == '__main__':
  main()