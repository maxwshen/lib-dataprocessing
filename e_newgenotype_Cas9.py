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
inp_dir = _config.OUT_PLACE + 'c6_polish_12kChar/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
lib_design = pd.read_csv(_config.DATA_DIR + 'library_12kChar.csv')

crispr_cutsite = None

##
# Data manipulation
##
def parse_header(header):
  # count = int(header.replace('>', '').replace('_', ''))
  count = int(header.split('_')[0].replace('>', ''))
  ulmi_count = int(header.split('_')[2])
  return count, ulmi_count

def get_indel_length(seq):
  # Assume exactly one indel in sequence
  central_indel_len = trim_start_end_dashes(seq).count('-')
  if central_indel_len != 0:
    return central_indel_len 
  else:
    # 3' end gap deletion
    return seq[::-1].index(seq[::-1].replace('-', '')[:1])

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
# Pandas dataframe manipulation
##
def compress(df):
  # Compress a dataframe along all columns except 'Count', summing it.
  cols = list(df.columns)
  cols.remove('Count')
  cols.remove('ULMI count')
  if len(cols) == 0:
    new_df = pd.DataFrame({'Count': [df['Count'].sum()], 'ULMI count': [df['ULMI count'].sum()]})
  else:
    grouped = df.groupby(cols)
    g = grouped.aggregate(sum)
    new_df = g.reset_index()
  return new_df

##
# Alignment properties
##
def calc_deletion_start_position(read, genome, del_len):
  # assumes exactly one deletion in alignment
  idx = genome.index(genome.replace('-', '')[:1]) + crispr_cutsite - del_len
  read_start = read.index(read.replace('-', '')[:1])
  for jdx in range(read_start, len(read)):
    if read[jdx] == '-':
      break
  del_start = jdx - idx
  return del_start, jdx

def calc_insertion_start_position(read, genome, ins_len):
  # assumes exactly one insertion in alignment
  genome_start_idx = genome.index(genome.replace('-', '')[:1])
  counter = 0
  for idx in range(genome_start_idx, len(genome)):
    if counter == crispr_cutsite:
      break
    if genome[idx] in 'ACGTN':
      counter += 1
  for jdx in range(genome_start_idx, len(genome)):
    if genome[jdx] == '-':
      break
  ins_start = jdx
  ins_end = ins_start + ins_len
  if abs(idx - ins_start) < abs(idx - ins_end):
    ins_position = ins_start - idx
  else:
    ins_position = ins_end - idx
  return ins_position, jdx

def check_mismatches(read, genome, del_start, del_len):
  # Check that 3 bp on both sides of deletion match perfectly
  match_5side = match(read[del_start - 3 : del_start], genome[del_start - 3 : del_start])
  del_end = del_start + del_len
  match_3side = match(read[del_end : del_end + 3], genome[del_end : del_end + 3])
  if bool(match_5side and match_3side):
    # is ok
    return 'no'
  else:
    return 'yes'

def match(s1, s2):
  for c1, c2 in zip(s1, s2):
    if c1 != c2:
      return False
  return True

def has_mh(read, genome, del_len, gt_pos):
  # assumes single deletion that arises from mmej mechanism
  if gt_pos < 0 or gt_pos > del_len:
    return 'na'
  genome_start_idx = genome.index(genome.replace('-', '')[:1])
  cut_loci = genome_start_idx + crispr_cutsite
  read_base_5side = read[cut_loci - del_len + gt_pos - 1]
  genome_base_3side = genome[cut_loci + gt_pos - 1]
  try:
    read_base_3side = read[cut_loci + gt_pos]
  except:
    read_base_3side = '*'
  genome_base_5side = genome[cut_loci -del_len + gt_pos]
  if gt_pos not in [0, del_len]:
    if bool(read_base_5side == genome_base_3side) or bool(read_base_3side == genome_base_5side):
      return 'yes'
  if gt_pos == 0:
    if bool(read_base_3side == genome_base_5side):
      return 'yes'
  if gt_pos == del_len:
    if bool(read_base_5side == genome_base_3side):
      return 'yes'
  return 'no'

def check_ins_templated(read, genome, is_pos, ins_len):
  # if the inserted sequence and some of the neighboring sequence is present in the wildtype sequence context, it's templated.

  def find_all_instances(query, seq):
    idxs = []
    for i in range(len(seq)):
      if seq[i : i + len(query)] == query:
        idxs.append(i)
    return idxs

  imer = read[is_pos : is_pos + ins_len]
  designed_genome = genome.replace('-', '')
  rc_designed_genome = compbio.reverse_complement(designed_genome)
  if imer not in designed_genome and imer not in rc_designed_genome:
    return 0, 'na', '', ''

  # try extending 5' side
  for idx in range(is_pos - 1, -1, -1):
    new_imer = read[idx : is_pos + ins_len]
    if new_imer not in designed_genome and new_imer not in rc_designed_genome:
      break
    # Template cannot be only where we are
    if new_imer in designed_genome and new_imer not in rc_designed_genome:
      inst = find_all_instances(new_imer, designed_genome)
      if len(inst) == 1 and idx in inst:
        break
  fiveside = idx + 1

  # try extending 3' side
  for idx in range(is_pos + ins_len + 1, len(read)):
    new_imer = read[fiveside : idx]
    if new_imer not in designed_genome and new_imer not in rc_designed_genome:
      break
    # Template cannot be only where we are
    if new_imer in designed_genome and new_imer not in rc_designed_genome:
      inst = find_all_instances(new_imer, designed_genome)
      if len(inst) == 1 and fiveside in inst:
        break
  threeside = idx - 1

  fiveside_seq = read[fiveside : is_pos]
  threeside_seq = read[is_pos + ins_len : threeside]

  # If no neighboring sequence is included in template, it's not templated.
  if len(fiveside_seq) == 0 or len(threeside_seq) == 0:
    return 0, 'na', '', ''

  template = read[fiveside : threeside]

  # get p2 and mh2
  if template in genome[:is_pos] or template in compbio.reverse_complement(genome[:is_pos].replace('-', '')):
    p2 = fiveside_seq
    mh2 = threeside_seq
  else:
    p2 = threeside_seq
    mh2 = fiveside_seq

  # Get template orientation
  if template in designed_genome and template not in rc_designed_genome:
    template_orientation = '+'
  if template not in designed_genome and template in rc_designed_genome:
    template_orientation = '-'
  if template in designed_genome and template in rc_designed_genome:
    template_orientation = 'both'

  # a random 5mer occurs in 55 bp at 5% rate. To threshold at various false positive rates, defer decision and just return length of longest template match.
  return len(template), template_orientation, p2, mh2

def check_ins_fivehomo(read, genome, is_pos, ins_len):
  imer = read[is_pos : is_pos + ins_len]
  five_template = read[is_pos - 1]
  if len(set(imer)) == 1 and five_template in set(imer):
    return 'yes'
  return 'no'

def standardized_del_gtpos(genome, del_len, del_start):
  genome_start_idx = genome.index(genome.replace('-', '')[:1])
  cutsite = genome_start_idx + crispr_cutsite
  if del_start < 0 or del_start > del_len:
    return del_start
  left = genome[cutsite - del_len : cutsite]
  right = genome[cutsite : cutsite + del_len]
  start_idx = max(len(right) - len(left), 0)
  mhs = []
  mh = [start_idx]
  for idx in range(min(len(right), len(left))):
    if left[idx] == right[start_idx + idx]:
      mh.append(start_idx + idx + 1)
    else:
      mhs.append(mh)
      mh = [start_idx + idx +1]
  mhs.append(mh)
  for mh in mhs:
    if del_start in mh:
      return max(mh)

##
# File-by-file subordinates for deletions and insertions
##
def get_deletions_in_file(inp_fn):
  global crispr_cutsite
  counts = []
  ulmi_counts = []
  mhs, mismatches = [], []
  gt_pos = []
  del_lens = []
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        header = line.strip()
        count, ulmi_count = parse_header(header)
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 2:
        genome = line.strip()
        del_len = get_indel_length(read)
        del_start, ds_pos = calc_deletion_start_position(read, genome, del_len)

        counts.append(count)
        ulmi_counts.append(ulmi_count)
        mh_basis = has_mh(read, genome, del_len, del_start)
        mhs.append(mh_basis)
        has_mismatch = check_mismatches(read, genome, ds_pos, del_len)
        mismatches.append(has_mismatch)
        gt_pos.append(standardized_del_gtpos(genome, del_len, del_start))
        del_lens.append(del_len)
  return counts, ulmi_counts, mhs, mismatches, gt_pos, del_lens

def get_insertions_in_file(inp_fn):
  global crispr_cutsite
  counts = []
  ulmi_counts = []
  ins_bases = []
  mismatches = []
  template_lens, fivehomopolymers = [], []
  template_orientations, p2s, mh2s = [], [], []
  ins_lens = []
  gt_pos = []
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        header = line.strip()
        count, ulmi_count = parse_header(header)
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 2:
        genome = line.strip()
        ins_len = get_indel_length(genome)
        ins_start, is_pos = calc_insertion_start_position(read, genome, ins_len)

        counts.append(count)
        ulmi_counts.append(ulmi_count)
        ins_bases.append(read[is_pos : is_pos + ins_len])
        has_mismatch = check_mismatches(read, genome, is_pos, ins_len)
        template_len, template_orientation, p2, mh2 = check_ins_templated(read, genome, is_pos, ins_len)
        template_lens.append(template_len)
        template_orientations.append(template_orientation)
        p2s.append(p2)
        mh2s.append(mh2)
        fivehomopolymer = check_ins_fivehomo(read, genome, is_pos, ins_len)
        fivehomopolymers.append(fivehomopolymer)
        mismatches.append(has_mismatch)
        gt_pos.append(ins_start)
        ins_lens.append(ins_len)
  return counts, ulmi_counts, ins_bases, mismatches, gt_pos, ins_lens, template_lens, fivehomopolymers, template_orientations, p2s, mh2s

##
# Main Subordinates for each Category
##
# noise filtered
def get_homopolymer(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'homopolymer.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'homopolymer'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

def get_hasN(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'hasN.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'hasN'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

def get_pcr_recombination(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'pcr_recombination.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'pcr_recombination'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

def get_poormatches(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'poormatches.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'poormatches'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

def get_cutsite_not_sequenced(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'cutsite_not_sequenced.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'cutsite_not_sequenced'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

def get_read_too_short(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'read_too_short.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'read_too_short'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

# indels
def get_deletions(master_df, exp, exp_dir):
  # Iterate over all insertion files. Build a temporary dataframe to add to the master dataframe.
  counts = []
  ulmi_counts = []
  mhs = []
  mismatches = []
  gt_poses = []
  del_lengths = []
  by_crispr_cuts = []
  categories = []
  
  for fn in os.listdir(exp_dir):
    if fnmatch.fnmatch(fn, 'del*.txt'):
      if 'not' not in fn:
        category = 'del'
      elif 'notcrispr' in fn:
        category = 'del_notcrispr'
      else:
        category = 'del_notatcut'

      inp_fn = exp_dir + fn
      ans = get_deletions_in_file(inp_fn)
      count, ulmi_count, mh, mismatch, gt_pos, del_lens = ans
      counts += count
      ulmi_counts += ulmi_count
      mhs += mh
      mismatches += mismatch
      gt_poses += gt_pos
      del_lengths += del_lens
      categories += [category] * len(del_lens)

  temp_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts, 'Indel with Mismatches': mismatches, 'Microhomology-Based': mhs, 'Genotype Position': gt_poses, 'Length': del_lengths, 'Category': categories})
  new_df = compress(temp_df)
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

def get_insertions(master_df, exp, exp_dir):
  # Iterate over all insertion files. Build a temporary dataframe to add to the master dataframe.
  counts = []
  ulmi_counts = []
  ins_bases = []
  mismatches = []
  gt_poses = []
  ins_lengths = []
  by_crispr_cuts = []
  categories = []
  template_lens, fivehomopolymers = [], []
  template_orientations, p2s, mh2s = [], [], []
  
  for fn in os.listdir(exp_dir):
    if fnmatch.fnmatch(fn, 'ins*.txt'):
      if 'not' not in fn:
        category = 'ins'
      elif 'notcrispr' in fn:
        category = 'ins_notcrispr'
      else:
        category = 'ins_notatcut'

      inp_fn = exp_dir + fn
      ans = get_insertions_in_file(inp_fn)
      count, ulmi_count, ins_base, mismatch, gt_pos, ins_lens, template_len, fivehomopolymer, template_orientation, p2, mh2 = ans
      ins_bases += ins_base
      counts += count
      ulmi_counts += ulmi_count
      mismatches += mismatch
      template_lens += template_len
      fivehomopolymers += fivehomopolymer
      template_orientations += template_orientation
      p2s += p2
      mh2s += mh2
      gt_poses += gt_pos
      ins_lengths += ins_lens
      categories += [category] * len(ins_lens)

  temp_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts, 'Inserted Bases': ins_bases, 'Indel with Mismatches': mismatches, 'Genotype Position': gt_poses, 'Length': ins_lengths, 'Category': categories, 'Ins Template Length': template_lens, 'Ins Fivehomopolymer': fivehomopolymers, 'Ins Template Orientation': template_orientations, 'Ins p2': p2s, 'Ins mh2': mh2s})
  new_df = compress(temp_df)
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

def get_combination_indels(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'combination_indel.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  temp_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df = compress(temp_df)
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'combination_indel'
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

# secondary
def get_forgiven_indels(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'forgiven_indel.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'forgiven_indel'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df  

def get_forgiven_combination_indels(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'forgiven_combination_indel.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'forgiven_combination_indel'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

# other
def get_combination_indels_notcrispr(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'combination_indel_notcrispr.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  new_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Length'] = -1
  new_df['Category'] = 'combination_indel_notcrispr'
  new_df = compress(new_df)
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

def get_other(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'other.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)

  temp_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df = compress(temp_df)
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'other'
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

# wildtypes
def get_wildtype(master_df, exp, exp_dir):
  counts = []
  ulmi_counts = []
  inp_fn = exp_dir + 'wildtype.txt'
  if os.path.isfile(inp_fn):
    with open(inp_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          count, ulmi_count = parse_header(header)
          counts.append(count)
          ulmi_counts.append(ulmi_count)
  temp_df = pd.DataFrame({'Count': counts, 'ULMI count': ulmi_counts})
  new_df = compress(temp_df)
  new_df['_ExpDir'] = exp_dir
  new_df['_Experiment'] = exp
  new_df['Category'] = 'wildtype'
  master_df = pd.concat([master_df, new_df], ignore_index = True, sort = True)
  return master_df

##
# main 
##
def genotype_data(inp_dir, out_dir, nm, start, end):
  print(nm)
  start, end = int(start), int(end)
  master_df = pd.DataFrame()

  global crispr_cutsite
  crispr_cutsite = len('GATGGGTGCGACGCGTCAT') + 39

  nms = lib_design['Name (unique)'][start : end + 1]

  timer = util.Timer(total = len(nms))
  for exp in nms:
    exp_dir = '%s%s/' % (inp_dir, exp)
    if not os.path.isdir(exp_dir):
      return

    # Noise categories
    master_df = get_homopolymer(master_df, exp, exp_dir)
    master_df = get_hasN(master_df, exp, exp_dir)
    master_df = get_pcr_recombination(master_df, exp, exp_dir)
    master_df = get_poormatches(master_df, exp, exp_dir)
    master_df = get_cutsite_not_sequenced(master_df, exp, exp_dir)
    master_df = get_read_too_short(master_df, exp, exp_dir)

    # Primary categories
    master_df = get_deletions(master_df, exp, exp_dir)
    master_df = get_insertions(master_df, exp, exp_dir)
    master_df = get_combination_indels(master_df, exp, exp_dir)

    # Secondary categories
    master_df = get_forgiven_indels(master_df, exp, exp_dir)
    master_df = get_forgiven_combination_indels(master_df, exp, exp_dir)

    # Other categories
    master_df = get_combination_indels_notcrispr(master_df, exp, exp_dir)
    master_df = get_other(master_df, exp, exp_dir)

    # Wildtypes
    master_df = get_wildtype(master_df, exp, exp_dir)

    timer.update()

  seq_contexts = []
  for s in master_df['_Experiment']:
    crit = (lib_design['Name (unique)'] == s)
    seq = lib_design[crit]['Sequence context (56nt)'].iloc[0]
    seq_contexts.append(seq)
  master_df['_Sequence Context'] = seq_contexts

  master_df['_Cutsite'] = crispr_cutsite

  master_df.to_csv(out_dir + '%s_genotypes_%s.csv' % (nm, start))
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
    if 'Cas9' not in bc and 'UT' not in bc:
      continue

    for start in range(0, 12000, 1000):
      end = start + 999
      command = 'python %s.py %s %s %s' % (NAME, bc, start, end)
      script_id = NAME.split('_')[0]

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, bc, start)
      with open(sh_fn, 'w') as f:
        f.write('#!/bin/bash\n%s\n' % (command))
      num_scripts += 1

      # Write qsub commands
      qsub_commands.append('qsub -V -l h_rt=14:00:00 -wd %s %s &' % (_config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))

  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return

@util.time_dec
def main(nm = '', start = '', end = ''):
  print(NAME)

  # Function calls
  genotype_data('%s%s/' % (inp_dir, nm), out_dir, nm, start, end)
  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1], start = sys.argv[2], end = sys.argv[3])
  else:
    gen_qsubs()