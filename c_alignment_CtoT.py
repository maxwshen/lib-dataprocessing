# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, copy
import numpy as np
from collections import defaultdict
sys.path.append('/home/unix/maxwshen/')
from mylib import util
from mylib import compbio
import pickle
import pandas as pd

# Default params
inp_dir = _config.OUT_PLACE + 'a_split/'
NAME = util.get_fn(__file__)
out_place = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_place)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')

lib_design = pd.read_csv(_config.DATA_DIR + 'library_CtoT.csv')

names_targets = dict()
for idx, row in lib_design.iterrows():
  names_targets[row['Name (unique)']] = row['Sequence context (56nt)']

all_names = list(lib_design['Name (unique)'])

##
# Alignments
##
def alignment(read, cand_idxs):
  seq_align_tool = '/home/unix/maxwshen/tools/seq-align/bin/needleman_wunsch'
  targets = [names_targets[nm] for nm in cand_idxs]
  aligns = []
  for target_seq in targets:
    try:
      targetse = 'GATGGGTGCGACGCGTCAT' + target_seq
      # targetse = 'GTCTGTGTTCCGTTGTCCGTGCTGTAACGAAAGGTGCAGTNNNNNNNNNNNNNNNGATGGGTGCGACGCGTCAT' + target_seq
      align = subprocess.check_output(seq_align_tool + ' --match 1 --mismatch -1 --gapopen -5 --gapextend -0 --freestartgap ' + read + ' ' + targetse, shell = True)
      aligns.append(align)
    except:
      pass

  if len(aligns) == 0:
    return None, None

  if len(aligns) > 1:
    best_align = pick_best_alignment(aligns)
    best_idx = cand_idxs[aligns.index(best_align)]
  else:
    best_align = aligns[0]
    best_idx = cand_idxs[0]
  if type(best_align) == bytes:
    best_align = best_align.decode('utf-8')
  # assert best_align[-2:] == '\n\n'
  best_align = best_align[:-2]
  return best_idx, best_align

def pick_best_alignment(aligns):
  scores = []
  for align in aligns:
    w = align.split()
    s1, s2 = w[0], w[1]
    score = 0
    for i in range(len(s1)):
      if s1[i] == s2[i] and s1[i] != '-':
        score += 1
    scores.append(score)
  best_idx = scores.index(max(scores))
  return aligns[best_idx]

##
# Locality sensitive hashing
##
def build_targets_better_lsh():
  lsh_dict = defaultdict(list)
  for nm in names_targets:
    target = names_targets[nm]
    kmers = get_lsh_kmers(target)
    for kmer in kmers:
      lsh_dict[kmer].append(nm)
  return lsh_dict

def get_lsh_kmers(target):
  kmer_len = 7
  kmers = []
  for idx in range(len(target) - kmer_len):
    kmer = target[idx : idx + kmer_len]
    kmers.append(kmer)
  return kmers

def find_best_designed_target(read, lsh_dict):
  kmers = get_lsh_kmers(read)
  scores = dict()
  for kmer in kmers:
    for exp in lsh_dict[kmer]:
      if exp not in scores:
        scores[exp] = 0
      scores[exp] += 1

  if len(scores) == 0:
    return []

  sorted_scores = sorted(scores, key = scores.get, reverse = True)
  best_score = scores[sorted_scores[0]]
  # cand_idxs = []
  # for exp in sorted_scores:
  #   if scores[exp] + 5 < best_score:
  #     break
  #   cand_idxs.append(exp)
  cand_idxs = [sorted_scores[0]]
  return cand_idxs

##
# Match gRNA to target
##
def compare_target_to_grna(cand_idxs, grna_read):
  # Only match best target site
  # If too many false negatives too few true positives), relax this
  # Slightly prefer false positives, which are cleaned up in polishing and genotyping

  idx = all_names.index(cand_idxs[0])
  expected_grna = lib_design['gRNA (20nt)'][idx]

  def get_match_score(query, ref):
    return sum([1 for idx in range(len(query)) if query[idx] != ref[idx]])

  scores = []
  for start_idx in range(18, 22 + 1):
    obs_grna = grna_read[start_idx : start_idx + 20]
    match_score = get_match_score(obs_grna, expected_grna)
    scores.append(match_score)

  if min(scores) > 2:
    return 'bad match'
  else:
    return 'ok'

def get_grna_from_nm(nm):
  idx = all_names.index(nm)
  return lib_design['gRNA (20nt)'][idx]

##
# IO
##
def store_alignment(alignment_buffer, idx, align_header, align, read_q):
  # Place quality scores in same configuration as read alignment
  read_align = align.split()[0]
  rq = ''
  curr_idx = 0
  for ch in read_align:
    if ch == '-':
      rq += ch
    else:
      rq += read_q[curr_idx]
      curr_idx += 1

  align_string = '%s\n%s\n%s\n' % (align_header, align, rq)
  alignment_buffer[idx].append(align_string)
  return

def init_alignment_buffer():
  alignment_buffer = defaultdict(list)
  return alignment_buffer

def flush_alignments(alignment_buffer, out_dir):
  print('Flushing... \n%s' % (datetime.datetime.now()))
  for exp in alignment_buffer:
    with open(out_dir + '%s.txt' % (exp), 'a') as f:
      for align in alignment_buffer[exp]:
        f.write(align)
  alignment_buffer = init_alignment_buffer()
  print('Done flushing.\n%s' % (datetime.datetime.now()))
  return

def prepare_outfns(out_dir):
  for exp in names_targets:
    out_fn = out_dir + '%s.txt' % (exp)
    util.exists_empty_fn(out_fn)
  return

def find_ulmi(line1):
  constant = 'ATGACGCGTCGCACCCATC'

  def get_match_score(query, ref):
    return sum([1 for idx in range(len(query)) if query[idx] != ref[idx]])

  ulmi_idx = 56 + len(constant)
  best_ulmi = line1[ulmi_idx : ulmi_idx + 15]
  for start_pos in range(56 - 5, 56 + 5):
    query = line1[start_pos : start_pos + 19]
    if get_match_score(query, constant) <= 3:
      best_ulmi = line1[ulmi_idx : ulmi_idx + 15]
      ulmi_idx = start_pos + len(constant)
      break

  return compbio.reverse_complement(best_ulmi), ulmi_idx

##
# Main
##
def matchmaker(nm, split):
  print(nm, split)
  stdout_fn = _config.SRC_DIR + 'nh_c_%s_%s.out' % (nm, split)
  util.exists_empty_fn(stdout_fn)
  out_dir = out_place + nm + '/' + split + '/'
  util.ensure_dir_exists(out_dir)

  read1_fn = inp_dir + '%s_R1_%s.fq' % (nm, split)
  read2_fn = inp_dir + '%s_R2_%s.fq' % (nm, split)
  
  lsh_dict = build_targets_better_lsh()
  alignment_buffer = init_alignment_buffer()

  prepare_outfns(out_dir)

  num_bad_matches = 0
  quality_pass = 0

  tot_lines = util.line_count(read1_fn)
  timer = util.Timer(total = tot_lines)
  with open(read1_fn) as f1, open(read2_fn) as f2:
    for i, (line1, line2) in enumerate(zip(f1, f2)):
      if i % 4 == 0:
        h1 = line1.strip()
        h2 = line2.strip()
      if i % 4 == 1:
        # RC of l1 contains target
        line1 = line1.strip()
        target_read = compbio.reverse_complement(line1[:56])
        ulmi, ulmi_idx = find_ulmi(line1)

        # l2 contains gRNA
        grna_read = line2.strip()

      if i % 4 == 3:

        q1, q2 = line1.strip(), line2.strip()
        read_q = q1[:56][::-1]
        ulmi_q = q1[ulmi_idx : ulmi_idx + len(ulmi)][::-1]
        grna_q = q2[18 : 22 + 20]

        qs = [ord(s)-33 for s in read_q + ulmi_q + grna_q]
        if np.mean(qs) >= 28:
          quality_pass += 1

          align_header = '>1_%s_%s' % (ulmi, ulmi_q)

          # Try to find designed target from LSH
          cand_idxs = find_best_designed_target(target_read, lsh_dict)
          if len(cand_idxs) > 0:

            bad_match = compare_target_to_grna(cand_idxs, grna_read)
            if bad_match == 'ok':
              # Run alignment and store in buffer
              best_idx, align = alignment(target_read, cand_idxs)
              if align is None:
                continue
              store_alignment(alignment_buffer, best_idx, align_header, align, read_q)
            else:
              num_bad_matches += 1
          else:
            num_bad_matches += 1


      if i % int(tot_lines / 100) == 1 and i > 1:
        # Flush alignment buffer
        flush_alignments(alignment_buffer, out_dir)
        alignment_buffer = init_alignment_buffer()

        # Stats for the curious
        with open(stdout_fn, 'a') as outf:
          outf.write('Time: %s\n' % (datetime.datetime.now()))
          outf.write('Progress: %s\n' % (i / int(tot_lines / 100)) )
          outf.write('Num. mismatched gRNA/target pairs: %s\n' % (num_bad_matches))
          outf.write('Frac. mismatched gRNA/target pairs: %s\n' % (num_bad_matches / quality_pass))

      timer.update()
  
  # Final flush
  flush_alignments(alignment_buffer, out_dir)

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
  for _nm in exp_design['Name']:
    if 'CtoT' not in _nm:
      continue
    for _split in range(60):
      command = 'python %s.py %s %s' % (NAME, _nm, _split)
      script_id = NAME.split('_')[0]

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, _nm, _split)
      with open(sh_fn, 'w') as f:
        f.write('#!/bin/bash\n%s\n' % (command))
      num_scripts += 1

      # Write qsub commands
      qsub_commands.append('qsub -V -l h_rt=14:00:00,h_vmem=2G -wd %s %s &' % (_config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands_CtoT.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))

  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return

@util.time_dec
def main(nm = '', split = ''):
  print(NAME)

  # Function calls
  matchmaker(nm, split) 
  return


if __name__ == '__main__':
  if len(sys.argv) > 2:
    main(nm = sys.argv[1], split = sys.argv[2])
  else:
    gen_qsubs()