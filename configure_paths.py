#!/usr/bin/env python3
"""Explicitly adapt copied historical HPCC paths. Run only in a working copy."""
from pathlib import Path
import argparse
p=argparse.ArgumentParser();p.add_argument('--apply',action='store_true');a=p.parse_args()
root=Path(__file__).resolve().parent
old='/bigdata/koeniglab/oluwafea/project/data/cowpealandraces_490'
for f in (root/'manuscript_reanalysis').rglob('*'):
 if f.suffix not in {'.R','.py','.sh','.slurm'}:continue
 s=f.read_text();t=s.replace(old+'/go_window20kb_topGO',str(root/'manuscript_reanalysis/go_window20kb/topGO'))
 t=t.replace(old,str(root/'manuscript_reanalysis'))
 if t!=s:
  print(f.relative_to(root))
  if a.apply:f.write_text(t)
print('Applied' if a.apply else 'Preview only; use --apply in a working copy. Slurm modules and other source paths still need review.')
