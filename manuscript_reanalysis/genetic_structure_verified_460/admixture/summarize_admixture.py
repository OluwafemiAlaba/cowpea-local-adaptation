"""Validate replicated ADMIXTURE outputs and align ancestry labels within each K.

Run after copying complete run logs, Q matrices and the input FAM from HPCC.
Cluster numbers are arbitrary; minimum-cost alignment is within K only.
Replicate variation is an optimization diagnostic, not biological uncertainty.
"""
from pathlib import Path
import csv
import json
import re
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'results'
OUT = ROOT / 'summary'
OUT.mkdir(exist_ok=True)
ids = [line.split()[1] for line in (RAW/'input/primary460.fam').read_text().splitlines()]
assert len(ids) == len(set(ids)) == 460
with (ROOT.parents[1]/'analysis_panel_460_ids.csv').open() as f:
    approved = [row['Accession'] for row in csv.DictReader(f)]
assert set(ids) == set(approved)

def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def align(q, reference):
    cost = ((q[:, :, None] - reference[:, None, :])**2).sum(axis=0)
    # Exact subset dynamic programming for the small K<=10 assignment problem.
    # State: original columns already used; assign the next reference column.
    states = {0: (0.0, ())}
    for target in range(q.shape[1]):
        next_states = {}
        for mask, (value, perm) in states.items():
            for original in range(q.shape[1]):
                if mask & (1 << original):
                    continue
                newmask = mask | (1 << original)
                candidate = (value + cost[original, target], perm + (original,))
                if newmask not in next_states or candidate < next_states[newmask]:
                    next_states[newmask] = candidate
        states = next_states
    permutation = np.array(states[(1 << q.shape[1])-1][1])
    return q[:, permutation], permutation

def assignment(q):
    return np.where(q.max(axis=1) >= .70, q.argmax(axis=1)+1, 0)

run_rows, k_rows, checksums = [], [], {}
display_label_alignment = {}
for k in range(1, 11):
    runs = []
    for rep in range(1, 11):
        folder = RAW/f'K{k}'/f'rep{rep}'
        assert (folder/'COMPLETE').read_text().strip() == 'COMPLETE'
        log = (folder/'admixture.log').read_text()
        assert 'Converged in ' in log and 'Writing output files.' in log
        assert 'Size of G: 460x4194' in log
        parameters = (folder/'run_parameters.txt').read_text()
        assert parameters == f'K={k}\nreplicate={rep}\nseed={460000+k*100+rep}\nCV_folds=5\n'
        cv = float(re.search(rf'CV error \(K={k}\): ([\d.eE+-]+)', log)[1])
        ll = float(re.findall(r'^Loglikelihood: ([\d.eE+-]+)', log, re.M)[-1])
        q_path = folder/f'primary460.{k}.Q'
        q = np.loadtxt(q_path).reshape(460, k)
        assert np.isfinite(q).all() and (q >= 0).all() and (q <= 1).all()
        assert np.max(np.abs(q.sum(axis=1)-1)) < 2e-5
        checksums[str(q_path.relative_to(ROOT))] = hashlib.sha256(q_path.read_bytes()).hexdigest()
        runs.append(dict(k=k, replicate=rep, seed=460000+k*100+rep, cv=cv,
                         loglikelihood=ll, q=q))
    # Highest likelihood point estimate, with deterministic tie-breaking.
    best = max(runs, key=lambda x: (x['loglikelihood'], -x['replicate']))
    ref = best['q']
    # Preserve previously approved colours by relabelling columns only.
    # No old genotype, Q value or assignment contributes to the new estimates.
    if k == 3:
        with (ROOT/'input/historical_K3_display_label_reference.csv').open() as f:
            old = {row['Accession']:row for row in csv.DictReader(f)}
        assert set(old) == set(ids)
        anchor = np.array([[float(old[i][f'Cluster{j}']) for j in range(1,4)] for i in ids])
        ref, display_perm = align(ref, anchor)
        display_label_alignment = dict(
            K=3, source='input/historical_K3_display_label_reference.csv',
            original_best_run_columns_1based=[int(x+1) for x in display_perm],
            note='Column relabelling solely for colour continuity; all estimates from corrected genotypes.')
    ref_assign = assignment(ref)
    aligned = []
    for run in runs:
        q, perm = align(run['q'], ref)
        aligned.append(q)
        run_rows.append({key: value for key, value in run.items() if key != 'q'} | {
            'reference_replicate': best['replicate'],
            'Q_RMSE_to_reference': float(np.sqrt(np.mean((q-ref)**2))),
            'Q070_assignment_agreement': float(np.mean(assignment(q) == ref_assign)),
            'aligned_columns_original_1based': ';'.join(str(i+1) for i in perm)})
    cv = np.array([r['cv'] for r in runs])
    rmse = [np.sqrt(np.mean((a-b)**2)) for i, a in enumerate(aligned) for b in aligned[i+1:]]
    k_rows.append(dict(K=k, replicates=10, mean_CV=float(cv.mean()), sd_CV=float(cv.std(ddof=1)),
                       min_CV=float(cv.min()), max_CV=float(cv.max()),
                       best_likelihood_replicate=best['replicate'],
                       best_loglikelihood=best['loglikelihood'],
                       mean_pairwise_aligned_Q_RMSE=float(np.mean(rmse)),
                       max_pairwise_aligned_Q_RMSE=float(np.max(rmse)),
                       minimum_assignment_agreement_to_reference=min(r['Q070_assignment_agreement'] for r in run_rows if r['k']==k),
                       assigned_Q070=int((ref_assign > 0).sum()), admixed_Q070=int((ref_assign==0).sum())))
    rows = [dict(Accession=i, **{f'Cluster{c+1}': float(q[c]) for c in range(k)},
                 maxQ=float(q.max()), Assignment=f'Cluster{a}' if a else 'Admixed')
            for i, q, a in zip(ids, ref, ref_assign)]
    write_csv(OUT/f'K{k}_best_Q_assignments.csv', rows)

write_csv(OUT/'replicate_diagnostics.csv', run_rows)
write_csv(OUT/'K_summary.csv', k_rows)
(OUT/'validation.json').write_text(json.dumps(dict(
    sample_count=460, SNP_count=4194, total_runs=100,
    display_label_alignment=display_label_alignment,
    all_samples_match_approved_panel=True, Q_threshold=.70,
    minimum_mean_CV_K=min(k_rows, key=lambda r: r['mean_CV'])['K'],
    selection_note='Minimum mean CV is a diagnostic; K is not a proven number of biological populations.',
    stability_note='Q RMSE and assignment agreement measure replicate optimization stability, not biological confidence.',
    reference_rule='Highest full-data loglikelihood within K; earliest replicate breaks ties.',
    q_sha256=checksums), indent=2)+'\n')
print(json.dumps(k_rows, indent=2))
