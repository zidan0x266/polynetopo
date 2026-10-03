# Cluster validation checklist: examples/shortest_path

Loop detection (`igraph.minimum_cycle_basis()`) on this dataset repeatedly drove a
24GB sandboxed environment into memory thrashing - three attempts (unbounded, and
`cutoff=90`) never completed. Strands and danglings are already validated exactly on
that smaller machine. This checklist runs the remaining stages (loops, Scanlan-Case)
on a machine with real headroom (128GB) and compares against the existing ground
truth in this folder.

## 1. Setup

- [ ] `git pull` (or clone `https://github.com/zidan0x266/polynetopo.git`) to get the
      latest code - needs at least commit `8fcd393` (full LAMMPS data file support).
- [ ] Confirm Python >= 3.10: `python3 --version`
- [ ] Install dependencies: `pip install -e .` from the repo root, or confirm
      `networkx`, `numpy`, `igraph` (the `python-igraph` package), `tqdm` are already
      available: `python3 -c "import networkx, numpy, igraph, tqdm"`

## 2. Sanity check before the big run

- [ ] Confirm the topology file loads and has the expected size:
  ```bash
  cd <repo root>
  python3 -c "
  from graph_analysis.network import PolymerNetwork
  n = PolymerNetwork.from_lammps('examples/shortest_path/atrp_s3_sys3_0.data')
  print(n)  # expect: PolymerNetwork(nodes=104704, edges=114029)
  "
  ```

## 3. Run loop detection (the step that failed locally)

```bash
python3 src/main.py --function topology --stage loops \
    --input-file examples/shortest_path/atrp_s3_sys3_0.data \
    --output-file igloops_new.pkl
```

- [ ] Completes without error
- [ ] If you can, note wall-clock time and peak memory (e.g. `/usr/bin/time -v ...` on
      Linux, or your job scheduler's reported peak RSS) - useful to know for next time
- [ ] Sanity check the loop count matches the old ground truth (9326):
  ```bash
  python3 -c "
  import pickle
  with open('igloops_new.pkl','rb') as f: loops, orders = pickle.load(f)
  print(len(loops))
  "
  ```

## 4. Run Scanlan-Case refinement

```bash
python3 src/main.py --function topology --stage scanlan-case \
    --input-file examples/shortest_path/atrp_s3_sys3_0.data \
    --strands-file examples/shortest_path/strands.pkl \
    --danglings-file examples/shortest_path/danglings.pkl \
    --loops-file igloops_new.pkl \
    --output-file scanlan_case_new.pkl --raw
```

- [ ] Completes without error (this step is fast - seconds, not minutes)

## 5. Compare against the ground truth in this folder

```bash
cd examples/shortest_path
python3 compare_new_vs_old.py --loops ../../igloops_new.pkl --scanlan-case ../../scanlan_case_new.pkl
```

Expected output:
- [ ] `strands`, `danglings`, `loops`, `loop_orders` (raw): **exact match: True**
- [ ] Scanlan-Case `stats exact match: True`
- [ ] `refined danglings exact match: True`, `refined loop_orders exact match: True`
- [ ] `refined strands`: same count, and a **small** number differing (the
      intentional strand-merge bug fix - see `src/graph_analysis/scanlan_case.py`'s
      module docstring; on `examples/atrp` this was 13 of 4632, expect a similarly
      small fraction here, not zero and not a large fraction)

Anything outside those expectations (a `False` where `True` is expected, or a large
fraction of strands differing) means something is actually wrong - stop and report
back rather than continuing to step 6.

## 6. (Optional) Full one-shot pipeline as an independent cross-check

Since this machine has plenty of memory, running the whole pipeline fresh in one shot
is a good independent check that doesn't reuse anything computed earlier:

```bash
python3 src/main.py --function topology \
    --input-file examples/shortest_path/atrp_s3_sys3_0.data \
    --output-file results_new.pkl --raw --ncpus <N> --max-chunk-size 10000

cd examples/shortest_path
python3 compare_new_vs_old.py --analyze-output ../../results_new.pkl --analyze-raw ../../results_new_raw.pkl
```

(`<N>`: pick based on cores allocated to your job; omit `--ncpus` to default to
`cpu_count() - 2`.)

## 7. Report back

Share:
- [ ] The full printed output of `compare_new_vs_old.py`
- [ ] Wall-clock time and peak memory for the `loops` step, if available
- [ ] Any warnings or errors encountered
