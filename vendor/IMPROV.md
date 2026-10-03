# Building the in-house igraph on LCRC Improv

Step-by-step for building `vendor/igraph`/`vendor/python-igraph` on Argonne LCRC's
Improv cluster. Verified against the current LCRC user guide (PBS Pro job syntax,
queue names, module system); a couple of module *names/versions* are flagged below
as things to confirm yourself with `module avail` since they can change over time -
don't take those two as gospel, everything else here is from the current docs.

## 0. Before you start

- You need an Improv account and an allocation project name (for `-A` below) - find
  it with `sbank-list-allocations -p <project>` once logged in, or ask your PI/LCRC
  if you don't have one.
- **Don't build on the login node.** Compiling igraph is real CPU work; LCRC's
  etiquette (like most HPC centers) is that login nodes are for editing/light work
  only. Step 2 below gets you an interactive compute node first.

## 1. Log in and get the code

```bash
ssh <your-username>@improv.lcrc.anl.gov
git clone https://github.com/zidan0x266/polynetopo.git   # or `git pull` if already cloned
cd polynetopo
git submodule update --init vendor/igraph vendor/python-igraph
```

## 2. Request an interactive compute node

Improv uses **PBS Pro**, not SLURM. A short interactive session on the standard
`compute` queue (128 cores/node, 256GB RAM) is enough to build on:

```bash
qsub -I -A <project> -q compute -l select=1:ncpus=16,walltime=01:00:00
```

- `-A <project>`: your allocation project name (mandatory).
- `-l select=1:ncpus=16`: one node, 16 cores - plenty for compiling; you don't need a
  full 128-core node just to build (save the full node for an actual production run,
  step 6).
- `-q compute` is the default/standard queue; `debug` (1 hour max) also works fine for
  a build this size. `walltime=01:00:00` (1 hour) is generous headroom - the whole
  build+test suite took ~15-20 minutes on a Mac laptop; expect similar or faster here.

Wait for the prompt to return on a compute node (hostname changes), then continue
from there - everything below runs on that compute node, in your `polynetopo` clone.

## 3. Load the toolchain

```bash
module avail gcc                 # confirm a recent gcc, e.g. gcc/11.4.0 is known to exist on Improv
module load gcc                  # or `module load gcc/<version>` for a specific one
module avail cmake bison flex pkg-config python anaconda3
```

The build needs Python (3.11+; this project has been validated with 3.13.15), a C/C++
toolchain, CMake >= 3.18, Make, Bison, Flex and pkg-config. `gcc` is confirmed
available as a module on Improv; **confirm the exact module names/versions for
cmake/bison/flex/pkg-config/python yourself** with `module avail <name>` - LCRC's
module set changes over time and these weren't individually confirmed when this
tutorial was written. Load whatever versions `module avail` shows you, e.g.:

```bash
module load cmake bison flex pkg-config
module load anaconda3   # or whichever module provides your Python 3
```

Verify before proceeding:

```bash
gcc --version && cmake --version && bison --version && flex --version && pkg-config --version
python3 --version
```

If any of these aren't found, they may need `spack load <package>` instead of
`module load` (LCRC uses Spack under the hood) - try `spack find <name>` if a plain
`module avail` comes up empty.

## 4. Build

From the repo root (on the compute node):

```bash
# The fork has no git tags; igraph's CMake build needs one or this file instead:
echo "1.0.0" > vendor/igraph/IGRAPH_VERSION

cd vendor
SPATIAL_BUILD_JOBS=16 bash python-igraph/scripts/build_spatial.sh
cd ..
```

`SPATIAL_BUILD_JOBS=16` matches the 16 cores requested in step 2 (compilation
parallelism only - the centrality calculation itself is always single-threaded,
regardless of this setting). This creates `vendor/.spatial-venv` and
`vendor/.spatial-build`, builds the C core, runs its test suite (including the new
`igraph_edge_betweenness_spatial` tests), builds the Python extension against it, and
runs the Python test suite. Expect it to print `64 passed, ...` (or similar) at the
end with no errors.

If you hit a *different* CMake version error than the tag one already fixed above,
or a missing-header/library error, it's almost certainly a module gap from step 3 -
re-check `module avail`/`spack find` for the missing tool rather than assuming the
build script itself is broken (it's already validated working, just not yet
specifically on Improv as of this writing).

## 5. Finish the environment and validate

```bash
vendor/.spatial-venv/bin/python -m pip install --no-build-isolation --no-deps networkx tqdm
vendor/.spatial-venv/bin/python -m pip install -e .

vendor/.spatial-venv/bin/python -c 'import igraph as ig; print(hasattr(ig.Graph, "edge_betweenness_spatial"))'
# must print True
```

Quick functional check with a real file:

```bash
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/betweenness/atrp_s1_sys1_0.data --inspect
```

Once this all works, `exit` the interactive session (`qsub -I` sessions count
against your walltime/allocation the whole time they're open).

## 6. Running a real analysis as a batch job

For anything beyond a quick interactive check, submit a proper batch job rather than
tying up an interactive session. Example PBS script
(`run_betweenness.pbs`) requesting a full node:

```bash
#!/bin/bash
#PBS -N betweenness
#PBS -A <project>
#PBS -q compute
#PBS -l select=1:ncpus=128
#PBS -l walltime=02:00:00
#PBS -j oe

cd "$PBS_O_WORKDIR"
vendor/.spatial-venv/bin/python src/main.py --function betweenness \
    --input-file examples/betweenness/atrp_s1_sys1_0.data \
    --direction x --validate-stock
```

```bash
qsub run_betweenness.pbs
pbsq -f <your-username>       # or: qstat -was1
```

Adjust `--l select=1:ncpus=...` and `walltime` to the job - `--function betweenness`
is single-threaded (more cores don't help it), while `--function topology`/
`shortest-path` do use multiple processes (`--ncpus`) for the crosslinker-pair search.

See `vendor/README.md` for what the build actually does and the underlying
`vendor/python-igraph/SPATIAL_BETWEENNESS.md` for the full scientific/build reference,
and the three `examples/*/TUTORIAL.md` files for what to actually run once this is
built.

Sources:
- [Getting Started on LCRC's Improv Cluster](https://www.anl.gov/event/getting-started-on-lcrcs-improv-cluster)
- [Running Jobs on Improv - LCRC User Guide](https://argonne-lcrc.github.io/user-guides/improv/running-jobs-improv/)
- [Running Jobs on PBS Pro Clusters - LCRC User Guide](https://argonne-lcrc.github.io/user-guides/running-jobs-at-lcrc/pbs-pro/)
