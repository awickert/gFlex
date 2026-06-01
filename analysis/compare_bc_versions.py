#!/usr/bin/env python3
"""
Compare zero_displacement_zero_slope BC accuracy across two git commits.

Uses git worktree to check out an older commit alongside the current
working tree, runs the same MMS convergence study with both, and
reports the error improvement quantitatively.

Usage
-----
    python analysis/compare_bc_versions.py [--old HASH]

    --old HASH   Git commit to treat as "old" (default: parent of the
                 fix commit 6f68270, i.e. the last commit before the
                 ghost-node correction landed).

The script must be run from the repository root.  It creates a
temporary git worktree, runs analysis/_mms_runner.py twice (once with
the current gflex on sys.path, once with the old one prepended via
PYTHONPATH), then removes the worktree and reports.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

# Hash of the commit that introduced the BC fix.  Its parent is the
# last state where the bug was present.
FIX_COMMIT = "6f68270"
DEFAULT_OLD = f"{FIX_COMMIT}^"

RUNNER = os.path.join(os.path.dirname(__file__), "_mms_runner.py")


def run_runner(extra_env=None):
    """Run _mms_runner.py and return parsed JSON dict."""
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)
    result = subprocess.run(
        [sys.executable, RUNNER],
        capture_output=True,
        text=True,
        env=env,
    )
    if result.returncode != 0:
        print("Runner stderr:", result.stderr, file=sys.stderr)
        sys.exit(f"Runner failed with exit code {result.returncode}")
    return json.loads(result.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--old", default=DEFAULT_OLD, metavar="HASH",
                        help=f"Old commit to compare against (default: {DEFAULT_OLD})")
    args = parser.parse_args()

    # Resolve the old commit to a full hash for display
    old_hash = subprocess.check_output(
        ["git", "rev-parse", "--short", args.old], text=True
    ).strip()
    new_hash = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], text=True
    ).strip()

    print(f"Comparing:")
    print(f"  old  {old_hash}  ({args.old})")
    print(f"  new  {new_hash}  (HEAD)")
    print()

    # --- Run with current (new) gflex ---
    print("Running MMS with current (new) gflex ...", flush=True)
    data_new = run_runner()

    # --- Set up old worktree and run ---
    worktree = tempfile.mkdtemp(prefix="gflex-old-")
    try:
        subprocess.run(
            ["git", "worktree", "add", "--detach", worktree, args.old],
            check=True, capture_output=True,
        )
        print(f"Running MMS with old gflex ({old_hash}) ...", flush=True)
        pythonpath = worktree
        if "PYTHONPATH" in os.environ:
            pythonpath += os.pathsep + os.environ["PYTHONPATH"]
        data_old = run_runner(extra_env={"PYTHONPATH": pythonpath})
    finally:
        subprocess.run(
            ["git", "worktree", "remove", "--force", worktree],
            capture_output=True,
        )
        shutil.rmtree(worktree, ignore_errors=True)

    # --- Report ---
    nx_vals  = data_new["nx_vals"]
    dx_km    = data_new["dx_km"]
    errs_new = data_new["errors"]
    errs_old = data_old["errors"]

    print()
    print("=" * 66)
    print(f"  MMS error: zero_displacement_zero_slope")
    print(f"  old={old_hash} ({data_old['gflex_version']})  "
          f"vs  new={new_hash} ({data_new['gflex_version']})")
    print("=" * 66)
    print(f"  {'nx':>5}  {'dx [km]':>8}  {'old error':>12}  {'new error':>12}  {'factor':>8}")
    print(f"  {'-'*5}  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*8}")
    for nx, dx, e_old, e_new in zip(nx_vals, dx_km, errs_old, errs_new):
        factor = e_old / e_new
        print(f"  {nx:>5}  {dx:>8.2f}  {e_old:>12.4e}  {e_new:>12.4e}  {factor:>8.1f}x")
    print()
    print(f"  Convergence order:")
    print(f"    old  O(dx^{data_old['convergence_slope']:.2f})")
    print(f"    new  O(dx^{data_new['convergence_slope']:.2f})")
    print()
    print(f"  w[0] at finest grid (should be 0):")
    print(f"    old  {data_old['w0_at_finest']:+.4e} m")
    print(f"    new  {data_new['w0_at_finest']:+.4e} m")
    print()

    # --- Optional convergence plot ---
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.loglog(dx_km, errs_old, "r^--", label=f"old  {old_hash}  "
                  f"O(dx^{data_old['convergence_slope']:.1f})")
        ax.loglog(dx_km, errs_new, "bs-",  label=f"new  {new_hash}  "
                  f"O(dx^{data_new['convergence_slope']:.1f})")
        ref = np.array([min(dx_km), max(dx_km)])
        ax.loglog(ref, 5e-3 * (ref / ref[-1]) ** 1, "k--", lw=0.8, label="O(dx¹)")
        ax.loglog(ref, 5e-5 * (ref / ref[-1]) ** 2, "k:",  lw=0.8, label="O(dx²)")
        ax.set_xlabel("dx [km]")
        ax.set_ylabel("L-inf relative error")
        ax.set_title("zero_displacement_zero_slope: MMS convergence\n"
                     "old vs corrected ghost-node treatment")
        ax.legend(fontsize=8)
        fig.tight_layout()
        out = "analysis/results/compare_clamped_bc.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"  Figure saved to {out}")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
