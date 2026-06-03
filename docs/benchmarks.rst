Performance Benchmarks
======================

The figures on this page were produced by running ``benchmarks/bench_solvers.py``
from the repository root.  To regenerate the figures from a fresh benchmark run::

   python benchmarks/bench_solvers.py
   python benchmarks/plot_results.py

Hardware context for the results shown here:

* **CPU**: AMD Ryzen 7 7840U (16 logical cores)
* **RAM**: 46.4 GiB
* **OS**: Linux 6.17, Ubuntu 24.04
* **Python** 3.12.7 · **NumPy** 1.26.4 · **SciPy** 1.16.3

Absolute timings will differ on other hardware; relative performance between
methods and cache modes is broadly representative.

Solver scaling
--------------

.. figure:: _static/bench_scaling.png
   :width: 100%
   :alt: Solve time vs. grid size for FD, FFT, and SAS (1-D and 2-D)

The **FD direct solver** uses a sparse LU factorisation and scales roughly as
:math:`O(N^{1.5})` in two dimensions, where :math:`N` is the total number of
grid cells.  At 400×400 cells (:math:`N^2 = 160\,000`) a single FD solve takes
about 5.5 s; at 200×200 it takes about 0.7 s.

The **FFT** method solves the problem in the spectral domain and scales as
:math:`O(N \log N)`.  In two dimensions it is roughly two orders of magnitude
faster than FD at the same grid size, but requires a uniform (scalar) elastic
thickness and imposes periodic boundary conditions (or zero-padded behaviour
with ``no_outside_loads``).

The **SAS** method superposes analytical deflection kernels using
``fftconvolve``, scaling as :math:`O(N^2 \log N)` in two dimensions.  It is
load-pattern-independent and faster than FD for small to moderate grids,
but becomes comparable to or slower than FD at very large grids.  It also
requires a uniform elastic thickness.

The Te profile has a modest effect on FD timing: the shaded band (min to max
across all profiles tested — constant, sinusoidal, abrupt step, tanh sigmoidal,
correlated noise, wide dynamic range) is narrow, typically within ±15 % of the
mean.  Spatially noisy profiles show slightly higher times at large grid sizes,
likely due to increased fill-in during sparse LU factorisation.

LU factorisation cache
-----------------------

.. figure:: _static/bench_lu_cache.png
   :width: 100%
   :alt: LU cache speedup via the run() path

   **Speedup of the ``True`` and ``"no_check"`` cache modes over uncached
   (``False``) via the full** :meth:`~gflex.F1D.run` **path.**
   The band shows the range across all Te profiles; the line is the mean.
   The dotted line at 1× marks no speedup.

When the load :math:`q_s` changes between calls but the grid, elastic
thickness, and boundary conditions remain fixed, the coefficient matrix does
not change.  Caching the LU factorisation eliminates the cost of re-factorising
on every call, reducing the per-solve work to a single triangular solve.

In two dimensions the **``"no_check"``** mode reaches **6–10× speedup** over
uncached at grid sizes of 50×50 to 200×200 cells.  The speedup grows with
grid size because the factorisation cost (eliminated by caching) scales as
:math:`O(N^{1.5})` while the cached solve (triangular back-substitution) scales
as :math:`O(N)`.

The **``True``** mode (hash-validated cache reuse) provides a similar trend
but at a lower speedup because it computes and compares a matrix hash on every
call.  For large 2-D grids the hash cost is a small fraction of the total, so
``True`` and ``"no_check"`` converge.  For 1-D problems, where solve times are
short, the hash overhead is comparatively larger.

The timings include all overhead from a ``run()`` call: :meth:`bc_check`,
coordinate setup, warning checks, and the :meth:`_solve_fd` cache-bypass logic.
This is more representative of real coupling-loop performance than timing
:meth:`fd_solve` in isolation.

See :doc:`api` for usage details, including the smart invalidation mechanism
that automatically clears the cache when any matrix-determining input (``te``,
``dx``, boundary conditions, etc.) is reassigned.

Cost of changing :math:`T_e`
-----------------------------

.. figure:: _static/bench_te_sweep.png
   :width: 100%
   :alt: Per-solve cost: load-only vs Te change

   **Per-solve cost when only the load changes (load-only, ``"no_check"``
   mode) versus when the scalar elastic thickness :math:`T_e` changes on
   every call (Te change, any cache mode).**

Reassigning :math:`T_e` triggers smart cache invalidation: the coefficient
matrix is cleared and the LU factorisation is discarded.  The next
:meth:`~gflex.F1D.run` call rebuilds the matrix from scratch and re-factorises.
All three cache modes (``False``, ``True``, ``"no_check"``) pay essentially the
same cost when :math:`T_e` changes, because the rebuild dominates the per-call
budget.

At 200×200 cells the load-only cost is roughly **70 ms** per solve while a
:math:`T_e` change costs roughly **550 ms** per solve — about an **8× penalty**.
The gap widens with grid size because factorisation scales as :math:`O(N^{1.5})`
and the incremental triangular solve scales as :math:`O(N)`.

The practical implication: in a coupling loop where :math:`T_e` is fixed and
only :math:`q_s` varies (e.g., a transient ice-sheet or sediment-loading model),
``cache_factorization = "no_check"`` or ``True`` provides substantial speedup.
In a parameter-sweep or inversion where :math:`T_e` changes on every iteration,
all modes are equivalent and the rebuild cost is unavoidable.
