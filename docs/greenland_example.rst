Greenland Glacial Isostasy
==========================

This example applies gFlex to a real geophysical problem: the deflection of
the Greenland lithosphere under the present-day ice sheet.  It demonstrates
how to integrate external datasets, handle coordinate reprojection, apply
domain padding, and interpret the role of spatially variable elastic thickness.

The full script is at ``examples/greenland_flexure.py``.

Datasets
--------

**Ice thickness** — BedMachine Greenland v6 (Morlighem et al. 2017), available
from NSIDC (Earthdata login required).  The native 150 m grid is subsampled to
~10 km before passing to gFlex; at that resolution each FD solve takes a few
seconds.

**Elastic thickness** — Steffen et al. (2018, 2026), freely available under
CC-BY-4.0 from Zenodo
(`doi:10.5281/zenodo.18403685 <https://doi.org/10.5281/zenodo.18403685>`_).
The script downloads this file automatically on first run.  The dataset provides
:math:`T_e` in km at 0.2° geographic spacing over the greater Greenland region;
values range from 5 to 87 km, reflecting the contrast between the weak western
margin and the cratonic interior.

Setup
-----

The BedMachine grid is already in EPSG:3413 (NSIDC Sea Ice Polar Stereographic
North).  The Steffen et al. data is in geographic lon/lat and must be
reprojected onto the same EPSG:3413 grid using :mod:`pyproj` and
:class:`scipy.interpolate.RegularGridInterpolator`.  Cells outside the Te
dataset coverage are assigned the domain median.

Only **grounded ice** (BedMachine mask == 2) contributes to the surface load;
floating ice is hydrostatically supported by the ocean and is excluded.  The
load is :math:`q_s = \rho_\text{ice}\,g\,H_\text{ice}` with
:math:`\rho_\text{ice} = 917` kg m⁻³.

Boundary conditions and domain padding
---------------------------------------

A clamped boundary condition (``zero_displacement_zero_slope``) is applied on
all four sides.  This requires that deflections approach zero before the domain
edge — enforced here by calling :func:`~gflex.pad_domain` with one flexural
wavelength of padding (~650 km for the mean :math:`T_e` of 46 km).  The
padding ring tapers :math:`T_e` smoothly to the domain mean and carries zero
load; after the solve the padding is trimmed from :attr:`~gflex.F2D.w` with
``w = flex.w[p:-p, p:-p]``.

Without padding the forebulge reached 61 m; with padding it is 36 m,
confirming that free-end boundary artefacts were inflating the unpadded result.

.. code-block:: python

   from gflex import F2D, pad_domain

   te_pad, qs_pad, p = pad_domain(te_proj, qs, dx=dx, dy=dy,
                                   E=65e9, nu=0.25, rho_m=3300., g=9.81)
   flex = F2D()
   flex.te, flex.qs = te_pad, qs_pad
   flex.dx, flex.dy = dx, dy
   flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = \
       'zero_displacement_zero_slope'
   flex.initialize()
   flex.run()
   w = flex.w[p:-p, p:-p]   # trim padding

Results
-------

.. figure:: _static/greenland_flexure.png
   :width: 100%
   :alt: Greenland Te, ice thickness, and lithospheric deflection

   *Left*: effective elastic thickness :math:`T_e` (Steffen et al.;
   roma colormap).  *Centre*: ice thickness from BedMachine v6 (grounded ice
   only; devon colormap).  *Right*: computed lithospheric deflection (vik
   colormap, asymmetric scale: blue = subsidence, red = uplift).

The maximum depression is **880 m** beneath the central ice dome.  The
forebulge reaches **36 m** and is notably sharp along the east coast, where
:math:`T_e` drops steeply from the cratonic interior toward the passive margin.
That abrupt change in plate stiffness concentrates the bending and narrows the
forebulge wavelength.

Effect of spatially variable :math:`T_e`
------------------------------------------

Running a second solve with :math:`T_e` fixed at the ice-sheet mean
(46 km everywhere) isolates the contribution of lateral heterogeneity.

.. figure:: _static/greenland_flexure_comparison.png
   :width: 100%
   :alt: Uniform vs variable Te deflection and their difference

   *Left*: deflection under uniform :math:`T_e` = 46 km.
   *Centre*: deflection under the spatially variable Steffen et al.
   :math:`T_e`.  *Right*: difference (variable − uniform); red = variable
   :math:`T_e` produces more subsidence (locally weak lithosphere), blue =
   less subsidence (locally strong lithosphere).

The two total depressions agree within 1 m (both ~880 m), because the load
magnitude is dominated by the thick central ice, which sits on near-average
:math:`T_e`.  The differences are largest at the margins: up to **+105 m**
(more subsidence) where :math:`T_e` dips well below the mean along the western
coast, and up to **−85 m** (less subsidence) over the strong cratonic keel in
the southeast.  For a regional model used to infer ice-loading history or
present-day rebound rates, assuming a uniform :math:`T_e` would introduce
errors of this magnitude at the margins.

References
----------

Morlighem M. et al. (2017), BedMachine v3: Complete bed topography and ocean
bathymetry mapping of Greenland from multi-beam echo sounding combined with
mass conservation, *Geophys. Res. Lett.*, 44, 11051–11061,
`doi:10.1002/2017GL074954 <https://doi.org/10.1002/2017GL074954>`_.

Steffen R., Audet P., and Lund B. (2018), Weakened lithosphere beneath
Greenland inferred from effective elastic thickness: A hot spot effect?,
*Geophys. Res. Lett.*, 45(10), 4733–4742,
`doi:10.1029/2017GL076885 <https://doi.org/10.1029/2017GL076885>`_.

Steffen R., Audet P., Strykowski G., and Lund B. (2026), Greenland — Moho and
effective elastic thickness models based on satellite gravity data, Zenodo,
`doi:10.5281/zenodo.18403685 <https://doi.org/10.5281/zenodo.18403685>`_.
