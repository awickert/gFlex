"""Blender mesh export for gFlex deflection results.

Writes a pure-Python mesh file that can be consumed by a companion
Blender script to build a 3-D render of a flexural deflection surface.
Because Blender's embedded Python interpreter does not have NumPy or
gFlex available, all data are serialised as plain Python lists.

See ``docs/examples/blender_flexure.py`` for a ready-to-use Blender
scene script that reads the exported file.
"""

import numpy as np


def export_for_blender(
    w,
    dx,
    dy=None,
    path="gflex_blender_mesh.py",
    z_exaggeration=1000.0,
    plate_size_bu=10.0,
    qs=None,
    Te=None,
    rho_m=3300.0,
    g=9.8,
):
    """Export a gFlex deflection surface as a Blender-ready mesh file.

    Writes a pure-Python file (no NumPy imports) containing vertex
    coordinates, quad-face indices, and optional per-vertex load and
    elastic-thickness attributes.  The file is designed to be
    ``exec()``'d directly inside Blender's embedded Python interpreter,
    which does not have NumPy or gFlex available.

    Pair this file with the companion Blender script in
    ``docs/examples/blender_flexure.py`` to produce publication-quality
    3-D renders of flexure surfaces.

    Parameters
    ----------
    w : array_like, shape (ny, nx)
        Deflection array [m].  Negative values are downward (the
        standard gFlex sign convention).
    dx : float
        Grid spacing in the x-direction [m].
    dy : float, optional
        Grid spacing in the y-direction [m].  Defaults to *dx*.
    path : str, optional
        Output file path for the mesh file.
        Default ``'gflex_blender_mesh.py'``.
    z_exaggeration : float, optional
        Vertical exaggeration factor.  The deflection is scaled so that
        a displacement of ``plate_width_m / z_exaggeration`` metres maps
        to 1 Blender unit.  Default 1000, appropriate for
        continental-scale flexure where deflections are of order km and
        the domain is of order 1000 km.  Increase for subtler signals;
        decrease for dramatic relief.
    plate_size_bu : float, optional
        Physical width of the plate along x in Blender units.
        Default 10.0.
    qs : array_like, optional
        Surface load array [Pa], same shape as *w*.  When supplied,
        per-vertex ``qs_val`` (raw Pa) and ``qs_norm`` (normalised to
        [0, 1]) attributes are written to the mesh file.  If *rho_m*
        and *g* are also provided, the load cylinder geometry
        (``load_radius_bu``, ``load_height_bu``, ``load_height_m_eq``,
        ``base_z_bu``) is written for use by the companion Blender script.
    Te : array_like or float, optional
        Elastic thickness [m], either a scalar or an array of the same
        shape as *w*.  When supplied, per-vertex ``te_val`` attributes
        and the scalar ``te_min`` / ``te_max`` bounds are written to the
        mesh file.
    rho_m : float, optional
        Mantle density [kg m⁻³] used to convert load pressure to
        mantle-equivalent height for the extruded load cylinder.
        Default 3300.
    g : float, optional
        Gravitational acceleration [m s⁻²].  Default 9.8.

    Returns
    -------
    None
        Writes *path* to disk and prints a confirmation line.

    Notes
    -----
    **Coordinate system** — the plate is centred on the Blender origin.
    One Blender unit equals ``plate_width_m / plate_size_bu`` metres in
    the horizontal plane.  The vertical axis is scaled by the same
    factor multiplied by *z_exaggeration*, making features that are
    invisible at true scale perceptible in the render.

    **Variables written to the mesh file:**

    * ``nrows``, ``ncols`` — grid dimensions.
    * ``dx_m``, ``dy_m`` — original grid spacing [m].
    * ``w_min_bu``, ``w_max_bu`` — deflection extremes in Blender units.
    * ``clip_bu`` — symmetric colormap clip
      (``max(|w_min_bu|, |w_max_bu|)``); use as ``vmin = -clip_bu`` and
      ``vmax = clip_bu`` when driving a diverging (e.g. vik) colour ramp.
    * ``plate_size_bu`` — x-width of the plate in Blender units.
    * ``vertices`` — list of ``(bx, by, bz)`` tuples.
    * ``faces`` — list of ``(v0, v1, v2, v3)`` quad-face index tuples.
    * ``qs_verts``, ``qs_norm`` — per-vertex load values *(if qs given)*.
    * ``load_radius_bu`` — footprint radius of the load in Blender units
      *(if qs given and load is approximately circular)*.
    * ``load_height_m_eq`` — mantle-equivalent height of the load [m]
      (``qs_max / (rho_m * g)``).
    * ``load_height_bu`` — load height in Blender units.
    * ``base_z_bu`` — z-coordinate of the plate surface at its centre [BU],
      used as the base of the extruded load cylinder.
    * ``te_verts``, ``te_min``, ``te_max`` — per-vertex Te *(if Te given)*.

    Examples
    --------
    Run gFlex, then export for Blender::

        import numpy as np
        from gflex import F2D, export_for_blender

        flex = F2D()
        flex.quiet = True
        flex.method = 'fd'
        flex.solver = 'direct'
        flex.g = 9.8;  flex.E = 65e9;  flex.nu = 0.25
        flex.rho_m = 3300.;  flex.rho_fill = 0.
        flex.T_e = 35e3 * np.ones((100, 100))
        flex.qs = np.zeros((100, 100));  flex.qs[40:60, 40:60] = 1e6
        flex.dx = flex.dy = 5000.
        flex.bc_west = flex.bc_east = flex.bc_north = flex.bc_south = 'zero_moment_zero_shear'
        flex.initialize();  flex.run();  flex.finalize()

        export_for_blender(
            flex.w, flex.dx,
            path='/tmp/flexure_mesh.py',
            qs=flex.qs,
        )
        # → open Blender and run docs/examples/blender_flexure.py
    """
    if dy is None:
        dy = dx

    w = np.asarray(w, dtype=float)
    ny, nx = w.shape

    # Horizontal scale: the physical x-width maps to plate_size_bu BU
    plate_width_m = nx * dx
    m_per_bu = plate_width_m / plate_size_bu
    half_x = plate_size_bu / 2.0
    # Preserve aspect ratio in y
    half_y = (ny * dy) / m_per_bu / 2.0

    # Vertical scale: z_exaggeration applied on top of horizontal scale
    z_scale = z_exaggeration / m_per_bu

    w_bu = w * z_scale
    w_min_bu = float(w_bu.min())
    w_max_bu = float(w_bu.max())

    # Asymmetric clip: when the forebulge is small compared to the subsidence
    # give the positive arm enough range so it is visible on the colour ramp.
    # When the two sides are comparable, fall back to symmetric.
    if w_max_bu > 0 and w_max_bu < 0.2 * abs(w_min_bu):
        vmax_bu = max(w_max_bu, 0.1 * abs(w_min_bu))
    else:
        vmax_bu = max(w_max_bu, abs(w_min_bu))
    clip_bu = max(abs(w_min_bu), vmax_bu)   # kept for backward compat

    # Pre-normalise deflection to [0, 1] with zero mapped to 0.5
    # (TwoSlopeNorm equivalent; stored as a vertex attribute in the mesh file)
    def _w_norm(val):
        if val <= 0.0:
            return 0.5 * (val - w_min_bu) / (0.0 - w_min_bu) if w_min_bu < 0 else 0.5
        else:
            return 0.5 + 0.5 * val / vmax_bu

    # Cell-centre coordinates in Blender units
    xi = (np.arange(nx) + 0.5) * dx
    yi = (np.arange(ny) + 0.5) * dy
    bx = (xi / plate_width_m) * plate_size_bu - half_x
    by = (yi / (ny * dy)) * (2.0 * half_y) - half_y

    # Vertex list: row-major (j = row = y, i = col = x)
    verts = [
        (float(bx[i]), float(by[j]), float(w_bu[j, i]))
        for j in range(ny)
        for i in range(nx)
    ]
    w_norm = [_w_norm(float(w_bu[j, i])) for j in range(ny) for i in range(nx)]

    # Quad faces with consistent counter-clockwise winding
    faces = [
        (j * nx + i, j * nx + i + 1,
         (j + 1) * nx + i + 1, (j + 1) * nx + i)
        for j in range(ny - 1)
        for i in range(nx - 1)
    ]

    # Optional per-vertex load attribute and cylinder geometry
    qs_verts = qs_norm = None
    load_radius_bu = load_height_bu = load_height_m_eq = base_z_bu = None
    if qs is not None:
        qs_arr = np.asarray(qs, dtype=float)
        qs_max = float(qs_arr.max()) if qs_arr.max() > 0.0 else 1.0
        qs_verts = [float(qs_arr[j, i]) for j in range(ny) for i in range(nx)]
        qs_norm  = [v / qs_max for v in qs_verts]

        # Load cylinder geometry — infer radius from the loaded footprint
        loaded = qs_arr > 0.0
        if loaded.any():
            # Approximate circular radius from loaded area
            load_area_m2 = float(loaded.sum()) * dx * dy
            load_radius_m = np.sqrt(load_area_m2 / np.pi)
            load_radius_bu = float(load_radius_m / m_per_bu)
        else:
            load_radius_bu = 0.0

        load_height_m_eq = float(qs_max / (rho_m * g))
        load_height_bu   = float(load_height_m_eq * z_scale)

        # Centre of the plate in Blender units — z of nearest vertex to (0, 0)
        cx_bu = 0.0;  cy_bu = 0.0
        cx_idx = int(nx / 2);  cy_idx = int(ny / 2)
        base_z_bu = float(w_bu[cy_idx, cx_idx])

    # Optional per-vertex elastic-thickness attribute
    te_verts = te_min_val = te_max_val = None
    if Te is not None:
        Te_arr = np.broadcast_to(np.asarray(Te, dtype=float), (ny, nx))
        te_verts   = [float(Te_arr[j, i]) for j in range(ny) for i in range(nx)]
        te_min_val = float(Te_arr.min())
        te_max_val = float(Te_arr.max())

    # Write pure-Python file (no imports — exec()-safe inside Blender)
    with open(path, "w") as fh:
        fh.write(f"# gFlex Blender mesh export\n")
        fh.write(f"# Grid: {ny} rows × {nx} cols  |  "
                 f"dx={dx:.1f} m, dy={dy:.1f} m\n")
        fh.write(f"# Deflection: {w.min():+.2f} to {w.max():+.2f} m  "
                 f"(z_exaggeration={z_exaggeration:.0f})\n\n")
        fh.write(f"nrows         = {ny}\n")
        fh.write(f"ncols         = {nx}\n")
        fh.write(f"dx_m          = {dx}\n")
        fh.write(f"dy_m          = {dy}\n")
        fh.write(f"w_min_bu      = {w_min_bu}\n")
        fh.write(f"w_max_bu      = {w_max_bu}\n")
        fh.write(f"vmax_bu       = {vmax_bu}\n")
        fh.write(f"clip_bu       = {clip_bu}\n")
        fh.write(f"plate_size_bu = {plate_size_bu}\n")
        fh.write(f"vertices      = {verts!r}\n")
        fh.write(f"faces         = {faces!r}\n")
        fh.write(f"w_norm        = {w_norm!r}\n")
        if qs_verts is not None:
            fh.write(f"qs_verts      = {qs_verts!r}\n")
            fh.write(f"qs_norm       = {qs_norm!r}\n")
            if load_radius_bu is not None:
                fh.write(f"load_radius_bu    = {load_radius_bu}\n")
                fh.write(f"load_height_m_eq  = {load_height_m_eq}\n")
                fh.write(f"load_height_bu    = {load_height_bu}\n")
                fh.write(f"base_z_bu         = {base_z_bu}\n")
        if te_verts is not None:
            fh.write(f"te_verts      = {te_verts!r}\n")
            fh.write(f"te_min        = {te_min_val}\n")
            fh.write(f"te_max        = {te_max_val}\n")

    print(f"Blender mesh written → {path}  "
          f"({len(verts)} verts, {len(faces)} quads)")
