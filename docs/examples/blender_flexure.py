"""Blender 4 scene script for gFlex deflection renders.

Run headless::

    blender --background --python docs/examples/blender_flexure.py

Reads a mesh file produced by :func:`gflex.export_for_blender` and
builds a 3-D scene: deflection surface coloured with the *vik* diverging
colour ramp, an optional wireframe grid coloured by elastic thickness,
and an optional load surface.  Renders to a PNG with a transparent
background suitable for compositing into papers or presentations.

Quick-start example (reproduces the gFlex logo scene)::

    # Step 1 — generate the flexure data and mesh file
    python docs/logo/generate_logo_data.py

    # Step 2 — export to Blender format
    python - <<'EOF'
    import numpy as np
    from gflex import export_for_blender
    data = np.load("/tmp/gflex_logo_data.npz")
    export_for_blender(
        data["w"], float(data["dx"]),
        path="/tmp/gflex_blender_mesh.py",
        z_exaggeration=1000.0,
        qs=data["qs"],
        Te=data["Te"],
    )
    EOF

    # Step 3 — render
    blender --background --python docs/examples/blender_flexure.py

This script is intentionally kept as a self-contained template.  Copy
it, adjust the USER CONFIG section at the top, and adapt the scene
builders to suit your figure.
"""

import math
import os

import bpy
import mathutils

# ── USER CONFIG ───────────────────────────────────────────────────────────────
# Path to the mesh file written by gflex.export_for_blender()
MESH_FILE    = "/tmp/gflex_blender_mesh.py"

# Output PNG (RGBA, transparent background)
OUTPUT       = "/tmp/gflex_render.png"

# Render settings
RESOLUTION_X = 1600
RESOLUTION_Y = 900
SAMPLES      = 256        # Cycles samples; 64 for drafts, 512+ for final

# Scene options
SHOW_GRID    = True       # wireframe grid coloured by Te (requires te_verts)
SHOW_LOAD    = True       # load surface above the plate (requires qs_norm)
GRID_STEP    = 12         # subsample every N rows/cols for the wireframe

# Camera — adjust these to frame your scene
# (x forward, y left-right, z up in Blender world space)
CAM_OFFSET   = (0.0, -1.8, 0.75)    # centred in x; pulled back and up to frame scene
CAM_LENS     = 35                    # slightly wider to capture more vertical extent
CAM_SHIFT_Y  = 0.0                   # no film shift — target computed from scene bounds

# ── Load mesh data ─────────────────────────────────────────────────────────────
_ns = {}
exec(open(MESH_FILE).read(), _ns)

vertices      = _ns["vertices"]
faces         = _ns["faces"]
w_min_bu      = _ns["w_min_bu"]
w_max_bu      = _ns["w_max_bu"]
clip_bu       = _ns["clip_bu"]
vmax_bu       = _ns.get("vmax_bu", clip_bu)   # asymmetric positive clip for forebulge
nrows         = _ns["nrows"]
ncols         = _ns["ncols"]
plate_size_bu = _ns["plate_size_bu"]

# Optional attributes (None if not written by export_for_blender)
qs_norm          = _ns.get("qs_norm")
te_verts         = _ns.get("te_verts")
te_min           = _ns.get("te_min")
te_max           = _ns.get("te_max")
load_radius_bu   = _ns.get("load_radius_bu")
load_height_bu   = _ns.get("load_height_bu")
load_height_m_eq = _ns.get("load_height_m_eq")
base_z_bu        = _ns.get("base_z_bu", 0.0)

print(f"Mesh: {len(vertices)} verts, {len(faces)} quads")
print(f"Z range: {w_min_bu:.3f} … {w_max_bu:.3f} BU  |  clip: ±{clip_bu:.3f} BU")


# ── Scene helpers ──────────────────────────────────────────────────────────────

def clear_scene():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)
    for block in (bpy.data.meshes, bpy.data.lights, bpy.data.cameras,
                  bpy.data.materials, bpy.data.curves):
        for item in list(block):
            block.remove(item)


def _vik_ramp(cramp):
    """Apply Crameri vik colour ramp to a colour-ramp node."""
    cr = cramp.color_ramp
    cr.color_mode    = "RGB"
    cr.interpolation = "LINEAR"
    cr.elements[0].position = 0.00;  cr.elements[0].color = (0.10, 0.16, 0.51, 1)
    cr.elements[1].position = 1.00;  cr.elements[1].color = (0.52, 0.04, 0.06, 1)
    e1 = cr.elements.new(0.25);  e1.color = (0.24, 0.53, 0.80, 1)
    e2 = cr.elements.new(0.40);  e2.color = (0.53, 0.77, 0.92, 1)
    e3 = cr.elements.new(0.50);  e3.color = (0.97, 0.97, 0.97, 1)
    e4 = cr.elements.new(0.60);  e4.color = (0.91, 0.53, 0.37, 1)
    e5 = cr.elements.new(0.75);  e5.color = (0.82, 0.25, 0.16, 1)


# ── Deflection surface ─────────────────────────────────────────────────────────

def _vik_color(fac):
    """Return linear-interpolated (R, G, B) from the vik colour ramp."""
    stops = [
        (0.00, (0.10, 0.16, 0.51)),
        (0.25, (0.24, 0.53, 0.80)),
        (0.40, (0.53, 0.77, 0.92)),
        (0.50, (0.97, 0.97, 0.97)),
        (0.60, (0.91, 0.53, 0.37)),
        (0.75, (0.82, 0.25, 0.16)),
        (1.00, (0.52, 0.04, 0.06)),
    ]
    fac = max(0.0, min(1.0, fac))
    for k in range(len(stops) - 1):
        p0, c0 = stops[k]
        p1, c1 = stops[k + 1]
        if p0 <= fac <= p1:
            t = (fac - p0) / (p1 - p0)
            return tuple(c0[j] + t * (c1[j] - c0[j]) for j in range(3))
    return stops[-1][1]


def _tsn_fac(z):
    """TwoSlopeNorm: z → [0, 1] with z=0 → 0.5 and asymmetric arms."""
    if z <= 0.0:
        return (0.5 * (z - w_min_bu) / (-w_min_bu)) if w_min_bu < 0.0 else 0.5
    return (0.5 + 0.5 * z / vmax_bu) if vmax_bu > 0.0 else 0.5


def build_plate():
    """Deflection surface coloured with vik (TwoSlopeNorm, z=0 → white).

    Colours are pre-computed in Python and stored as a FLOAT_COLOR vertex
    attribute, which Cycles reads reliably via ShaderNodeAttribute → Color.
    This avoids issues with geometry-position node access in headless mode.
    """
    mesh = bpy.data.meshes.new("FlexureMesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for p in mesh.polygons:
        p.use_smooth = True

    # Pre-compute per-vertex vik colours with TwoSlopeNorm
    col_attr = mesh.color_attributes.new("w_col", "FLOAT_COLOR", "POINT")
    for k, v in enumerate(vertices):
        r, g, b = _vik_color(_tsn_fac(v[2]))
        col_attr.data[k].color = (r, g, b, 1.0)

    plate = bpy.data.objects.new("FlexurePlate", mesh)
    bpy.context.collection.objects.link(plate)

    mat = bpy.data.materials.new("PlateMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    vcol = nt.nodes.new("ShaderNodeVertexColor")
    vcol.layer_name = "w_col"
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    out  = nt.nodes.new("ShaderNodeOutputMaterial")

    for node, loc in zip([vcol, bsdf, out], [(-300, 0), (0, 0), (300, 0)]):
        node.location = loc

    nt.links.new(vcol.outputs["Color"], bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"],  out.inputs["Surface"])

    bsdf.inputs["Roughness"].default_value          = 0.55
    bsdf.inputs["Specular IOR Level"].default_value = 0.4
    plate.data.materials.append(mat)
    return plate


# ── Load cylinder (optional) ───────────────────────────────────────────────────

def build_load_cylinder():
    """Extruded cylinder representing the surface load in mantle-equivalent height.

    The cylinder base follows the deflected plate surface so the edifice sits
    naturally in the bowl it creates.  The top is offset uniformly by
    ``load_height_bu`` above each base vertex (same vertical exaggeration as
    the plate).  Its radius matches the inferred circular load footprint.
    """
    if load_radius_bu is None or load_height_bu is None:
        return None

    # Grid cell pitch in Blender units (from first two vertices in each direction)
    dx_bu = plate_size_bu / ncols
    dy_bu = vertices[ncols][1] - vertices[0][1]

    def plate_z_at(bx, by):
        """Nearest-neighbour deflection lookup at a Blender-unit (bx, by) position."""
        col = int(round((bx - vertices[0][0]) / dx_bu))
        row = int(round((by - vertices[0][1]) / dy_bu))
        col = max(0, min(ncols - 1, col))
        row = max(0, min(nrows - 1, row))
        return vertices[row * ncols + col][2]

    n_seg  = 64
    a_step = 2.0 * math.pi / n_seg
    r      = load_radius_bu

    # Bottom ring: base follows the deflected plate surface
    bot_ring_xy = [(r * math.cos(k * a_step), r * math.sin(k * a_step))
                   for k in range(n_seg)]
    bot_ring = [(x, y, plate_z_at(x, y)) for x, y in bot_ring_xy]

    # Top ring: uniformly above each bottom vertex by load_height_bu
    top_ring = [(x, y, z + load_height_bu) for x, y, z in bot_ring]

    bot_z_c = plate_z_at(0.0, 0.0)
    top_c   = (0.0, 0.0, bot_z_c + load_height_bu)
    bot_c   = (0.0, 0.0, bot_z_c)

    # Index layout: top_ring[0..n-1], top_c=n, bot_ring[n+1..2n], bot_c=2n+1
    i_tc  = n_seg
    i_br0 = n_seg + 1
    i_bc  = 2 * n_seg + 1
    all_verts = top_ring + [top_c] + bot_ring + [bot_c]

    all_faces = []
    for k in range(n_seg):
        k1 = (k + 1) % n_seg
        # Side quad — winding gives outward normal (validated analytically)
        all_faces.append((k1, k, i_br0 + k, i_br0 + k1))
        # Top cap — normal points +Z
        all_faces.append((i_tc, k, k1))
        # Bottom cap — normal points -Z
        all_faces.append((i_bc, i_br0 + k1, i_br0 + k))

    mesh = bpy.data.meshes.new("LoadCylMesh")
    mesh.from_pydata(all_verts, [], all_faces)
    mesh.validate()
    mesh.update()
    for p in mesh.polygons:
        p.use_smooth = True

    obj = bpy.data.objects.new("LoadCylinder", mesh)
    bpy.context.collection.objects.link(obj)

    mat = bpy.data.materials.new("LoadCylMat")
    mat.use_nodes = True
    pb = mat.node_tree.nodes["Principled BSDF"]

    # Dark basalt — clearly distinct from both the blue subsidence and vik forebulge reds
    pb.inputs["Base Color"].default_value         = (0.18, 0.12, 0.08, 1)
    pb.inputs["Roughness"].default_value          = 0.85
    pb.inputs["Specular IOR Level"].default_value = 0.05
    obj.data.materials.append(mat)
    return obj


# ── Plate body (solid slab) ────────────────────────────────────────────────────

def _attach_float_attr(mesh, name, values):
    attr = mesh.attributes.new(name, "FLOAT", "POINT")
    attr.data.foreach_set("value", values)


def build_te_slab():
    """Solid slab: bottom surface + perimeter side walls, exactly as in the logo.

    Each bottom vertex is offset from its corresponding top (deflection) vertex
    by a Te-proportional depth, so the bottom surface follows the same curvature
    as the deflection and the plate has physically variable thickness.  Falls back
    to a uniform thickness equal to 12 % of plate_size_bu when te_verts is None.
    """
    if te_verts is not None and te_max is not None:
        te_z_scale = 2.0 * abs(w_min_bu) / te_max
        depth = [t * te_z_scale for t in te_verts]
    else:
        uniform = plate_size_bu * 0.12
        depth = [uniform] * len(vertices)

    n_top     = len(vertices)
    bot_verts = [(v[0], v[1], v[2] - depth[k]) for k, v in enumerate(vertices)]
    all_verts = list(vertices) + bot_verts

    # Bottom faces — reversed winding so normals point downward
    bot_faces = [(f[3] + n_top, f[2] + n_top, f[1] + n_top, f[0] + n_top)
                 for f in faces]

    # Perimeter side faces
    side_faces = []
    for i in range(ncols - 1):           # south edge  (j = 0)
        t0, t1 = i, i + 1
        side_faces.append((t0, t1, t1 + n_top, t0 + n_top))
    for i in range(ncols - 1):           # north edge  (j = nrows-1)
        t0 = (nrows - 1) * ncols + i;  t1 = t0 + 1
        side_faces.append((t1, t0, t0 + n_top, t1 + n_top))
    for j in range(nrows - 1):           # west edge   (i = 0)
        t0 = j * ncols;  t1 = (j + 1) * ncols
        side_faces.append((t1, t0, t0 + n_top, t1 + n_top))
    for j in range(nrows - 1):           # east edge   (i = ncols-1)
        t0 = j * ncols + (ncols - 1);  t1 = (j + 1) * ncols + (ncols - 1)
        side_faces.append((t0, t1, t1 + n_top, t0 + n_top))

    mesh = bpy.data.meshes.new("TeSlabMesh")
    mesh.from_pydata(all_verts, [], bot_faces + side_faces)
    mesh.update()
    for p in mesh.polygons:
        p.use_smooth = True

    obj = bpy.data.objects.new("TeSlab", mesh)
    bpy.context.collection.objects.link(obj)

    mat = bpy.data.materials.new("TeSlabMat")
    mat.use_nodes = True
    pb = mat.node_tree.nodes["Principled BSDF"]
    pb.inputs["Base Color"].default_value         = (0.72, 0.70, 0.67, 1)
    pb.inputs["Roughness"].default_value          = 0.80
    pb.inputs["Specular IOR Level"].default_value = 0.10
    obj.data.materials.append(mat)
    return obj


# ── Te wireframe grid (optional) ───────────────────────────────────────────────

def build_te_grid(step=GRID_STEP):
    """Subsampled wireframe on the deflection surface, coloured by Te."""
    if te_verts is None:
        return None

    j_vals = list(range(0, nrows, step))
    if j_vals[-1] != nrows - 1:
        j_vals.append(nrows - 1)
    i_vals = list(range(0, ncols, step))
    if i_vals[-1] != ncols - 1:
        i_vals.append(ncols - 1)

    idx_map  = {}
    gverts   = []
    gte_vals = []
    for jj, j in enumerate(j_vals):
        for ii, i in enumerate(i_vals):
            idx_map[(j, i)] = len(gverts)
            gverts.append(vertices[j * ncols + i])
            gte_vals.append(te_verts[j * ncols + i])

    gfaces = []
    for jj in range(len(j_vals) - 1):
        for ii in range(len(i_vals) - 1):
            j0, j1 = j_vals[jj], j_vals[jj + 1]
            i0, i1 = i_vals[ii], i_vals[ii + 1]
            gfaces.append((idx_map[(j0, i0)], idx_map[(j0, i1)],
                           idx_map[(j1, i1)], idx_map[(j1, i0)]))

    mesh = bpy.data.meshes.new("TeGridMesh")
    mesh.from_pydata(gverts, [], gfaces)
    mesh.update()
    _attach_float_attr(mesh, "te_val", gte_vals)

    obj = bpy.data.objects.new("TeGrid", mesh)
    bpy.context.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type="WIREFRAME")
    obj.modifiers["Wireframe"].thickness       = 0.009
    obj.modifiers["Wireframe"].use_even_offset = True

    mat = bpy.data.materials.new("TeGridMat")
    mat.use_nodes = True
    pb = mat.node_tree.nodes["Principled BSDF"]
    pb.inputs["Base Color"].default_value         = (0.04, 0.04, 0.08, 1)
    pb.inputs["Roughness"].default_value          = 0.90
    pb.inputs["Specular IOR Level"].default_value = 0.05
    obj.data.materials.append(mat)
    return obj


# ── Camera ─────────────────────────────────────────────────────────────────────

def build_camera():
    cam_data = bpy.data.cameras.new("Camera")
    cam_data.lens    = CAM_LENS
    cam_data.shift_y = CAM_SHIFT_Y
    cam_obj = bpy.data.objects.new("Camera", cam_data)
    bpy.context.collection.objects.link(cam_obj)
    bpy.context.scene.camera = cam_obj

    s = plate_size_bu
    cam_loc = mathutils.Vector((
        s * CAM_OFFSET[0],
        s * CAM_OFFSET[1],
        s * CAM_OFFSET[2],
    ))
    cam_obj.location = cam_loc
    # Target the geometric centre of the scene (between slab bottom and cylinder top).
    # Cylinder top = deepest deflection + load height (base deforms with plate).
    # Slab bottom  ≈ 2 × |w_min| below z=0 (thickest-Te edge, w≈0).
    cyl_top_z  = (w_min_bu + load_height_bu) if load_height_bu is not None else 0.0
    slab_bot_z = -2.0 * abs(w_min_bu)
    target_z   = 0.5 * (cyl_top_z + slab_bot_z)
    target     = mathutils.Vector((0.0, 0.0, target_z))
    direction = target - cam_loc
    cam_obj.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
    return cam_obj


# ── Lighting ───────────────────────────────────────────────────────────────────

def build_lights():
    def add_sun(name, energy, color, rotation_euler):
        d = bpy.data.lights.new(name, "SUN")
        d.energy = energy
        d.color  = color
        o = bpy.data.objects.new(name, d)
        bpy.context.collection.objects.link(o)
        o.rotation_euler = rotation_euler
        return o

    add_sun("Key",  4.0, (1.00, 0.93, 0.80),
            (math.radians(55), math.radians(0),  math.radians(35)))
    add_sun("Fill", 1.2, (0.65, 0.80, 1.00),
            (math.radians(40), math.radians(0),  math.radians(-120)))
    add_sun("Rim",  0.6, (1.00, 0.95, 0.85),
            (math.radians(20), math.radians(0),  math.radians(160)))


# ── World background ───────────────────────────────────────────────────────────

def build_world():
    world = bpy.data.worlds["World"]
    world.use_nodes = True
    bg = world.node_tree.nodes["Background"]
    bg.inputs["Color"].default_value    = (0.01, 0.02, 0.06, 1)
    bg.inputs["Strength"].default_value = 0.25


# ── Render ─────────────────────────────────────────────────────────────────────

def render_to(path):
    scene = bpy.context.scene
    scene.render.engine                          = "CYCLES"
    scene.cycles.samples                         = SAMPLES
    scene.cycles.use_denoising                   = False
    scene.render.resolution_x                    = RESOLUTION_X
    scene.render.resolution_y                    = RESOLUTION_Y
    scene.render.film_transparent                = True
    scene.render.image_settings.file_format      = "PNG"
    scene.render.image_settings.color_mode       = "RGBA"
    scene.render.filepath                        = path
    bpy.ops.render.render(write_still=True)
    print(f"✓  Rendered → {path}")


# ── Build and render ───────────────────────────────────────────────────────────

clear_scene()
build_plate()
build_te_slab()
build_camera()
build_lights()
build_world()

if SHOW_LOAD:
    build_load_cylinder()
if SHOW_GRID:
    build_te_grid()

render_to(OUTPUT)
