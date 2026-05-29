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
SHOW_TE_FLOOR = True      # flat plane below the plate coloured by Te
SHOW_LOAD    = True       # load surface above the plate (requires qs_norm)
GRID_STEP    = 12         # subsample every N rows/cols for the wireframe

# Camera — adjust these to frame your scene
# (x forward, y left-right, z up in Blender world space)
CAM_OFFSET   = (0.15, -1.4, 0.40)   # as fractions of plate_size_bu
CAM_LENS     = 40                    # focal length [mm]; larger = more telephoto
CAM_SHIFT_Y  = -0.07                 # sensor shift; negative pushes plate up

# ── Load mesh data ─────────────────────────────────────────────────────────────
_ns = {}
exec(open(MESH_FILE).read(), _ns)

vertices      = _ns["vertices"]
faces         = _ns["faces"]
w_min_bu      = _ns["w_min_bu"]
w_max_bu      = _ns["w_max_bu"]
clip_bu       = _ns["clip_bu"]
nrows         = _ns["nrows"]
ncols         = _ns["ncols"]
plate_size_bu = _ns["plate_size_bu"]

# Optional attributes (None if not written by export_for_blender)
qs_norm  = _ns.get("qs_norm")
te_verts = _ns.get("te_verts")
te_min   = _ns.get("te_min")
te_max   = _ns.get("te_max")

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

def build_plate():
    """Deflection surface coloured by Z position via the vik ramp."""
    mesh = bpy.data.meshes.new("FlexureMesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for p in mesh.polygons:
        p.use_smooth = True

    plate = bpy.data.objects.new("FlexurePlate", mesh)
    bpy.context.collection.objects.link(plate)

    mat = bpy.data.materials.new("PlateMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    geom   = nt.nodes.new("ShaderNodeNewGeometry")
    sep    = nt.nodes.new("ShaderNodeSeparateXYZ")
    mrange = nt.nodes.new("ShaderNodeMapRange")
    cramp  = nt.nodes.new("ShaderNodeValToRGB")
    bsdf   = nt.nodes.new("ShaderNodeBsdfPrincipled")
    out    = nt.nodes.new("ShaderNodeOutputMaterial")

    for node, loc in zip([geom, sep, mrange, cramp, bsdf, out],
                         [(-900,0), (-700,0), (-500,0), (-300,0), (0,0), (300,0)]):
        node.location = loc

    nt.links.new(geom.outputs["Position"], sep.inputs["Vector"])
    nt.links.new(sep.outputs["Z"],         mrange.inputs["Value"])
    nt.links.new(mrange.outputs["Result"], cramp.inputs["Fac"])
    nt.links.new(cramp.outputs["Color"],   bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"],     out.inputs["Surface"])

    # Map [-clip_bu, +clip_bu] → [0, 1] for the vik ramp
    mrange.inputs["From Min"].default_value = -clip_bu
    mrange.inputs["From Max"].default_value =  clip_bu
    mrange.inputs["To Min"].default_value   = 0.0
    mrange.inputs["To Max"].default_value   = 1.0
    mrange.clamp = True

    _vik_ramp(cramp)

    bsdf.inputs["Roughness"].default_value          = 0.55
    bsdf.inputs["Specular IOR Level"].default_value = 0.4
    plate.data.materials.append(mat)
    return plate


# ── Load surface (optional) ────────────────────────────────────────────────────

def build_load_surface():
    """Flat surface slightly above the plate, opacity driven by qs_norm."""
    if qs_norm is None:
        return None

    # Flat surface at the plate's zero level, same xy footprint
    load_verts = [(v[0], v[1], 0.0) for v in vertices]
    mesh = bpy.data.meshes.new("LoadMesh")
    mesh.from_pydata(load_verts, [], faces)
    mesh.update()

    # Per-vertex qs_norm as a colour attribute
    col_attr = mesh.color_attributes.new("qs_col", "FLOAT_COLOR", "POINT")
    for i, val in enumerate(qs_norm):
        col_attr.data[i].color = (val, val, val, 1.0)

    obj = bpy.data.objects.new("LoadSurface", mesh)
    bpy.context.collection.objects.link(obj)

    mat = bpy.data.materials.new("LoadMat")
    mat.use_nodes = True
    mat.blend_method = "BLEND"
    nt = mat.node_tree
    nt.nodes.clear()

    attr   = nt.nodes.new("ShaderNodeAttribute")
    attr.attribute_name = "qs_col"
    attr.attribute_type = "GEOMETRY"
    bsdf   = nt.nodes.new("ShaderNodeBsdfPrincipled")
    out    = nt.nodes.new("ShaderNodeOutputMaterial")

    nt.links.new(attr.outputs["Fac"],  bsdf.inputs["Alpha"])
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    bsdf.inputs["Base Color"].default_value         = (0.06, 0.05, 0.04, 1)
    bsdf.inputs["Roughness"].default_value          = 0.80
    bsdf.inputs["Specular IOR Level"].default_value = 0.05

    obj.data.materials.append(mat)
    return obj


# ── Te floor (optional) ────────────────────────────────────────────────────────

def _attach_float_attr(mesh, name, values):
    attr = mesh.attributes.new(name, "FLOAT", "POINT")
    attr.data.foreach_set("value", values)


def build_te_floor():
    """Flat plane below the plate, grey-shaded by elastic thickness."""
    if te_verts is None:
        return None

    floor_z     = w_min_bu * 1.8
    floor_verts = [(v[0], v[1], floor_z) for v in vertices]

    mesh = bpy.data.meshes.new("TeFloorMesh")
    mesh.from_pydata(floor_verts, [], faces)
    mesh.update()
    _attach_float_attr(mesh, "te_val", te_verts)

    obj = bpy.data.objects.new("TeFloor", mesh)
    bpy.context.collection.objects.link(obj)

    mat = bpy.data.materials.new("TeFloorMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    attr   = nt.nodes.new("ShaderNodeAttribute");  attr.attribute_name = "te_val"
    mrange = nt.nodes.new("ShaderNodeMapRange")
    cramp  = nt.nodes.new("ShaderNodeValToRGB")
    bsdf   = nt.nodes.new("ShaderNodeBsdfPrincipled")
    out    = nt.nodes.new("ShaderNodeOutputMaterial")

    mrange.inputs["From Min"].default_value = te_min
    mrange.inputs["From Max"].default_value = te_max
    mrange.inputs["To Min"].default_value   = 0.0
    mrange.inputs["To Max"].default_value   = 1.0
    mrange.clamp = True

    cr = cramp.color_ramp
    cr.color_mode    = "RGB"
    cr.interpolation = "LINEAR"
    cr.elements[0].position = 0.0;  cr.elements[0].color = (0.85, 0.85, 0.85, 1)
    cr.elements[1].position = 1.0;  cr.elements[1].color = (0.20, 0.20, 0.20, 1)

    nt.links.new(attr.outputs["Fac"],      mrange.inputs["Value"])
    nt.links.new(mrange.outputs["Result"], cramp.inputs["Fac"])
    nt.links.new(cramp.outputs["Color"],   bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"],     out.inputs["Surface"])

    bsdf.inputs["Roughness"].default_value          = 0.85
    bsdf.inputs["Specular IOR Level"].default_value = 0.1

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
    target    = mathutils.Vector((0.0, 0.0, w_min_bu * 0.25))
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
build_camera()
build_lights()
build_world()

if SHOW_LOAD:
    build_load_surface()
if SHOW_TE_FLOOR:
    build_te_floor()
if SHOW_GRID:
    build_te_grid()

render_to(OUTPUT)
