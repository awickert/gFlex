"""Blender 4 Python script — renders the gFlex logo.

Run headless:
    blender --background --python scripts/blender_logo.py

Reads /tmp/gflex_logo_mesh.py (produced by generate_logo_data.py) and
writes docs/_static/logo.png (RGBA, 1200×900).
"""

import math
import os

import bpy
import mathutils

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MESH   = "/tmp/gflex_logo_mesh.py"
OUTDIR = os.path.join(REPO, "docs", "_static")

# ── Load pre-computed mesh (pure Python, no numpy needed) ─────────────────────
ns = {}
exec(open(MESH).read(), ns)

vertices    = ns["vertices"]
faces       = ns["faces"]
te_verts    = ns["te_verts"]
w_min_bu    = ns["w_min_bu"]
w_max_bu    = ns["w_max_bu"]
load_bx     = ns["load_bx"]
load_by     = ns["load_by"]
load_br     = ns["load_br"]
base_z      = ns["base_z"]
load_height = ns["load_height"]
clip_bu     = ns["clip_bu"]
nrows       = ns["nrows"]
ncols       = ns["ncols"]
te_min      = ns["te_min"]
te_max      = ns["te_max"]

print(f"Mesh: {len(vertices)} verts, {len(faces)} quads")
print(f"Z range: {w_min_bu:.3f} … {w_max_bu:.3f} BU")
print(f"Te range: {te_min/1e3:.0f} … {te_max/1e3:.0f} km")


# ── Scene helpers ──────────────────────────────────────────────────────────────

def clear_scene():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)
    for block in (bpy.data.meshes, bpy.data.lights, bpy.data.cameras,
                  bpy.data.materials, bpy.data.curves):
        for item in list(block):
            block.remove(item)


def build_plate():
    """Build the flexure plate mesh with deflection-coloured material."""
    mesh = bpy.data.meshes.new("FlexureMesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for p in mesh.polygons:
        p.use_smooth = True

    plate = bpy.data.objects.new("FlexurePlate", mesh)
    bpy.context.collection.objects.link(plate)

    mat = bpy.data.materials.new("Plate")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    geom   = nt.nodes.new("ShaderNodeNewGeometry")
    sep    = nt.nodes.new("ShaderNodeSeparateXYZ")
    mrange = nt.nodes.new("ShaderNodeMapRange")
    cramp  = nt.nodes.new("ShaderNodeValToRGB")
    bsdf   = nt.nodes.new("ShaderNodeBsdfPrincipled")
    out    = nt.nodes.new("ShaderNodeOutputMaterial")

    for n, loc in zip([geom, sep, mrange, cramp, bsdf, out],
                      [(-900,0),(-700,0),(-500,0),(-300,0),(0,0),(300,0)]):
        n.location = loc

    nt.links.new(geom.outputs["Position"],  sep.inputs["Vector"])
    nt.links.new(sep.outputs["Z"],          mrange.inputs["Value"])
    nt.links.new(mrange.outputs["Result"],  cramp.inputs["Fac"])
    nt.links.new(cramp.outputs["Color"],    bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"],      out.inputs["Surface"])

    clip_display = abs(w_max_bu) * 1.4
    mrange.inputs["From Min"].default_value = -clip_display
    mrange.inputs["From Max"].default_value =  clip_display
    mrange.inputs["To Min"].default_value   = 0.0
    mrange.inputs["To Max"].default_value   = 1.0
    mrange.clamp = True

    cr = cramp.color_ramp
    cr.color_mode     = "RGB"
    cr.interpolation  = "LINEAR"
    cr.elements[0].position = 0.00;  cr.elements[0].color = (0.02, 0.07, 0.45, 1)
    cr.elements[1].position = 1.00;  cr.elements[1].color = (0.55, 0.02, 0.06, 1)
    e1 = cr.elements.new(0.20);  e1.color = (0.08, 0.25, 0.75, 1)
    e2 = cr.elements.new(0.38);  e2.color = (0.45, 0.65, 0.90, 1)
    e3 = cr.elements.new(0.50);  e3.color = (0.96, 0.96, 0.96, 1)
    e4 = cr.elements.new(0.62);  e4.color = (0.90, 0.45, 0.30, 1)
    e5 = cr.elements.new(0.80);  e5.color = (0.80, 0.10, 0.08, 1)

    bsdf.inputs["Roughness"].default_value          = 0.55
    bsdf.inputs["Specular IOR Level"].default_value = 0.4
    plate.data.materials.append(mat)
    return plate


def build_cylinder():
    """Build the load cylinder."""
    height = abs(w_min_bu) * 1.0
    bpy.ops.mesh.primitive_cylinder_add(
        radius=load_br,
        depth=height,
        location=(load_bx, load_by, base_z + height / 2.0),
        vertices=64,
    )
    cyl = bpy.context.active_object
    cyl.name = "Load"
    for p in cyl.data.polygons:
        p.use_smooth = True

    cyl_mat = bpy.data.materials.new("Load")
    cyl_mat.use_nodes = True
    pb = cyl_mat.node_tree.nodes["Principled BSDF"]
    pb.inputs["Base Color"].default_value          = (0.06, 0.05, 0.04, 1)
    pb.inputs["Roughness"].default_value           = 0.80
    pb.inputs["Specular IOR Level"].default_value  = 0.05
    cyl.data.materials.append(cyl_mat)
    return cyl


def build_camera():
    cam_data = bpy.data.cameras.new("Camera")
    cam_data.lens    = 40
    cam_data.shift_y = -0.13   # push plate up, trimming dead space above
    cam_obj = bpy.data.objects.new("Camera", cam_data)
    bpy.context.collection.objects.link(cam_obj)
    bpy.context.scene.camera = cam_obj
    cam_loc = mathutils.Vector((2.0, -13.0, 4.5))
    cam_obj.location = cam_loc
    target    = mathutils.Vector((0.0, 0.5, w_min_bu * 0.25))
    direction = target - cam_loc
    cam_obj.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
    return cam_obj


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


def build_world():
    world = bpy.data.worlds["World"]
    world.use_nodes = True
    bg = world.node_tree.nodes["Background"]
    bg.inputs["Color"].default_value    = (0.01, 0.02, 0.06, 1)
    bg.inputs["Strength"].default_value = 0.25


def render_to(filename):
    scene = bpy.context.scene
    scene.render.engine            = "CYCLES"
    scene.cycles.samples           = 256
    scene.cycles.use_denoising     = False
    scene.render.resolution_x      = 1200
    scene.render.resolution_y      = 900
    scene.render.film_transparent  = True
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGBA"
    scene.render.filepath = os.path.join(OUTDIR, filename)
    bpy.ops.render.render(write_still=True)
    print(f"✓  Rendered → {scene.render.filepath}")


def remove_obj(obj):
    if obj is None:
        return
    mesh = obj.data if hasattr(obj.data, "name") else None
    bpy.data.objects.remove(obj, do_unlink=True)
    if mesh and mesh.users == 0:
        bpy.data.meshes.remove(mesh)


# ── Te-display variant builders ────────────────────────────────────────────────

def _te_material(name):
    """Grayscale material driven by a 'te_val' FLOAT vertex attribute."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    attr   = nt.nodes.new("ShaderNodeAttribute");  attr.attribute_name = "te_val"
    mrange = nt.nodes.new("ShaderNodeMapRange")
    cramp  = nt.nodes.new("ShaderNodeValToRGB")
    bsdf   = nt.nodes.new("ShaderNodeBsdfPrincipled")
    out    = nt.nodes.new("ShaderNodeOutputMaterial")

    for n, loc in zip([attr, mrange, cramp, bsdf, out],
                      [(-700,0),(-500,0),(-300,0),(0,0),(300,0)]):
        n.location = loc

    mrange.inputs["From Min"].default_value = te_min
    mrange.inputs["From Max"].default_value = te_max
    mrange.inputs["To Min"].default_value   = 0.0
    mrange.inputs["To Max"].default_value   = 1.0
    mrange.clamp = True

    cr = cramp.color_ramp
    cr.color_mode    = "RGB"
    cr.interpolation = "LINEAR"
    cr.elements[0].position = 0.0;  cr.elements[0].color = (0.85, 0.85, 0.85, 1)  # light (thin)
    cr.elements[1].position = 1.0;  cr.elements[1].color = (0.20, 0.20, 0.20, 1)  # dark  (thick)

    bsdf.inputs["Roughness"].default_value          = 0.85
    bsdf.inputs["Specular IOR Level"].default_value = 0.1

    nt.links.new(attr.outputs["Fac"],      mrange.inputs["Value"])
    nt.links.new(mrange.outputs["Result"], cramp.inputs["Fac"])
    nt.links.new(cramp.outputs["Color"],   bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"],     out.inputs["Surface"])
    return mat


def _te_slab_material():
    """
    Slab material: Te-based warm gray (thin=light, thick=darker) multiplied
    by a Z-position depth gradient (surface=white/no-effect, deep=dark ochre).
    The multiply keeps Te variation readable at every depth while the face
    darkens and warms toward the base of the lithosphere.
    """
    mat = bpy.data.materials.new("TeSlabMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    # ── Te horizontal variation ────────────────────────────────────────────────
    attr_te  = nt.nodes.new("ShaderNodeAttribute");  attr_te.attribute_name = "te_val"
    te_rng   = nt.nodes.new("ShaderNodeMapRange")
    te_cramp = nt.nodes.new("ShaderNodeValToRGB")

    te_rng.inputs["From Min"].default_value = te_min
    te_rng.inputs["From Max"].default_value = te_max
    te_rng.inputs["To Min"].default_value   = 0.0
    te_rng.inputs["To Max"].default_value   = 1.0
    te_rng.clamp = True

    cr = te_cramp.color_ramp
    cr.color_mode    = "RGB"
    cr.interpolation = "LINEAR"
    cr.elements[0].position = 0.0;  cr.elements[0].color = (1.00, 0.97, 0.93, 1)  # thin: pale warm
    cr.elements[1].position = 1.0;  cr.elements[1].color = (0.50, 0.48, 0.45, 1)  # thick: medium warm gray

    # ── Depth gradient from Z position ────────────────────────────────────────
    geom    = nt.nodes.new("ShaderNodeNewGeometry")
    sep     = nt.nodes.new("ShaderNodeSeparateXYZ")
    z_rng   = nt.nodes.new("ShaderNodeMapRange")
    z_cramp = nt.nodes.new("ShaderNodeValToRGB")

    # Z = 0 at plate surface, w_min_bu*3 at deepest slab point
    z_rng.inputs["From Min"].default_value = 0.0
    z_rng.inputs["From Max"].default_value = 3.0 * w_min_bu   # negative
    z_rng.inputs["To Min"].default_value   = 0.0
    z_rng.inputs["To Max"].default_value   = 1.0
    z_rng.clamp = True

    zc = z_cramp.color_ramp
    zc.color_mode    = "RGB"
    zc.interpolation = "LINEAR"
    zc.elements[0].position = 0.00;  zc.elements[0].color = (1.00, 1.00, 1.00, 1)  # surface: white (no tint)
    zc.elements[1].position = 1.00;  zc.elements[1].color = (0.14, 0.07, 0.03, 1)  # deep: dark ochre
    e1 = zc.elements.new(0.30);  e1.color = (0.88, 0.74, 0.56, 1)   # upper mid: warm sandy
    e2 = zc.elements.new(0.65);  e2.color = (0.42, 0.22, 0.10, 1)   # lower mid: terracotta

    # ── Multiply: depth darkens + warms the Te gray ────────────────────────────
    mul  = nt.nodes.new("ShaderNodeMixRGB")
    mul.blend_type = "MULTIPLY"
    mul.inputs["Fac"].default_value = 1.0

    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    out  = nt.nodes.new("ShaderNodeOutputMaterial")

    bsdf.inputs["Roughness"].default_value          = 0.80
    bsdf.inputs["Specular IOR Level"].default_value = 0.08

    nt.links.new(attr_te.outputs["Fac"],     te_rng.inputs["Value"])
    nt.links.new(te_rng.outputs["Result"],   te_cramp.inputs["Fac"])
    nt.links.new(geom.outputs["Position"],   sep.inputs["Vector"])
    nt.links.new(sep.outputs["Z"],           z_rng.inputs["Value"])
    nt.links.new(z_rng.outputs["Result"],    z_cramp.inputs["Fac"])
    nt.links.new(te_cramp.outputs["Color"],  mul.inputs["Color1"])
    nt.links.new(z_cramp.outputs["Color"],   mul.inputs["Color2"])
    nt.links.new(mul.outputs["Color"],       bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"],       out.inputs["Surface"])
    return mat


def _attach_te_attribute(mesh, te_values):
    attr = mesh.attributes.new("te_val", "FLOAT", "POINT")
    attr.data.foreach_set("value", te_values)


def build_te_floor():
    """Flat plane sitting below the plate, grey-shaded by Te."""
    floor_z = w_min_bu * 1.8
    floor_verts = [(v[0], v[1], floor_z) for v in vertices]

    mesh = bpy.data.meshes.new("TeFloorMesh")
    mesh.from_pydata(floor_verts, [], faces)
    mesh.update()
    for p in mesh.polygons:
        p.use_smooth = False

    _attach_te_attribute(mesh, te_verts)

    obj = bpy.data.objects.new("TeFloor", mesh)
    bpy.context.collection.objects.link(obj)
    obj.data.materials.append(_te_material("TeFloorMat"))
    return obj


def build_te_grid(step=12):
    """Subsampled wireframe sitting on the deflection surface."""
    j_vals = list(range(0, nrows, step))
    i_vals = list(range(0, ncols, step))

    idx_map = {}
    gverts  = []
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

    obj = bpy.data.objects.new("TeGrid", mesh)
    bpy.context.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='WIREFRAME')
    obj.modifiers['Wireframe'].thickness        = 0.009
    obj.modifiers['Wireframe'].use_even_offset  = True

    mat = bpy.data.materials.new("TeGridMat")
    mat.use_nodes = True
    pb = mat.node_tree.nodes["Principled BSDF"]
    pb.inputs["Base Color"].default_value         = (0.04, 0.04, 0.08, 1)  # near-black navy
    pb.inputs["Roughness"].default_value          = 0.90
    pb.inputs["Specular IOR Level"].default_value = 0.05
    obj.data.materials.append(mat)
    return obj


def build_te_slab():
    """
    Solid slab: bottom surface + perimeter side walls.
    Uniform gray material — thickness variation reads from the geometry and
    lighting alone, with no redundant color coding.
    """
    te_z_scale = 2.0 * abs(w_min_bu) / te_max

    n_top     = len(vertices)
    bot_verts = [(v[0], v[1], v[2] - te_verts[k] * te_z_scale)
                 for k, v in enumerate(vertices)]
    all_verts = list(vertices) + bot_verts

    # Bottom surface (reversed winding → normals point down)
    bot_faces = [(f[3] + n_top, f[2] + n_top, f[1] + n_top, f[0] + n_top)
                 for f in faces]

    # Perimeter side faces
    side_faces = []
    for i in range(ncols - 1):          # South edge (j = 0)
        t0, t1 = i, i + 1
        side_faces.append((t0, t1, t1 + n_top, t0 + n_top))
    for i in range(ncols - 1):          # North edge (j = nrows-1)
        t0 = (nrows - 1) * ncols + i;  t1 = t0 + 1
        side_faces.append((t1, t0, t0 + n_top, t1 + n_top))
    for j in range(nrows - 1):          # West edge (i = 0)
        t0 = j * ncols;  t1 = (j + 1) * ncols
        side_faces.append((t1, t0, t0 + n_top, t1 + n_top))
    for j in range(nrows - 1):          # East edge (i = ncols-1)
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
    pb.inputs["Base Color"].default_value         = (0.72, 0.70, 0.67, 1)  # warm mid-gray
    pb.inputs["Roughness"].default_value          = 0.80
    pb.inputs["Specular IOR Level"].default_value = 0.10
    obj.data.materials.append(mat)
    return obj


def build_base_plane():
    """
    Horizontal plane just below the deepest slab point, representing the
    base of the lithosphere / top of the asthenosphere.  The slab casts a
    shadow onto it; the shadow contrast gives the sense of depth.
    """
    z_plane = 3.4 * w_min_bu        # slightly below deepest slab (~3 × w_min_bu)
    half    = 8.0                   # wider than the plate so it peeks out at edges

    plane_verts = [
        (-half, -half, z_plane),
        ( half, -half, z_plane),
        ( half,  half, z_plane),
        (-half,  half, z_plane),
    ]
    plane_faces = [(0, 1, 2, 3)]

    mesh = bpy.data.meshes.new("BasePlaneMesh")
    mesh.from_pydata(plane_verts, [], plane_faces)
    mesh.update()

    obj = bpy.data.objects.new("BasePlane", mesh)
    bpy.context.collection.objects.link(obj)

    mat = bpy.data.materials.new("BasePlaneMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    # Principled BSDF (receives shadows) mixed with faint warm emission (mantle heat)
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    emit = nt.nodes.new("ShaderNodeEmission")
    mix  = nt.nodes.new("ShaderNodeMixShader")
    out  = nt.nodes.new("ShaderNodeOutputMaterial")

    bsdf.inputs["Base Color"].default_value         = (0.18, 0.09, 0.04, 1)  # dark warm rock
    bsdf.inputs["Roughness"].default_value          = 0.95
    bsdf.inputs["Specular IOR Level"].default_value = 0.0

    emit.inputs["Color"].default_value    = (1.0, 0.35, 0.08, 1)   # hot orange glow
    emit.inputs["Strength"].default_value = 0.30                    # subtle

    mix.inputs["Fac"].default_value = 0.15   # mostly BSDF (shadow-receiving), hint of glow

    nt.links.new(bsdf.outputs["BSDF"],     mix.inputs[1])
    nt.links.new(emit.outputs["Emission"], mix.inputs[2])
    nt.links.new(mix.outputs["Shader"],    out.inputs["Surface"])

    obj.data.materials.append(mat)
    return obj


# ── Build the static part of the scene ────────────────────────────────────────

clear_scene()
build_plate()
build_cylinder()
build_camera()
build_lights()
build_world()

# ── Render three variants ──────────────────────────────────────────────────────

build_te_slab()
build_te_grid()
render_to("logo.png")
