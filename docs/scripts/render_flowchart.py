#!/usr/bin/env python3
"""Render docs/scripts/flowchart.graphml to docs/_static/flowchart.svg (and PNG via inkscape)."""
import os
import subprocess
import xml.etree.ElementTree as ET
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DOCS_DIR   = os.path.dirname(SCRIPT_DIR)
GRAPHML    = os.path.join(SCRIPT_DIR, 'flowchart.graphml')
SVG_OUT    = os.path.join(DOCS_DIR, '_static', 'flowchart_new.svg')
PNG_OUT    = os.path.join(DOCS_DIR, '_static', 'flowchart_new.png')

GNS  = 'http://graphml.graphdrawing.org/graphml'
YNS  = 'http://www.yworks.com/xml/graphml'
DPI  = 150
ALEN = 11   # arrowhead length px

# ── Parse GraphML ─────────────────────────────────────────────────────────────

def parse():
    tree = ET.parse(GRAPHML)
    root = tree.getroot()

    nodes = {}
    for el in root.findall(f'.//{{{GNS}}}node'):
        nid = el.get('id')
        sn  = el.find(f'.//{{{YNS}}}ShapeNode')
        if sn is None:
            continue
        geom  = sn.find(f'{{{YNS}}}Geometry')
        fill  = sn.find(f'{{{YNS}}}Fill')
        bdr   = sn.find(f'{{{YNS}}}BorderStyle')
        lbl   = sn.find(f'{{{YNS}}}NodeLabel')
        shp   = sn.find(f'{{{YNS}}}Shape')
        style = lbl.get('fontStyle', 'plain') if lbl is not None else 'plain'
        nodes[nid] = dict(
            x          = float(geom.get('x')),
            y          = float(geom.get('y')),
            w          = float(geom.get('width')),
            h          = float(geom.get('height')),
            fill       = fill.get('color', '#E8E8E8') if fill is not None else '#E8E8E8',
            lw         = float(bdr.get('width', 1.5))  if bdr  is not None else 1.5,
            ec         = bdr.get('color', '#000000')    if bdr  is not None else '#000000',
            label      = lbl.text or ''                 if lbl  is not None else '',
            fs         = int(lbl.get('fontSize', 12))   if lbl  is not None else 12,
            bold       = ('bold'   in style),
            italic     = ('italic' in style),
            text_color = lbl.get('textColor', '#000000') if lbl is not None else '#000000',
            shape      = shp.get('type', 'rectangle')   if shp  is not None else 'rectangle',
            background = nid.startswith('bg_'),
        )

    edges = []
    for el in root.findall(f'.//{{{GNS}}}edge'):
        poly   = el.find(f'.//{{{YNS}}}PolyLineEdge')
        ls     = poly.find(f'{{{YNS}}}LineStyle') if poly is not None else None
        dashed = (ls.get('type') == 'dashed') if ls is not None else False
        bends = [(float(pt.get('x')), float(pt.get('y')))
                 for pt in (poly.findall(f'{{{YNS}}}Point') if poly is not None else [])]
        edges.append(dict(src=el.get('source'), dst=el.get('target'), dashed=dashed, bends=bends))

    return nodes, edges

# ── SVG helpers ───────────────────────────────────────────────────────────────

def q(v):
    return f'{v:.2f}'


def svg_node(nd, background=False):
    x, y, w, h = nd['x'], nd['y'], nd['w'], nd['h']
    cx, cy = x + w/2, y + h/2
    fill = nd['fill'];  ec = nd['ec'];  sw = nd['lw']
    lines = []

    if nd['shape'] == 'ellipse':
        lines.append(f'  <ellipse cx="{q(cx)}" cy="{q(cy)}" rx="{q(w/2)}" ry="{q(h/2)}"'
                     f' fill="{fill}" stroke="{ec}" stroke-width="{sw}"/>')
    elif nd['shape'] == 'roundrectangle':
        r = min(10, h * 0.14)
        lines.append(f'  <rect x="{q(x)}" y="{q(y)}" width="{w}" height="{h}"'
                     f' rx="{q(r)}" ry="{q(r)}" fill="{fill}" stroke="{ec}" stroke-width="{sw}"/>')
    elif nd['shape'] == 'diamond':
        pts = f'{q(cx)},{q(y)} {q(x+w)},{q(cy)} {q(cx)},{q(y+h)} {q(x)},{q(cy)}'
        lines.append(f'  <polygon points="{pts}" fill="{fill}" stroke="{ec}" stroke-width="{sw}"/>')
    elif nd['shape'] == 'hexagon':
        d = h * 0.3
        pts = (f'{q(x+d)},{q(y)} {q(x+w-d)},{q(y)} {q(x+w)},{q(cy)} '
               f'{q(x+w-d)},{q(y+h)} {q(x+d)},{q(y+h)} {q(x)},{q(cy)}')
        lines.append(f'  <polygon points="{pts}" fill="{fill}" stroke="{ec}" stroke-width="{sw}"/>')

    label_lines = nd['label'].split('\n')
    nl = len(label_lines)
    fs = nd['fs']
    lh = fs * 1.35
    fw = 'bold' if nd['bold'] else 'normal'
    tc = nd.get('text_color', '#000000')
    fi = 'italic' if nd.get('italic') else 'normal'

    if background:
        # Rotated swim-lane label: centered vertically on the 20 px left strip
        sx = x + 10          # horizontal center of the label strip
        sy = y + h / 2       # vertical center of the band
        lines.append(f'  <text text-anchor="middle" dominant-baseline="central"'
                     f' font-family="Arial,Helvetica,sans-serif"'
                     f' font-size="{fs}" font-weight="{fw}" font-style="{fi}" fill="{tc}"'
                     f' transform="rotate(-90,{q(sx)},{q(sy)})">')
        for i, txt in enumerate(label_lines):
            ty = sy + (i - (nl - 1) / 2) * lh
            lines.append(f'    <tspan x="{q(sx)}" y="{q(ty)}" dominant-baseline="central">{txt}</tspan>')
        lines.append('  </text>')
    else:
        lines.append(f'  <text text-anchor="middle" font-family="Arial,Helvetica,sans-serif"'
                     f' font-size="{fs}" font-weight="{fw}">')
        for i, txt in enumerate(label_lines):
            ty = cy - (nl - 1) * lh / 2 + i * lh
            lines.append(f'    <tspan x="{q(cx)}" y="{q(ty)}" dominant-baseline="central">{txt}</tspan>')
        lines.append('  </text>')
    return lines


def svg_edge(nodes, src_id, dst_id, dashed, bends=None):
    sn = nodes[src_id];  dn = nodes[dst_id]
    x0 = sn['x'] + sn['w']/2;  y0 = sn['y'] + sn['h']
    x1 = dn['x'] + dn['w']/2;  y1 = dn['y']
    da = ' stroke-dasharray="13,8"' if dashed else ''

    if bends:
        lx, ly = bends[-1]
        ex, ey = x1 - lx, y1 - ly
        el = np.hypot(ex, ey)
        ux, uy = ex / el, ey / el
        xe, ye = x1 - ux * ALEN, y1 - uy * ALEN
        pts = ' '.join(f'{q(px)},{q(py)}' for px, py in [(x0, y0)] + list(bends) + [(xe, ye)])
        return (f'  <polyline points="{pts}" fill="none"'
                f' stroke="#000000" stroke-width="1.5"{da} marker-end="url(#arr)"/>')

    ex, ey = x1 - x0, y1 - y0
    el = np.hypot(ex, ey)
    ux, uy = ex / el, ey / el
    xe, ye = x1 - ux * ALEN, y1 - uy * ALEN
    return (f'  <line x1="{q(x0)}" y1="{q(y0)}" x2="{q(xe)}" y2="{q(ye)}"'
            f' stroke="#000000" stroke-width="1.5"{da} marker-end="url(#arr)"/>')


def main():
    nodes, edges = parse()

    max_x = max(nd['x'] + nd['w'] for nd in nodes.values())
    max_y = max(nd['y'] + nd['h'] for nd in nodes.values())
    # infer margin from the minimum node origin (already baked in by Graphviz layout)
    margin = min(min(nd['x'] for nd in nodes.values()),
                 min(nd['y'] for nd in nodes.values()))
    W = int(max_x + margin)
    H = int(max_y + margin)

    o = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '  <defs>',
        '    <marker id="arr" markerWidth="11" markerHeight="8" refX="0" refY="4"',
        '            orient="auto" markerUnits="userSpaceOnUse">',
        '      <polygon points="0,0 11,4 0,8" fill="#000000"/>',
        '    </marker>',
        '  </defs>',
        f'  <rect width="{W}" height="{H}" fill="white"/>',
    ]

    for nd in nodes.values():               # background bands first
        if nd['background']:
            o.extend(svg_node(nd, background=True))
    for e in edges:
        o.append(svg_edge(nodes, e['src'], e['dst'], e['dashed'], e['bends'] or None))
    for nd in nodes.values():               # regular nodes on top
        if not nd['background']:
            o.extend(svg_node(nd))

    o.append('</svg>')

    with open(SVG_OUT, 'w') as f:
        f.write('\n'.join(o))
    print(f'Wrote {SVG_OUT}  ({W}×{H} px)')

    r = subprocess.run(['inkscape', '--export-type=png', f'--export-dpi={DPI}',
                        f'--export-filename={PNG_OUT}', SVG_OUT], capture_output=True)
    if r.returncode == 0:
        print(f'Wrote {PNG_OUT}')
    else:
        print('inkscape error:', r.stderr.decode()[:300])


if __name__ == '__main__':
    main()
