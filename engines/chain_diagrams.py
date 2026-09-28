"""Multi-vertex (chain) diagrams: the order-N amplitude drawn with N vertices.

Each vertex is one application of L_int; the bundles between vertices are the
intermediate labels, summed over.  The single-vertex tables of the paper give the
first-order factors that sit at each vertex.
"""
DY = 0.46
def f(x): return ("%.3f" % x).rstrip('0').rstrip('.') if abs(x) > 1e-9 else "0"
def ys(n): return [(n-1)/2*DY - k*DY for k in range(n)]

def bundle(x0, x1, styles, label=None, mark=0.5, lab_dy=0.0):
    """lines fanning out of the vertex at x0 and converging into the vertex at x1"""
    out = []
    n = len(styles)
    b = min(0.85, 0.35*(x1-x0))
    for st, y in zip(styles, ys(n)):
        if abs(y) < 1e-9:
            path = f"({f(x0)},0) -- ({f(x1)},0)"
        else:
            path = (f"({f(x0)},0) .. controls ({f(x0+0.45)},0) and ({f(x0+0.5)},{f(y)}) .. "
                    f"({f(x0+b)},{f(y)}) -- ({f(x1-b)},{f(y)}) .. controls "
                    f"({f(x1-0.5)},{f(y)}) and ({f(x1-0.45)},0) .. ({f(x1)},0)")
        out.append(f"  \\draw[{st},ld mark at={mark}] {path};")
    if label:
        out.append(f"  \\node[ld label,anchor=south,inner sep=1pt] at "
                   f"({f(0.5*(x0+x1))},{f(max(ys(n))+0.12+lab_dy)}) {{\\scriptsize {label}}};")
    return out

def openbundle(x_vertex, L, styles, label=None, side='in'):
    out = []
    n = len(styles)
    for st, y in zip(styles, ys(n)):
        if side == 'in':
            x0, x1 = x_vertex, x_vertex - L
            if abs(y) < 1e-9: path = f"({f(x0)},0) -- ({f(x1)},0)"
            else: path = (f"({f(x0)},0) .. controls ({f(x0-0.42)},0) and ({f(x0-0.5)},{f(y)}) .. "
                          f"({f(x0-0.9)},{f(y)}) -- ({f(x1)},{f(y)})")
            mark = 0.5 if abs(y) < 1e-9 else 0.68
        else:
            x0, x1 = x_vertex + L, x_vertex
            if abs(y) < 1e-9: path = f"({f(x0)},0) -- ({f(x1)},0)"
            else: path = (f"({f(x0)},{f(y)}) -- ({f(x1+0.9)},{f(y)}) .. controls "
                          f"({f(x1+0.5)},{f(y)}) and ({f(x1+0.42)},0) .. ({f(x1)},0)")
            mark = 0.5 if abs(y) < 1e-9 else 0.32
        out.append(f"  \\draw[{st},ld mark at={mark}] {path};")
    if label:
        x = x_vertex - L - 0.12 if side == 'in' else x_vertex + L + 0.12
        anc = "east" if side == 'in' else "west"
        out.append(f"  \\node[ld ext,anchor={anc}] at ({f(x)},0) {{{label}}};")
    return out

def chain(ext_in, segments, ox=0.0, oy=0.0, right=None, dx=2.3, Lext=1.7,
          seg_labels=None, vlabels=None, dots_after=None):
    """ext_in: styles of the observable lines; segments: list of style-lists for the
    bundles between consecutive vertices (len = #vertices - 1)."""
    nv = len(segments) + 1
    out = [f"\\begin{{scope}}[shift={{({f(ox)},{f(oy)})}}]"]
    out += openbundle(0.0, Lext, ext_in, None, 'in')
    for k, st in enumerate(segments):
        x0, x1 = k*dx, (k+1)*dx
        lab = seg_labels[k] if seg_labels else None
        if dots_after is not None and k == dots_after:
            out.append(f"  \\node[anchor=center] at ({f(0.5*(x0+x1))},0) {{$\\cdots$}};")
            continue
        out += bundle(x0, x1, st, lab)
    for k in range(nv):
        x = k*dx
        out.append(f"  \\node[ld vtx] at ({f(x)},0) {{}};")
        if vlabels and vlabels[k]:
            out.append(f"  \\node[ld order,fill=white,fill opacity=0.85,text opacity=1,"
                       f"rounded corners=1pt] at ({f(x)},0.15) {{${vlabels[k]}$}};")
    if right:
        out.append(f"  \\node[ld ext,anchor=west] at ({f((nv-1)*dx+0.22)},0) {{{right}}};")
    out.append("\\end{scope}")
    return "\n".join(out)

def txt(s, x, y, size="normalsize"):
    return f"\\node[anchor=center,font=\\{size}] at ({f(x)},{f(y)}) {{{s}}};"
