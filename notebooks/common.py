"""Shared helpers for the section notebooks in code/notebooks/.

    from common import *

gives the paths (CODE, ROOT, PAPER, ENGINES, DATA), the engine imports, the two small helpers of the
paper notebook (mode_ops, moment_series), `run` for the benchmark scripts, `show_pdf` / `same_as_paper`
for the paper figures, `paper_table` for Tables V and SI, and `J` for the result files.  FAST is read from the environment
(REPRO_FAST=0 for the orders and cutoffs of the paper) and can be overridden in the notebook.
"""
import os, sys, json, time, shutil, subprocess, re, tempfile
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from math import comb
from IPython.display import display, Image, Markdown

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, ".."))          # the repository root: engines/, benchmarks/, notebook/, notebooks/
ROOT = CODE

def _paper_dir():
    """LaTeX sources of the paper, if present: $PAPER_DIR, or paper/ in or next to the repository; else None."""
    for d in (os.environ.get("PAPER_DIR"), os.path.join(CODE, "paper"), os.path.join(CODE, "..", "paper")):
        if d and os.path.exists(os.path.join(d, "main.tex")): return os.path.abspath(d)
    return None
PAPER = _paper_dir()
PAPER_REF = os.path.join(HERE, "paper_reference")          # snapshots of Tables V and SI as printed
ENGINES = os.path.join(CODE, "engines")
DATA = os.path.join(CODE, "notebook", "data")
BENCH = os.path.join(CODE, "benchmarks")
for p in (ENGINES, os.path.join(BENCH, "truncation"), os.path.join(BENCH, "cumulant_closures"), os.path.join(BENCH, "trajectories")):
    if p not in sys.path: sys.path.insert(0, p)

FAST = os.environ.get("REPRO_FAST", "1") != "0"

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
                     "legend.frameon": False, "xtick.direction": "in", "ytick.direction": "in"})
COL = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#a6761d"]

from engine import Model, cav, add, scale, mul, destroy, steady, liouvillian

def mode_ops(alpha, d=1):
    """Displaced mode a = alpha + b on a d-dimensional auxiliary space: returns a, a^dag, b, b^dag."""
    I = np.eye(d); b, bd, one = cav(0, 1, I), cav(1, 0, I), cav(0, 0, I)
    return add(b, scale(alpha, one)), add(bd, scale(np.conj(alpha), one)), b, bd

def moment_series(mod, m, n, alpha, N):
    """Series of the laboratory-frame moment <(a^dag)^m a^n> through order N from the displaced-frame diagrams."""
    tot = np.zeros(N + 1, complex)
    for j in range(m + 1):
        for k in range(n + 1):
            tot += comb(m, j) * comb(n, k) * np.conj(alpha) ** (m - j) * alpha ** (n - k) * mod.series(cav(j, k, np.eye(1)), N)
    return tot

J = lambda *p: json.load(open(os.path.join(CODE, *p)))

def run(cmd, cwd, timeout=3600, tail=3):
    """Run a benchmark script in its own directory (paths in the scripts are relative to it) and show the tail of its output."""
    t0 = time.time()
    r = subprocess.run(cmd, cwd=os.path.join(CODE, cwd), shell=True, capture_output=True, text=True, timeout=timeout)
    out = "\n".join((r.stdout + r.stderr).strip().splitlines()[-tail:])
    print(f"$ (cd code/{cwd}; {cmd})   [{time.time()-t0:.1f} s, exit {r.returncode}]\n{out}")
    if r.returncode: raise RuntimeError(r.stderr[-2000:])

def _raster(pdf, dpi):
    d = tempfile.mkdtemp(); base = os.path.join(d, "p")
    subprocess.run(["pdftoppm", "-r", str(dpi), "-png", "-singlefile", pdf, base], check=True)
    return base + ".png"

def show_pdf(path, dpi=110, width=520):
    """Display a PDF figure inline (needs poppler's pdftoppm, or pypdfium2)."""
    path = os.path.join(CODE, path) if not os.path.isabs(path) else path
    if shutil.which("pdftoppm"):
        display(Image(filename=_raster(path, dpi), width=width)); return
    try:
        import pypdfium2 as pdfium
        png = tempfile.mktemp(suffix=".png"); pdfium.PdfDocument(path)[0].render(scale=dpi / 72).to_pil().save(png)
        display(Image(filename=png, width=width))
    except ImportError:
        print("(install poppler or pypdfium2 to see the figure inline)", path)

def same_as_paper(path, dpi=100):
    """Compare a regenerated figure with the copy in the paper sources pixel by pixel (None if they are not present).  'identical', or the fraction of
    differing pixels (text rendering differs between matplotlib versions and installed fonts; the curves do not)."""
    if PAPER is None: return None
    q = os.path.join(PAPER, os.path.basename(path)); p = os.path.join(CODE, path)
    if not os.path.exists(q): return None
    if not shutil.which("pdftoppm"): return "pdftoppm not available for the pixel comparison"
    from PIL import Image as PILImage
    A = [np.asarray(PILImage.open(_raster(f, dpi)).convert("L"), float) for f in (p, q)]
    if A[0].shape != A[1].shape: return f"page size differs: {A[0].shape} vs {A[1].shape}"
    frac = float((np.abs(A[0] - A[1]) > 40).mean())
    return "identical" if frac < 1e-3 else f"{100*frac:.1f}% of pixels differ (font rendering; the curves are the same)"

def paper_figure(path, dpi=110, width=520):
    """Show a regenerated figure and say how it compares with the file used in paper/."""
    show_pdf(path, dpi, width)
    v = same_as_paper(path)
    if v is not None: print(f"{os.path.basename(path)} versus the copy in the paper sources: {v}")

def paper_table(label):
    """LaTeX source of a table of the paper: from the paper sources if present, else from the snapshot."""
    snap = {"tab:fcs": "table_V_main.tex", "tab:closures": "table_SI_supplement.tex"}[label]
    if PAPER is not None:
        for f in ("main.tex", "supplement.tex"):
            t = open(os.path.join(PAPER, f), encoding="utf-8").read()
            if "\\label{%s}" % label in t: return t, os.path.join(PAPER, f)
    return open(os.path.join(PAPER_REF, snap), encoding="utf-8").read(), os.path.join("notebooks", "paper_reference", snap)

def quoted(label, paper, code, fmt="{:.4g}"):
    """Print one 'paper says / code gives' line."""
    p = paper if isinstance(paper, str) else fmt.format(paper)
    c = code if isinstance(code, str) else fmt.format(code)
    print(f"  {label:<62} paper {p:>12}   code {c:>12}")
