#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Draw the two documentation figures as SVG and PNG.

1. workflow: every analysis step, from the inputs to the key numbers.
2. how-it-works: how one run moves between Python, R and the output directory.

Every label is taken from docs/methods.md, src/irbdmirna/pipeline.py, src/irbdmirna/de.py and
the R scripts in src/irbdmirna/r. The two small plots in the workflow figure are schematic. They
are drawn from fixed synthetic shapes, and no data of any kind enters either figure.

    pip install cairosvg
    python tools/make_doc_figures.py

On macOS with cairo from Homebrew, prefix the command with
DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib.
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "docs" / "img"
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"
MONO = "Menlo, DejaVu Sans Mono, Consolas, monospace"
MATH = "STIXGeneral, STIX Two Text, Cambria Math, Times New Roman, serif"

INK, MUTED, LINE, PANEL, CARD = "#1B2230", "#5B6474", "#D6DBE4", "#F3F5F9", "#FFFFFF"
R_COL, PY_COL, FILE_COL, RED, OK = "#2E3A8C", "#B85A0D", "#6B7385", "#B4233A", "#0F7B6C"
# Group colours match src/irbdmirna/figures.py.
G_CTRL, G_PD, G_IRBD = "#0072B2", "#D55E00", "#009E73"
UP, DOWN, GREY = "#C0392B", "#0072B2", "#C8CDD6"

W, H = 2560, 1540
BODY, LH, ITEM_GAP = 15.5, 21.5, 7.0


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def t(x, y, s, size=BODY, *, weight=400, fill=INK, anchor="start", family=FONT, italic=False):
    style = ' font-style="italic"' if italic else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="{family}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{style}>{esc(s)}</text>'
    )


def box(x, y, w, h, *, fill=CARD, stroke=LINE, rx=10, dash=False, width=1.2):
    d = ' stroke-dasharray="6 4"' if dash else ""
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"{d}/>'
    )


def badge(x, y, label, colour=INK, r=14):
    return f'<circle cx="{x}" cy="{y}" r="{r}" fill="{colour}"/>' + t(
        x, y + 4.5, label, 12.5, weight=700, fill="#FFFFFF", anchor="middle"
    )


def pill(x, y, label, colour, *, anchor="end"):
    """Small filled tag. With anchor 'end', x is the right edge."""
    w = 12 + len(label) * 7.6
    x0 = x - w if anchor == "end" else x
    return box(x0, y - 14, w, 21, fill=colour, stroke=colour, rx=10.5) + t(
        x0 + w / 2, y + 1.5, label, 12, weight=700, fill="#FFFFFF", anchor="middle"
    )


def code(x, y, s, *, size=14):
    w = 16 + len(s) * size * 0.62
    return box(x, y - 17, w, 26, fill="#EEF1F6", stroke="#EEF1F6", rx=5) + t(
        x + 8, y + 1, s, size, family=MONO, fill=INK
    )


def arrow(x1, y1, x2, y2, colour="#6B7385", width=2.4):
    return (
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{colour}" '
        f'stroke-width="{width}" marker-end="url(#{"arrowR" if colour == R_COL else "arrow"})"/>'
    )


def wrap(text: str, width_px: float, size: float = BODY) -> list[str]:
    chars = max(10, int(width_px / (size * 0.52)))
    return textwrap.wrap(text, chars, break_on_hyphens=False)


def card(
    x,
    y,
    w,
    title,
    items,
    *,
    tag=None,
    extra=0.0,
    dash=False,
    border=None,
    title_fill=INK,
    step=None,
    step_colour=INK,
):
    """A white card with a title, wrapped body lines and an optional language tag.

    `items` is a list of sentences. Each sentence starts on a new line and wraps inside the card.
    Returns the SVG and the height, so that cards can be stacked.
    """
    body: list[str] = []
    yy = y + 68
    for item in items:
        for line in wrap(item, w - 40):
            body.append(t(x + 18, yy, line, BODY, fill=MUTED))
            yy += LH
        yy += ITEM_GAP
    end = yy - LH - ITEM_GAP + 10
    h = end - y + extra + 16
    svg = [box(x, y, w, h, dash=dash, stroke=border or LINE, width=1.8 if border else 1.2)]
    if step:
        svg.append(badge(x + 30, y + 32, step, step_colour, r=12))
        svg.append(t(x + 50, y + 38, title, 18.5, weight=700, fill=title_fill))
    else:
        svg.append(t(x + 18, y + 38, title, 18.5, weight=700, fill=title_fill))
    if tag:
        svg.append(pill(x + w - 14, y + 32, tag, R_COL if tag == "R" else PY_COL))
    svg.extend(body)
    return "".join(svg), h, end


def panel(x, y, w, h, number, title, subtitle):
    return (
        box(x, y, w, h, fill=PANEL, stroke="#E3E7EE", rx=18)
        + badge(x + 32, y + 44, number)
        + t(x + 56, y + 51, title, 21, weight=700)
        + t(x + 22, y + 84, subtitle, 16, fill=MUTED)
    )


def count_grid(x, y, cols=12, rows=4, cell=17):
    shades = ["#E4E8F2", "#C6CEE3", "#9AA7CB", "#6C7DB0", "#3F5194"]
    out = []
    for r in range(rows):
        for c in range(cols):
            k = int((math.sin(1.7 * c + 2.3 * r) + math.cos(0.9 * c * r + 1.1) + 2) / 4 * 4.99)
            out.append(
                f'<rect x="{x + c * (cell + 3)}" y="{y + r * (cell + 3)}" width="{cell}" '
                f'height="{cell}" rx="3" fill="{shades[k]}"/>'
            )
    return "".join(out)


def mini_volcano(x, y, w, h):
    """Schematic volcano: grey cloud, coloured hits beyond both thresholds."""
    out = [
        f'<line x1="{x}" y1="{y + h}" x2="{x + w}" y2="{y + h}" stroke="{LINE}" stroke-width="1.2"/>'
    ]
    cx = x + w / 2
    fc = w * 0.17
    yline = y + h * 0.58
    for xx in (cx - fc, cx + fc):
        out.append(
            f'<line x1="{xx:.1f}" y1="{y}" x2="{xx:.1f}" y2="{y + h}" stroke="#9AA0A6" '
            'stroke-width="1" stroke-dasharray="4 3"/>'
        )
    out.append(
        f'<line x1="{x}" y1="{yline:.1f}" x2="{x + w}" y2="{yline:.1f}" stroke="#9AA0A6" '
        'stroke-width="1" stroke-dasharray="4 3"/>'
    )
    for i in range(70):
        u = math.sin(i * 12.9898) * 43758.5453
        u -= math.floor(u)
        v = math.sin(i * 78.233) * 12345.678
        v -= math.floor(v)
        dx = (u - 0.5) * 2
        lift = abs(dx) ** 1.6 * (0.55 + 0.6 * v)
        px = cx + dx * w * 0.47
        py = y + h - lift * h * 0.95
        hit = abs(dx) * w * 0.47 > fc and py < yline
        colour = (UP if dx > 0 else DOWN) if hit else GREY
        r = 3.6 if hit else 2.8
        out.append(
            f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{r}" fill="{colour}" stroke="#FFFFFF" stroke-width="0.6"/>'
        )
    out.append(t(x + w, y + h + 18, "log2 fold change", 12, fill=MUTED, anchor="end"))
    out.append(t(x, y + h + 18, "schematic", 12, fill=MUTED, italic=True))
    return "".join(out)


def chromosome(x, y, w):
    """Schematic 14q32 domain: DLK1 and DIO3 flank the cluster on the forward strand."""
    out = [
        f'<rect x="{x}" y="{y}" width="{w}" height="12" rx="6" fill="#E4E8F2" stroke="{LINE}"/>',
        f'<rect x="{x + 6}" y="{y}" width="26" height="12" rx="3" fill="{FILE_COL}"/>',
        f'<rect x="{x + w - 32}" y="{y}" width="26" height="12" rx="3" fill="{FILE_COL}"/>',
        t(x + 19, y + 34, "DLK1", 12, fill=MUTED, anchor="middle"),
        t(x + w - 19, y + 34, "DIO3", 12, fill=MUTED, anchor="middle"),
    ]
    for k in range(26):
        xx = x + 50 + k * (w - 100) / 25
        hh = 10 + (k % 3) * 4
        out.append(
            f'<line x1="{xx:.1f}" y1="{y - hh}" x2="{xx:.1f}" y2="{y}" stroke="{R_COL}" stroke-width="2"/>'
        )
    out.append(arrow(x + w * 0.35, y + 30, x + w * 0.65, y + 30, colour=R_COL, width=1.6))
    out.append(t(x + w / 2, y + 48, "forward strand", 12, fill=MUTED, anchor="middle"))
    return "".join(out)


def figure() -> str:
    s: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
        'role="img" aria-label="Workflow of the irbd-plasma-mirna pipeline">',
        "<title>Workflow of the irbd-plasma-mirna pipeline</title>",
        f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>',
        "<defs>"
        + "".join(
            f'<marker id="{mid}" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" '
            f'markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="{col}"/>'
            "</marker>"
            for mid, col in (("arrow", "#6B7385"), ("arrowR", R_COL))
        )
        + "</defs>",
    ]
    # Header.
    s.append(t(40, 62, "irbd-plasma-mirna workflow", 34, weight=700))
    s.append(
        t(
            40,
            96,
            "Plasma microRNA counts of iRBD, de novo PD and healthy controls, "
            "from the count matrix to the key numbers of the manuscript",
            18,
            fill=MUTED,
        )
    )
    s.append(code(40, 132, "irbd-mirna run --config config/config.example.yaml"))
    lx = 1950
    for label, colour, text in (
        ("R", R_COL, "runs in R"),
        ("Py", PY_COL, "runs in Python"),
        ("", FILE_COL, "files"),
    ):
        s.append(f'<circle cx="{lx}" cy="64" r="14" fill="{colour}"/>')
        if label:
            s.append(t(lx, 68.5, label, 12, weight=700, fill="#FFFFFF", anchor="middle"))
        s.append(t(lx + 24, 70, text, 18))
        lx += 52 + len(text) * 9.5

    top, ph = 170, 1040
    cols = [(40, 330), (400, 330), (760, 540), (1330, 820), (2180, 340)]
    titles = [
        ("1", "INPUTS", "plain tab-separated files"),
        ("2", "CHECKS", "io.py, before any model"),
        ("3", "DIFFERENTIAL EXPRESSION", "DESeq2 in R, one fit per comparison"),
        ("4", "DOWNSTREAM ANALYSES", "SciPy and statsmodels in Python"),
        ("5", "OUTPUTS", "tables, figures and key numbers"),
    ]
    for (x, w), (n, title, sub) in zip(cols, titles, strict=True):
        s.append(panel(x, top, w, ph, n, title, sub))
    for (x, w), (nx, _) in zip(cols[:-1], cols[1:], strict=True):
        s.append(arrow(x + w + 3, top + ph / 2, nx - 3, top + ph / 2, width=3))

    gap, y0 = 18, top + 110

    # 1. Inputs.
    x, w = cols[0][0] + 18, cols[0][1] - 36
    y = y0
    c, h, by = card(x, y, w, "Count matrix", ["miRNA by library, miRMaster 2.0 counts"], extra=84)
    s.append(c + count_grid(x + 18, by + 4))
    y += h + gap
    c, h, by = card(
        x,
        y,
        w,
        "Sample sheet",
        ["library_id, group, age, sex, batch, qc_pass; converter for the iRBD comparison"],
        extra=112,
    )
    s.append(c)
    groups = [
        (G_IRBD, True, "iRBD"),
        (G_PD, True, "de novo PD"),
        (G_CTRL, True, "high-PRS controls"),
        (G_CTRL, False, "low-PRS controls"),
    ]
    for k, (colour, filled, label) in enumerate(groups):
        gy = by + 18 + k * 27
        fill = colour if filled else "#FFFFFF"
        s.append(
            f'<circle cx="{x + 28}" cy="{gy - 5}" r="7" fill="{fill}" stroke="{colour}" stroke-width="2"/>'
        )
        s.append(t(x + 44, gy, label, 15.5))
    y += h + gap
    c, h, _ = card(
        x, y, w, "Clinical table", ["optional; continuous and categorical variables"], dash=True
    )
    s.append(c)
    y += h + gap
    c, h, _ = card(
        x,
        y,
        w,
        "Resources",
        [
            "30 predefined brain-enriched miRNAs",
            "miRBase v22 table of the DLK1-DIO3 domain, GRCh38",
        ],
    )
    s.append(c)

    # 2. Checks.
    x, w = cols[1][0] + 18, cols[1][1] - 36
    y = y0
    for title, items in (
        (
            "Input checks",
            [
                "required columns and known group labels",
                "numeric, unique count columns; every library in the count matrix",
            ],
        ),
        ("Quality control", ["only libraries with qc_pass enter the analysis set"]),
        ("Rounding", ["non-integer counts rounded, as DESeq2 models integer counts"]),
        (
            "Summaries",
            [
                "miRNAs detected in at least one and in at least half of the libraries",
                "median and range of miRNA reads per library",
                "reads per million mapped miRNA reads",
            ],
        ),
    ):
        c, h, _ = card(x, y, w, title, items, tag="Py")
        s.append(c)
        if title != "Summaries":
            s.append(arrow(x + w / 2, y + h + 1, x + w / 2, y + h + gap - 1, width=2))
        y += h + gap

    # 3. Differential expression.
    x, w = cols[2][0] + 18, cols[2][1] - 36
    y = y0
    rows = [
        ("High-PRS vs low-PRS controls", "edgeR filterByExpr"),
        ("iRBD vs controls", "edgeR filterByExpr"),
        ("iRBD vs de novo PD", "edgeR filterByExpr"),
        ("de novo PD vs controls", "≥ 10 counts in ≥ 3 libraries"),
    ]
    c, h, by = card(
        x,
        y,
        w,
        "Four comparisons and their pre-filter",
        ["controls pooled except in the first comparison"],
        tag="R",
        extra=len(rows) * 32 + 4,
    )
    s.append(c)
    for k, (cmp_, rule) in enumerate(rows):
        ry = by + 22 + k * 32
        s.append(
            f'<line x1="{x + 18}" y1="{ry - 22}" x2="{x + w - 18}" y2="{ry - 22}" stroke="#EEF1F6" stroke-width="1.2"/>'
        )
        s.append(t(x + 18, ry, cmp_, 15.5, weight=700))
        s.append(t(x + w - 18, ry, rule, 15, fill=MUTED, anchor="end"))
    y += h + gap
    s.append(arrow(x + w / 2, y - gap + 1, x + w / 2, y - 1, width=2))
    c, h, by = card(
        x,
        y,
        w,
        "Model",
        [
            "age centred and scaled; sequencing batch not in the model",
            "median-of-ratios size factors; dispersions shrunk towards a fitted trend",
        ],
        tag="R",
        extra=40,
    )
    s.append(c)
    s.append(t(x + 18, by + 28, "design  ~ age + sex + group", 21, family=MATH, italic=True))
    y += h + gap
    s.append(arrow(x + w / 2, y - gap + 1, x + w / 2, y - 1, width=2))
    c, h, _ = card(
        x,
        y,
        w,
        "Test and shrinkage",
        [
            "Wald test of the group coefficient",
            "independent filtering at alpha 0.1, then Benjamini-Hochberg adjustment",
            "apeglm shrinkage of the log2 fold change",
        ],
        tag="R",
    )
    s.append(c)
    y += h + gap
    s.append(arrow(x + w / 2, y - gap + 1, x + w / 2, y - 1, width=2))
    c, h, by = card(
        x,
        y,
        w,
        "Dysregulation rule",
        ["FDR below 0.05 and absolute shrunken log2 fold change of at least 0.5"],
        tag="Py",
        border=RED,
        title_fill=RED,
        extra=150,
    )
    s.append(c + mini_volcano(x + 40, by + 12, w - 80, 108))

    # 4. Downstream analyses, two sub-columns.
    x0, wc = cols[3][0] + 18, (cols[3][1] - 54) / 2
    s.append(t(x0 + 4, y0 + 4, "FROM THE DESEQ2 TABLES", 12.5, weight=700, fill=MUTED))
    s.append(t(x0 + wc + 22, y0 + 4, "FROM THE CHECKED INPUTS", 12.5, weight=700, fill=MUTED))
    y = y0 + 18
    c, h, by = card(
        x0,
        y,
        wc,
        "14q32 cluster",
        [
            "members by genomic position, not by name",
            "tested and dysregulated members per comparison",
            "hsa-mir-1247, on the reverse strand, reported apart",
        ],
        extra=70,
    )
    s.append(c + chromosome(x0 + 26, by + 26, wc - 52))
    y += h + gap
    c, h, _ = card(
        x0,
        y,
        wc,
        "Brain-enriched panel",
        [
            "30 miRNAs followed in the three main comparisons",
            "shared candidates: dysregulated in both iRBD comparisons, same direction",
            "one-sided Fisher's exact test, with and without the 14q32 cluster",
        ],
    )
    s.append(c)
    y += h + gap
    c, h, _ = card(
        x0,
        y,
        wc,
        "Haemolysis index and cluster score",
        [
            "index: log2 of miR-451a over miR-23a-3p, reads per million plus one",
            "score: mean standardised expression of the cluster miRNAs tested in iRBD vs controls",
            "two-sided Mann-Whitney U; Spearman correlation within each group",
            "linear models of the score: group, plus index, plus age and sex",
        ],
    )
    s.append(c)

    x1 = x0 + wc + 18
    y = y0 + 18
    c, h, _ = card(
        x1,
        y,
        wc,
        "Converters within iRBD",
        [
            "DESeq2 size factors of the iRBD libraries",
            "miRNAs detected in at least 20% of iRBD libraries, and in at least 3",
            "two-sided Mann-Whitney U, Benjamini-Hochberg, the same rule",
            "repeated on 200 datasets with permuted labels",
        ],
    )
    s.append(c)
    y += h + gap
    c, h, _ = card(
        x1,
        y,
        wc,
        "Batch analysis",
        [
            "chi-square test of batch by group, Cramer's V",
            "principal components of log2 expression, miRNAs detected in at least half the libraries",
            "Kruskal-Wallis test of the three leading components",
            "per-miRNA F test of group plus batch against group, Benjamini-Hochberg",
        ],
    )
    s.append(c)
    y += h + gap
    c, h, _ = card(
        x1,
        y,
        wc,
        "Clinical comparisons",
        [
            "Welch's t test for continuous and Fisher's exact test for categorical variables",
            "Benjamini-Hochberg within each pairwise comparison",
            "pairwise deletion; missing values listed per group",
        ],
        dash=True,
    )
    s.append(c)
    s.append(
        t(
            x0 + 4,
            top + ph - 22,
            "All six analyses run in Python. The clinical comparisons run when a clinical table is given.",
            14,
            fill=MUTED,
        )
    )

    # 5. Outputs.
    x, w = cols[4][0] + 18, cols[4][1] - 36
    y = y0
    figs = [
        ("1", "volcano plots, main comparisons"),
        ("2", "brain-enriched panel"),
        ("S1", "high-PRS vs low-PRS controls"),
        ("S2", "converters vs non-converters"),
        ("S3", "batch principal components"),
    ]
    c, h, by = card(
        x, y, w, "Figures", ["PNG and PDF, in figures/"], tag="Py", extra=len(figs) * 30
    )
    s.append(c)
    for k, (num, label) in enumerate(figs):
        fy = by + 22 + k * 30
        s.append(pill(x + 18, fy, num, FILE_COL, anchor="start"))
        s.append(t(x + 64, fy + 1, label, 14.5))
    y += h + gap
    c, h, _ = card(
        x,
        y,
        w,
        "Tables",
        [
            "one folder per analysis: de, brain_enriched, cluster_14q32, haemolysis, converters, batch, clinical"
        ],
    )
    s.append(c)
    y += h + gap
    c, h, by = card(
        x,
        y,
        w,
        "Key numbers",
        [
            "summary/key_numbers.tsv: every computed number beside the value in the manuscript",
            "with --strict, a difference fails the run",
        ],
        tag="Py",
        border=RED,
        title_fill=RED,
        extra=36,
    )
    s.append(c)
    sx = x + 18
    for label, colour in (("match", OK), ("differs", RED), ("not reported", FILE_COL)):
        s.append(pill(sx, by + 22, label, colour, anchor="start"))
        sx += 22 + len(label) * 7.6 + 8
    y += h + gap
    c, h, _ = card(
        x,
        y,
        w,
        "Run record",
        ["summary/run_info.yaml: software versions, input checksums and the configuration"],
    )
    s.append(c)

    # Reference environment band.
    by0, bh = top + ph + 30, 232
    s.append(box(40, by0, W - 80, bh, fill=PANEL, stroke="#E3E7EE", rx=18))
    s.append(t(76, by0 + 46, "REFERENCE ENVIRONMENT AND CHECKS", 21, weight=700))
    s.append(
        t(
            76,
            by0 + 78,
            "the versions of the manuscript, checked at build time and on every push",
            16,
            fill=MUTED,
        )
    )
    steps = [
        ("Docker image", "R 4.6.1, Bioconductor 3.23"),
        ("Version gate in R", "DESeq2 1.52, apeglm 1.34, edgeR"),
        ("Version gate in Python", "Python 3.10.12, SciPy 1.15.3, statsmodels 0.14.6"),
        ("Tests", "ruff, pytest and testthat"),
        ("Example run", "synthetic data in examples/data"),
        ("Continuous integration", "pip, conda and Docker on every push"),
    ]
    bw, bg = (W - 152 - 5 * 44) / 6, 44
    bx = 76
    for k, (title, sub) in enumerate(steps):
        xx = bx + k * (bw + bg)
        s.append(box(xx, by0 + 104, bw, 94))
        s.append(t(xx + 18, by0 + 138, title, 17, weight=700))
        for j, line in enumerate(wrap(sub, bw - 36, 15)):
            s.append(t(xx + 18, by0 + 164 + j * 20, line, 15, fill=MUTED))
        if k < len(steps) - 1:
            s.append(arrow(xx + bw + 6, by0 + 146, xx + bw + bg - 6, by0 + 146, width=2))
    s.append("</svg>")
    return "".join(s)


def elbow(points, colour="#6B7385", width=2.2, dash=False):
    """Arrow along a polyline, for connectors that change direction."""
    pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    d = ' stroke-dasharray="5 4"' if dash else ""
    return (
        f'<polyline points="{pts}" fill="none" stroke="{colour}" stroke-width="{width}"{d} '
        f'stroke-linejoin="round" marker-end="url(#{"arrowR" if colour == R_COL else "arrow"})"/>'
    )


def files_card(x, y, w, paths):
    h = 20 + len(paths) * 24 + 8
    out = [box(x, y, w, h, fill="#FFFFFF", stroke=LINE, dash=True)]
    for k, path in enumerate(paths):
        out.append(t(x + 16, y + 34 + k * 24, path, 13.5, family=MONO, fill=INK))
    return "".join(out), h


def how_it_works() -> tuple[str, int, int]:
    width, height = 2560, 1040
    s: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-label="How one run of irbd-mirna works">',
        "<title>How one run of irbd-mirna works</title>",
        f'<rect width="{width}" height="{height}" fill="#FFFFFF"/>',
        "<defs>"
        + "".join(
            f'<marker id="{mid}" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" '
            f'markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="{col}"/>'
            "</marker>"
            for mid, col in (("arrow", "#6B7385"), ("arrowR", R_COL))
        )
        + "</defs>",
        t(40, 62, "How one run works", 34, weight=700),
        t(
            40,
            96,
            "Python drives the run and calls R once for DESeq2. The two parts exchange plain "
            "tab-separated files and one YAML job file.",
            18,
            fill=MUTED,
        ),
    ]
    lx = 1950
    for label, colour, text in (
        ("R", R_COL, "runs in R"),
        ("Py", PY_COL, "runs in Python"),
        ("", FILE_COL, "files"),
    ):
        s.append(f'<circle cx="{lx}" cy="64" r="14" fill="{colour}"/>')
        if label:
            s.append(t(lx, 68.5, label, 12, weight=700, fill="#FFFFFF", anchor="middle"))
        s.append(t(lx + 24, 70, text, 18))
        lx += 52 + len(text) * 9.5

    col0, cw, cg = 300, 340, 22
    cx = [col0 + k * (cw + cg) for k in range(6)]
    lanes = [
        ("YOU", "the command line", FILE_COL, 130, 140),
        ("PYTHON", "irbd-mirna", PY_COL, 290, 210),
        ("R", "Rscript run_deseq2.R", R_COL, 520, 250),
        ("OUTPUT DIRECTORY", "output_dir of the configuration", FILE_COL, 790, 210),
    ]
    for title, sub, colour, y, h in lanes:
        s.append(box(40, y, width - 80, h, fill=PANEL, stroke="#E3E7EE", rx=18))
        s.append(f'<rect x="40" y="{y + 18}" width="6" height="{h - 36}" rx="3" fill="{colour}"/>')
        s.append(t(66, y + 46, title, 19, weight=700))
        for j, line in enumerate(wrap(sub, 200, 15)):
            s.append(t(66, y + 72 + j * 20, line, 15, fill=MUTED))

    # You.
    s.append(code(cx[0], 196, "irbd-mirna run --config config.yaml"))
    s.append(
        t(
            cx[0] + 360,
            201,
            "add --strict to stop on any difference from the manuscript",
            14.5,
            fill=MUTED,
        )
    )
    c, h_exit, _ = card(
        cx[5], 148, cw, "Exit", ["code 0 and the line Finished, or code 1 with the reason"]
    )
    s.append(c)

    # Python steps.
    py_y, r_y = 310, 540
    steps = {
        0: (
            "Read and check the inputs",
            [
                "count matrix, sample sheet and clinical table",
                "quality-control filter, rounding, detection summary",
            ],
        ),
        1: (
            "Hand over to R",
            [
                "write the analysis set and a job file with the four comparisons, covariates and "
                "pre-filters",
            ],
        ),
        3: (
            "Downstream analyses",
            [
                "dysregulation rule, 14q32 cluster, brain-enriched panel, haemolysis, converters, "
                "batch and clinical comparisons",
            ],
        ),
        4: ("Figures", ["Figures 1 and 2 and Supplementary Figures S1 to S3, as PNG and PDF"]),
        5: (
            "Key numbers",
            [
                "every number beside the value in the manuscript",
            ],
        ),
    }
    heights = {}
    for k, (title, items) in steps.items():
        c, h, _ = card(cx[k], py_y, cw, title, items, tag="Py", step=str(k + 1))
        s.append(c)
        heights[k] = h

    # R step.
    c, h_r, _ = card(
        cx[2],
        r_y,
        cw,
        "DESeq2 per comparison",
        [
            "pre-filter, fit of ~ age + sex + group, Wald test",
            "independent filtering, Benjamini-Hochberg, apeglm shrinkage",
            "then size factors over all libraries and within iRBD",
        ],
        tag="R",
        step="3",
        step_colour=R_COL,
    )
    s.append(c)

    # Output files.
    out_y = 812
    outputs = {
        0: ["qc/library_depth.tsv"],
        1: [
            "work/counts_analysis_set.tsv",
            "work/samples_analysis_set.tsv",
            "work/deseq2_job.yaml",
        ],
        2: [
            "de/<comparison>_all.tsv",
            "de/comparison_summary.tsv",
            "de/size_factors_*.tsv",
            "de/normalised_counts_all_samples.tsv",
            "de/R_session_info.txt",
        ],
        3: [
            "de/<comparison>_annotated.tsv",
            "de/<comparison>_significant.tsv",
            "brain_enriched/, cluster_14q32/",
            "haemolysis/, converters/",
            "batch/, clinical/",
        ],
        4: ["figures/*.png", "figures/*.pdf"],
        5: ["summary/key_numbers.tsv", "summary/run_info.yaml"],
    }
    for k, paths in outputs.items():
        c, _ = files_card(cx[k], out_y, cw, paths)
        s.append(c)
        top = (r_y + h_r) if k == 2 else (py_y + heights[k])
        s.append(
            elbow([(cx[k] + cw / 2, top + 2), (cx[k] + cw / 2, out_y - 4)], dash=True, width=1.6)
        )

    # Control flow.
    mid = cx[0] + 120
    s.append(elbow([(mid, 211), (mid, py_y - 4)]))
    s.append(elbow([(cx[0] + cw + 2, py_y + 60), (cx[1] - 4, py_y + 60)]))
    s.append(
        elbow(
            [
                (cx[1] + cw + 2, py_y + 60),
                (cx[2] + cw / 2, py_y + 60),
                (cx[2] + cw / 2, r_y - 4),
            ],
            colour=R_COL,
        )
    )
    s.append(t(cx[2] + cw / 2 + 12, py_y + 110, "job file", 14, fill=R_COL, weight=700))
    s.append(
        elbow(
            [
                (cx[2] + cw + 2, r_y + h_r / 2),
                (cx[3] + cw / 2, r_y + h_r / 2),
                (cx[3] + cw / 2, py_y + heights[3] + 4),
            ],
            colour=R_COL,
        )
    )
    s.append(t(cx[3] + cw / 2 + 12, r_y + 30, "DESeq2 tables", 14, fill=R_COL, weight=700))
    s.append(elbow([(cx[3] + cw + 2, py_y + 60), (cx[4] - 4, py_y + 60)]))
    s.append(elbow([(cx[4] + cw + 2, py_y + 60), (cx[5] - 4, py_y + 60)]))
    s.append(elbow([(cx[5] + cw / 2, py_y - 2), (cx[5] + cw / 2, 160 + h_exit + 4)]))
    s.append("</svg>")
    return "".join(s), width, height


def main() -> None:
    import cairosvg

    OUT.mkdir(parents=True, exist_ok=True)
    for name, (text, width) in {
        "workflow": (figure(), W),
        "how-it-works": how_it_works()[:2],
    }.items():
        (OUT / f"{name}.svg").write_text(text, encoding="utf-8")
        cairosvg.svg2png(
            bytestring=text.encode(), write_to=str(OUT / f"{name}.png"), output_width=width * 1.5
        )
        print(f"Wrote {OUT / name}.svg and .png")


if __name__ == "__main__":
    main()
