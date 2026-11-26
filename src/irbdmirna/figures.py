# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Figures: volcano plots, brain-enriched panel, converter comparison and batch PCA.

The figures draw the tables of record only. They run no test of their own: thresholds, adjusted
p values and the dysregulation calls come from the DESeq2 and downstream tables.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# Red is up and blue is down, as in the figure captions of the manuscript. Context stays grey.
UP, DOWN, NEUTRAL = "#C0392B", "#0072B2", "#BDBDBD"
INK, MUTED, RULE = "#1B2230", "#5B6474", "#9AA0A6"
GROUP_STYLE = {
    "controls": ("#0072B2", "o"),
    "de novo PD": ("#D55E00", "s"),
    "iRBD": ("#009E73", "^"),
}
CYCLE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#000000", "#F0E442"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]
# Batches get colours that carry no group meaning.
BATCH_CYCLE = ["#CC79A7", "#56B4E9", "#E69F00", "#000000", "#999999", "#F0E442"]


def _group_style(label: str, k: int) -> tuple[str, str, bool]:
    """Colour, marker and fill of a study group, consistent with the other figures."""
    low = label.lower()
    if "irbd" in low:
        return GROUP_STYLE["iRBD"][0], GROUP_STYLE["iRBD"][1], True
    if "pd" in low:
        return GROUP_STYLE["de novo PD"][0], GROUP_STYLE["de novo PD"][1], True
    if "low" in low:
        return GROUP_STYLE["controls"][0], "o", False
    if "high" in low or "control" in low or "hc" in low:
        return GROUP_STYLE["controls"][0], "o", True
    return CYCLE[k % len(CYCLE)], MARKERS[k % len(MARKERS)], True


STYLE = {
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "bold",
    "axes.labelsize": 8,
    "axes.labelcolor": INK,
    "axes.edgecolor": INK,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "irbd-plasma-mirna",
}


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", metadata={"CreationDate": None})
    plt.close(fig)


def _fdr_line(table: pd.DataFrame, fdr: float) -> float | None:
    """Largest unadjusted p value whose adjusted p value is below the FDR threshold."""
    hits = table.loc[table["padj"] < fdr, "pvalue"]
    return float(hits.max()) if len(hits) else None


def _spread(wanted: list[float], gap: float, low: float, high: float) -> list[float]:
    """Move sorted positions apart so neighbours are at least `gap` apart, inside [low, high]."""
    pos = list(wanted)
    for _ in range(50):
        moved = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] < gap:
                shift = (gap - (pos[i] - pos[i - 1])) / 2
                pos[i - 1] -= shift
                pos[i] += shift
                moved = True
        if pos:
            if pos[0] < low:
                pos = [p + (low - pos[0]) for p in pos]
            if pos[-1] > high:
                pos = [p - (pos[-1] - high) for p in pos]
        if not moved:
            break
    return pos


def _label_points(ax: plt.Axes, points: pd.DataFrame, names: list[str]) -> float:
    """Label points without overlap.

    Labels of each side stack in one column just outside the outermost labelled point, ordered
    by height, so that leader lines do not cross. Returns how far, in pixels, the widest label
    runs past the axes, so that the caller can widen the x range and draw again.
    """
    if not len(points):
        return 0.0
    fig = ax.figure
    fig.canvas.draw()
    to_px = ax.transData
    to_data = ax.transData.inverted()
    y0, y1 = ax.get_ylim()
    px_low = to_px.transform((0, y0))[1]
    px_high = to_px.transform((0, y1))[1]
    gap = 9.0 * fig.dpi / 72
    offset = 14 * fig.dpi / 72
    texts = []
    for side in (-1, 1):
        sel = [
            (name, x, y)
            for name, x, y in zip(names, points["x"], points["y"], strict=True)
            if (x >= 0) == (side > 0)
        ]
        if not sel:
            continue
        sel.sort(key=lambda r: r[2])
        px = [to_px.transform((x, y)) for _, x, y in sel]
        column = (max if side > 0 else min)(q[0] for q in px) + side * offset
        target = _spread([q[1] for q in px], gap, px_low + gap / 2, px_high - gap * 2.2)
        for (name, x, y), ty in zip(sel, target, strict=True):
            tx, tyd = to_data.transform((column, ty))
            texts.append(
                ax.annotate(
                    name,
                    xy=(x, y),
                    xytext=(tx, tyd),
                    fontsize=6.5,
                    color=INK,
                    ha="left" if side > 0 else "right",
                    va="center",
                    arrowprops={
                        "arrowstyle": "-",
                        "lw": 0.4,
                        "color": MUTED,
                        "shrinkA": 1,
                        "shrinkB": 2.5,
                    },
                )
            )
    fig.canvas.draw()
    frame = ax.get_window_extent()
    renderer = fig.canvas.get_renderer()
    over = 0.0
    for text in texts:
        ext = text.get_window_extent(renderer)
        over = max(over, frame.x0 - ext.x0, ext.x1 - frame.x1)
    return over


def volcano(
    table: pd.DataFrame,
    panel: list[str],
    title: str,
    path: Path,
    fdr: float = 0.05,
    abs_log2fc: float = 0.5,
    label_all_panel: bool = False,
) -> None:
    """Volcano plot. The y axis shows unadjusted p values; the dashed line marks FDR 0.05."""
    t = table.dropna(subset=["pvalue"]).copy()
    floor = t.loc[t["pvalue"] > 0, "pvalue"].min() if (t["pvalue"] > 0).any() else 1e-300
    t["y"] = -np.log10(t["pvalue"].clip(lower=floor))
    in_panel = t.index.isin(panel)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(4.6, 4.0))
        rest = t[t["direction"] == "none"]
        ax.scatter(
            rest["log2FoldChange"],
            rest["y"],
            s=6,
            c=NEUTRAL,
            alpha=0.55,
            linewidths=0,
            rasterized=True,
        )
        for direction, colour in (("up", UP), ("down", DOWN)):
            sub = t[(t["direction"] == direction) & ~in_panel]
            ax.scatter(sub["log2FoldChange"], sub["y"], s=10, c=colour, alpha=0.45, linewidths=0)
        if label_all_panel:
            quiet = t[in_panel & (t["direction"] == "none")]
            ax.scatter(
                quiet["log2FoldChange"],
                quiet["y"],
                s=22,
                facecolors="white",
                edgecolors=INK,
                linewidths=0.7,
            )
        for direction, colour in (("up", UP), ("down", DOWN)):
            sub = t[(t["direction"] == direction) & in_panel]
            ax.scatter(
                sub["log2FoldChange"], sub["y"], s=28, c=colour, edgecolors=INK, linewidths=0.7
            )
        for x in (-abs_log2fc, abs_log2fc):
            ax.axvline(x, ls=(0, (3, 2)), lw=0.6, color=RULE, zorder=0)
        line = _fdr_line(t, fdr)
        if line is not None:
            ax.axhline(-np.log10(line), ls=(0, (3, 2)), lw=0.6, color=RULE, zorder=0)
        lim = max(float(t["log2FoldChange"].abs().max()) * 1.25, abs_log2fc * 3)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(0, float(t["y"].max()) * 1.08 + 0.2)
        if line is not None:
            ax.text(
                1.01,
                -np.log10(line),
                f"FDR {fdr:g}",
                transform=ax.get_yaxis_transform(),
                fontsize=6.5,
                color=MUTED,
                ha="left",
                va="center",
                clip_on=False,
            )
        ax.set_xlabel("log$_2$ fold change")
        ax.set_ylabel("$-$log$_{10}$ p value (unadjusted)")
        ax.set_title(title, loc="left", color=INK)
        n_up, n_down = int((t["direction"] == "up").sum()), int((t["direction"] == "down").sum())
        ax.text(
            0.99,
            0.99,
            f"{n_up} up",
            transform=ax.transAxes,
            ha="right",
            va="top",
            color=UP,
            fontsize=8,
        )
        ax.text(
            0.01,
            0.99,
            f"{n_down} down",
            transform=ax.transAxes,
            ha="left",
            va="top",
            color=DOWN,
            fontsize=8,
        )
        labelled = t[in_panel] if label_all_panel else t[in_panel & (t["direction"] != "none")]
        points = pd.DataFrame({"x": labelled["log2FoldChange"], "y": labelled["y"]})
        names = [n.replace("hsa-", "") for n in labelled.index]
        for attempt in range(6):
            over = _label_points(ax, points, names)
            if over <= 0 or attempt == 5:
                break
            for text in list(ax.texts):
                if text.get_text() in names:
                    text.remove()
            lim *= 1 + 2.2 * (over + 4) / ax.get_window_extent().width
            ax.set_xlim(-lim, lim)
        handles = [
            Line2D(
                [],
                [],
                ls="",
                marker="o",
                ms=4.5,
                mfc=UP,
                mec=INK,
                mew=0.7,
                label="brain-enriched, up",
            ),
            Line2D(
                [],
                [],
                ls="",
                marker="o",
                ms=4.5,
                mfc=DOWN,
                mec=INK,
                mew=0.7,
                label="brain-enriched, down",
            ),
        ]
        if label_all_panel:
            handles.append(
                Line2D(
                    [],
                    [],
                    ls="",
                    marker="o",
                    ms=4.5,
                    mfc="white",
                    mec=INK,
                    mew=0.7,
                    label="brain-enriched, not called",
                )
            )
        ax.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.14),
            ncol=len(handles),
            fontsize=6.5,
        )
        _save(fig, path)


def _stars(table: pd.DataFrame, name: str) -> str | None:
    if name not in table.index or not bool(table.loc[name, "dysregulated"]):
        return None
    p = table.loc[name, "padj"]
    return "***" if p < 0.001 else "**" if p < 0.01 else "*"


def panel_figure(
    normalised: pd.DataFrame,
    group3: pd.Series,
    tables: dict[tuple[str, str], pd.DataFrame],
    panel: list[str],
    path: Path,
) -> None:
    """One panel per brain-enriched miRNA: log2(DESeq2-normalised count + 1) by group.

    Black lines show the median and the interquartile range. Asterisks mark comparisons
    that meet the dysregulation rule: * FDR below 0.05, ** below 0.01, *** below 0.001.
    """
    order = ["controls", "de novo PD", "iRBD"]
    ticks = ["controls", "de novo\nPD", "iRBD"]
    n_col = 6
    n_row = int(np.ceil(len(panel) / n_col))
    rng = np.random.RandomState(2)
    groups = group3.loc[normalised.columns].to_numpy()
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(n_row, n_col, figsize=(n_col * 1.9, n_row * 1.95 + 0.5))
        axes = np.atleast_1d(axes).ravel()
        for ax, name in zip(axes, panel, strict=False):
            ax.set_title(name.replace("hsa-", ""), fontsize=8, loc="left", color=INK)
            if name not in normalised.index:
                ax.text(
                    0.5,
                    0.5,
                    "not in count table",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=MUTED,
                )
                ax.set_axis_off()
                continue
            y = np.log2(normalised.loc[name].to_numpy(dtype=float) + 1.0)
            for j, g in enumerate(order):
                yy = y[groups == g]
                colour, marker = GROUP_STYLE[g]
                ax.scatter(
                    rng.normal(j, 0.07, len(yy)),
                    yy,
                    s=7,
                    color=colour,
                    marker=marker,
                    alpha=0.75,
                    linewidths=0,
                )
                if len(yy):
                    q1, q2, q3 = np.percentile(yy, [25, 50, 75])
                    ax.plot(
                        [j - 0.26, j + 0.26], [q2, q2], color=INK, lw=1.3, solid_capstyle="butt"
                    )
                    ax.plot([j, j], [q1, q3], color=INK, lw=0.9)
            top = float(y.max())
            step = max(top - float(y.min()), 1.0) * 0.09
            marks = [
                ((order.index(a), order.index(b)), s)
                for (a, b), t in tables.items()
                if (s := _stars(t, name)) is not None
            ]
            for k, ((a, b), s) in enumerate(sorted(marks, key=lambda m: abs(m[0][1] - m[0][0]))):
                yb = top + step * (k + 1)
                ax.plot([a, a, b, b], [yb - step * 0.3, yb, yb, yb - step * 0.3], color=INK, lw=0.6)
                ax.text(
                    (a + b) / 2,
                    yb + step * 0.02,
                    s,
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    color=INK,
                )
            ax.set_ylim(min(0.0, float(y.min()) - 0.2), top + step * (len(marks) + 1.6))
            ax.set_xlim(-0.6, 2.6)
            ax.set_xticks(range(3))
            ax.set_xticklabels(ticks, fontsize=6.5)
            ax.tick_params(axis="x", length=0)
        for ax in axes[len(panel) :]:
            ax.set_axis_off()
        handles = [
            Line2D([], [], ls="", marker=GROUP_STYLE[g][1], ms=5, color=GROUP_STYLE[g][0], label=g)
            for g in order
        ]
        handles.append(Line2D([], [], color=INK, lw=1.3, label="median and interquartile range"))
        fig.legend(
            handles=handles, loc="upper left", ncol=4, bbox_to_anchor=(0.01, 1.0), fontsize=7.5
        )
        fig.text(
            0.99,
            0.995,
            "* FDR below 0.05, ** below 0.01, *** below 0.001, for comparisons meeting the dysregulation rule",
            ha="right",
            va="top",
            fontsize=7,
            color=MUTED,
        )
        fig.supylabel("log$_2$(DESeq2-normalised count + 1)", fontsize=8.5, color=INK)
        fig.tight_layout(rect=(0, 0, 1, 0.975), h_pad=1.2, w_pad=0.8)
        _save(fig, path)


def batch_pca(
    scores: pd.DataFrame, groups: pd.Series, batch: pd.Series, components: pd.DataFrame, path: Path
) -> None:
    """First two principal components, coloured by batch and by study group."""
    if scores.shape[1] < 2:
        return
    var = components["variance_explained_percent"]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
        pad = 0.08 * float(np.ptp(scores[["PC1", "PC2"]].to_numpy()))
        xlim = (float(scores["PC1"].min()) - pad, float(scores["PC1"].max()) + pad)
        ylim = (float(scores["PC2"].min()) - pad, float(scores["PC2"].max()) + pad)
        for ax, labels, title in (
            (axes[0], batch, "By sequencing batch"),
            (axes[1], groups, "By study group"),
        ):
            values = labels.to_numpy()
            for k, level in enumerate(sorted(pd.unique(values), key=str)):
                sub = scores[values == level]
                if labels is groups:
                    colour, marker, filled = _group_style(str(level), k)
                else:
                    colour = BATCH_CYCLE[k % len(BATCH_CYCLE)]
                    marker, filled = MARKERS[k % len(MARKERS)], True
                ax.scatter(
                    sub["PC1"],
                    sub["PC2"],
                    s=16,
                    marker=marker,
                    facecolors=colour if filled else "white",
                    edgecolors="white" if filled else colour,
                    linewidths=0.3 if filled else 0.9,
                    alpha=0.9,
                    label=str(level),
                )
            ax.set_title(title, loc="left", color=INK)
            ax.set_xlabel(f"PC1 ({var.iloc[0]:.1f}%)")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_aspect("equal", adjustable="box")
            ax.axhline(0, lw=0.4, color="#E3E7EE", zorder=0)
            ax.axvline(0, lw=0.4, color="#E3E7EE", zorder=0)
            ax.legend(
                loc="upper center",
                bbox_to_anchor=(0.5, -0.21),
                ncol=4,
                fontsize=6.5,
                handletextpad=0.2,
            )
        axes[0].set_ylabel(f"PC2 ({var.iloc[1]:.1f}%)")
        fig.tight_layout()
        _save(fig, path)
