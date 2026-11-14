"""Drawing a 384-well source plate: one rounded cell per well, text inside, rings and marks.

Colour rules (dataviz skill, reference palette):
  * a plate map is an any-two-cells-can-touch form, so at most THREE categorical hues are
    used in one figure (validated all-pairs: worst CVD dE 9.2, normal-vision 24.0); more
    identities than that are carried by the text label inside the well, never by new hues;
  * magnitudes use one blue ramp, light -> dark, starting no lighter than step 250 so the
    lightest step still separates from the surface;
  * text is ink, never series colour: white inside dark fills, near-black inside light ones.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

ROWS = "ABCDEFGHIJKLMNOP"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
HAIRLINE = "#e1e0d9"
NEUTRAL = "#f0efec"          # "nothing to report" fill
OTHER = "#52514e"            # folded "other" identity: an ink step, not a hue
SLOT = ["#2a78d6", "#eb6834", "#1baf7a"]
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]   # blue 250..650
DARK_FILLS = {OTHER, SLOT[0], SLOT[1], RAMP[2], RAMP[3], RAMP[4]}

GAP = 0.10                   # surface gap between cells
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                     "pdf.fonttype": 42})


def text_ink(fill):
    return "#ffffff" if fill in DARK_FILLS else INK


def new_plate_figure(width=17.0, height=11.6):
    fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
    ax = fig.add_axes([0.035, 0.165, 0.955, 0.682])
    ax.set_facecolor(SURFACE)
    return fig, ax


def draw_plate(ax, cells):
    """cells: well -> dict(fill, edge, edge_w, hatch, top, mid, mid_size, mark, mark_color,
    dashed). Missing wells are drawn empty."""
    for r, row in enumerate(ROWS):
        for c in range(1, 25):
            w = f"{row}{c}"
            d = cells.get(w, {})
            fill = d.get("fill", SURFACE)
            x, y = c - 1 + GAP / 2, r + GAP / 2
            box = FancyBboxPatch((x, y), 1 - GAP, 1 - GAP,
                                 boxstyle="round,pad=0,rounding_size=0.14",
                                 facecolor=fill, edgecolor=d.get("edge", HAIRLINE),
                                 linewidth=d.get("edge_w", 0.6), hatch=d.get("hatch"),
                                 linestyle="--" if d.get("dashed") else "-")
            ax.add_patch(box)
            ink = d.get("ink", text_ink(fill))
            if d.get("top"):
                ax.text(c - 0.5, r + 0.30, d["top"], ha="center", va="center", fontsize=5.6,
                        color=ink if fill != SURFACE else INK_2)
            if d.get("mid"):
                ax.text(c - 0.5, r + 0.66, d["mid"], ha="center", va="center",
                        fontsize=d.get("mid_size", 6.4), fontweight="bold", color=ink)
            if d.get("mark"):
                ax.plot(c - 0.17, r + 0.20, marker=d["mark"], markersize=4.2,
                        markerfacecolor=d.get("mark_color", ink), markeredgecolor=d.get("mark_edge", ink),
                        markeredgewidth=0.7, linestyle="none")
    ax.set_xlim(-0.05, 24.05); ax.set_ylim(16.05, -0.05)
    ax.set_aspect("equal")
    ax.set_xticks([c - 0.5 for c in range(1, 25)], [str(c) for c in range(1, 25)])
    ax.set_yticks([r + 0.5 for r in range(16)], list(ROWS))
    ax.xaxis.tick_top()
    ax.tick_params(length=0, labelsize=8, colors=MUTED, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)


def connect(ax, wells, color=INK_2, lw=0.9):
    """Hairline connectors between the wells of one set. They start and stop at the cell
    edge and sit under the text, so they never cross a label inside an endpoint well."""
    import numpy as np
    pts = sorted(((int(w[1:]) - 0.5), ROWS.index(w[0]) + 0.5) for w in wells)
    for a, b in zip(pts, pts[1:]):
        v = np.subtract(b, a); L = np.hypot(*v); u = v / L
        trim = 0.47 / max(abs(u[0]), abs(u[1]))          # leave the square cell, any angle
        p0, p1 = np.add(a, u * trim), np.subtract(b, u * trim)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, alpha=0.6, zorder=2.5,
                solid_capstyle="round")


def header(fig, title, subtitle):
    fig.text(0.035, 0.965, title, fontsize=15, fontweight="bold", color=INK, va="top")
    fig.text(0.035, 0.925, subtitle, fontsize=9.6, color=INK_2, va="top", linespacing=1.45)


def legend(fig, items, ncol=4, y=0.152):
    """items: list of (kind, label, style) where kind is 'patch' or 'mark'."""
    handles = []
    for kind, label, st in items:
        if kind == "patch":
            handles.append(Patch(facecolor=st.get("fill", SURFACE), edgecolor=st.get("edge", HAIRLINE),
                                 linewidth=st.get("edge_w", 0.6), hatch=st.get("hatch"),
                                 linestyle="--" if st.get("dashed") else "-", label=label))
        else:
            handles.append(Line2D([], [], marker=st["mark"], linestyle="none",
                                  markersize=st.get("size", 6),
                                  markerfacecolor=st.get("mark_color", INK),
                                  markeredgecolor=st.get("mark_edge", INK), label=label))
    lg = fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.035, y), ncol=ncol,
                    frameon=False, fontsize=8.6, handlelength=1.6, handleheight=1.2,
                    columnspacing=2.2, labelcolor=INK)
    return lg


def footnote(fig, text, y=0.008):
    fig.text(0.035, y, text, fontsize=8.0, color=MUTED, va="bottom", linespacing=1.4)


def save(fig, out, stem):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{stem}.png", dpi=200, facecolor=SURFACE)
    fig.savefig(out / f"{stem}.pdf", facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out / stem}.png / .pdf")
