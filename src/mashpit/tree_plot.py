"""Shared adaptive tree rendering for CLI exports and the Streamlit report."""
from io import BytesIO, StringIO
from textwrap import wrap

from Bio import Phylo
import matplotlib

# Both the CLI and Streamlit render off-screen, including on macOS workers.
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from phytreeviz import TreeViz


def tree_layout(tip_count, font_size=10, spacing=1.0):
    # TreeViz height is inches PER TIP. Match line spacing to point size,
    # with enough room for small trees and no compression of large trees.
    height = max(2.4 / max(tip_count, 1), font_size * 1.8 / 72 * spacing)
    return {"height": height, "width": 7.5, "leaf_label_size": font_size}


def render_tree(newick, query_name=None, font_size=10, spacing=1.0):
    tree = Phylo.read(StringIO(newick), "newick")
    tips = len(tree.get_terminals())
    layout = tree_layout(tips, font_size, spacing)
    label_lines = max(
        (len(wrap(tip.name or "", width=64)) for tip in tree.get_terminals()), default=1
    )
    layout["height"] *= max(1, label_lines)
    tv = TreeViz(tree, **layout)
    names = {tip.name for tip in tree.get_terminals()}
    if query_name in names:
        tv.set_node_label_props(query_name, color="#d94841")
    elif query_name and query_name.replace(" ", "_") in names:
        tv.set_node_label_props(query_name.replace(" ", "_"), color="#d94841")
    # Bound raster height for large result sets; SVG retains readable vector text.
    dpi = min(150, max(30, int(12000 / (tips * layout["height"] + 1))))
    fig = tv.plotfig(dpi=dpi)
    for ax in fig.axes:
        for label in ax.texts:
            label.set_text("\n".join(wrap(label.get_text(), width=64)))
    try:
        png, svg = BytesIO(), BytesIO()
        fig.savefig(png, format="png", dpi=dpi, bbox_inches="tight", pad_inches=.2)
        fig.savefig(svg, format="svg", bbox_inches="tight", pad_inches=.2)
        return png.getvalue(), svg.getvalue(), tips
    finally:
        plt.close(fig)
