"""Depth-safe Newick handling and off-screen tree rendering."""

from io import BytesIO, StringIO
import re
from textwrap import wrap

from Bio import Phylo


def tree_nodes(tree):
    """Preorder traversal without BioPython's recursive find_clades helpers."""
    stack = [tree.root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.clades))


def annotate_newick(newick, annotations):
    """Annotate and serialize using an explicit stack, never deepcopy/recursion."""
    tree = Phylo.read(StringIO(newick), "newick")
    for node in tree_nodes(tree):
        if not node.clades and node.name in annotations:
            node.name += "_" + "_".join(str(annotations[node.name]).split())
    parts = []
    stack = [tree.root]
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            parts.append(item)
            continue
        label = item.name or ""
        if re.search(r"[\s\[\](),:;'\"]", label):
            label = "'" + label.replace("\\", "\\\\").replace("'", "\\'") + "'"
        suffix = label
        if item.branch_length is not None:
            suffix += ":" + format(item.branch_length, ".10g")
        if item.clades:
            parts.append("(")
            stack.append(")" + suffix)
            for i in range(len(item.clades) - 1, -1, -1):
                stack.append(item.clades[i])
                if i:
                    stack.append(",")
        else:
            parts.append(suffix)
    return "".join(parts) + ";\n"


def tree_layout(tip_count, font_size=10, spacing=1.0):
    height = max(2.4 / max(tip_count, 1), font_size * 1.8 / 72 * spacing)
    return {"height": height, "width": 7.5, "leaf_label_size": font_size}


def render_tree(newick, query_name=None, font_size=10, spacing=1.0):
    # Import plotting dependencies only when rendering is requested. Use the
    # Agg canvas directly so no process-global backend/recursion changes occur.
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.collections import LineCollection
    from matplotlib.figure import Figure

    tree = Phylo.read(StringIO(newick), "newick")  # BioPython's parser is iterative.
    nodes = list(tree_nodes(tree))
    leaves = [node for node in nodes if not node.clades]
    tips = len(leaves)
    labels = {node: "\n".join(wrap(node.name or "", width=64)) for node in leaves}
    layout = tree_layout(tips, font_size, spacing)
    layout["height"] *= max(label.count("\n") + 1 for label in labels.values())

    x = {tree.root: tree.root.branch_length or 0.0}
    for node in nodes:
        for child in node.clades:
            x[child] = x[node] + (child.branch_length or 0.0)
    # Zero-distance trees still need a visible topology. Unit edges are used
    # only for display, and the Newick distances remain untouched.
    unit_edges = max(x.values()) == min(x.values())
    if unit_edges:
        x[tree.root] = 0
        for node in nodes:
            for child in node.clades:
                x[child] = x[node] + 1
    y = {leaf: i for i, leaf in enumerate(leaves)}
    segments = []
    for node in reversed(nodes):
        if node.clades:
            y[node] = (y[node.clades[0]] + y[node.clades[-1]]) / 2
            segments.append(
                [(x[node], y[node.clades[0]]), (x[node], y[node.clades[-1]])]
            )
            for child in node.clades:
                segments.append([(x[node], y[child]), (x[child], y[child])])

    height = tips * layout["height"]
    dpi = min(150, 12000 / (height + 1))
    fig = Figure(figsize=(layout["width"], height), dpi=dpi)
    FigureCanvasAgg(fig)
    try:
        ax = fig.add_subplot(111)
        ax.add_collection(LineCollection(segments, colors="#333333", linewidths=0.8))
        span = max(x.values()) - min(x.values()) or 1
        query_names = {query_name, query_name.replace(" ", "_") if query_name else None}
        for leaf in leaves:
            ax.text(
                x[leaf] + span * 0.01,
                y[leaf],
                labels[leaf],
                fontsize=font_size,
                va="center",
                color="#d94841" if leaf.name in query_names else "#222222",
            )
        ax.set_xlim(min(x.values()) - span * 0.02, max(x.values()) + span * 0.02)
        ax.set_ylim(tips - 0.5, -0.5)
        ax.set_axis_off()
        if unit_edges:
            ax.set_title(
                "Zero branch distances · unit edges shown for readability", fontsize=9
            )
        png, svg = BytesIO(), BytesIO()
        fig.savefig(png, format="png", dpi=dpi, bbox_inches="tight", pad_inches=0.2)
        fig.savefig(svg, format="svg", bbox_inches="tight", pad_inches=0.2)
        return png.getvalue(), svg.getvalue(), tips
    finally:
        fig.clear()
