from matplotlib import cbook, cm, colors, rcParams
from cycler import cycler


# default matplotlib 2.0 palette slightly modified.
vega_10 = list(map(colors.to_hex, cm.tab10.colors))
vega_20 = list(map(colors.to_hex, cm.tab20.colors))
vega_20 = [
    *vega_20[0:14:2],
    *vega_20[16::2],
    *vega_20[1:15:2],
    *vega_20[17::2],
    "#ad494a",
    "#8c6d31",
]
vega_20[2], vega_20[4], vega_20[7] = (
    "#279e68",
    "#aa40fc",
    "#b5bd61",
)  # green, purple, khaki


def set_figure_params(fontsize=12, color_map=None, frameon=None, figsize=None, dpi=100, dpi_save=300):
    # dpi options
    if figsize is not None:
        rcParams["figure.figsize"] = figsize
    if dpi is not None:
        rcParams["figure.dpi"] = dpi
    if dpi_save is not None:
        rcParams["savefig.dpi"] = dpi_save
    rcParams['pdf.fonttype'] = 42

    rcParams["figure.subplot.left"] = 0.18
    rcParams["figure.subplot.right"] = 0.96
    rcParams["figure.subplot.bottom"] = 0.15
    rcParams["figure.subplot.top"] = 0.91

    rcParams["lines.linewidth"] = 1.5  # the line width of the frame
    rcParams["lines.markersize"] = 6
    rcParams["lines.markeredgewidth"] = 1

    # font
    rcParams["font.sans-serif"] = [
        "Arial",
        "Helvetica",
        "DejaVu Sans",
        "Bitstream Vera Sans",
        "sans-serif",
    ]

    fontsize = fontsize
    labelsize = 0.92 * fontsize

    rcParams["font.size"] = fontsize
    rcParams["legend.fontsize"] = labelsize
    rcParams["axes.titlesize"] = fontsize
    rcParams["axes.labelsize"] = fontsize

    # legend
    rcParams["legend.numpoints"] = 1
    rcParams["legend.scatterpoints"] = 1
    rcParams["legend.handlelength"] = 0.5
    rcParams["legend.handletextpad"] = 0.4

    # color cycle
    rcParams["axes.prop_cycle"] = cycler(color=vega_20)

    # lines
    rcParams["axes.linewidth"] = 0.8
    rcParams["axes.edgecolor"] = "black"
    rcParams["axes.facecolor"] = "white"

    # ticks
    rcParams["xtick.color"] = "k"
    rcParams["ytick.color"] = "k"
    rcParams["xtick.labelsize"] = fontsize
    rcParams["ytick.labelsize"] = fontsize

    # axes grid
    rcParams["axes.grid"] = True
    rcParams["grid.color"] = ".8"
    rcParams["axes.spines.bottom"] = "on"
    rcParams["axes.spines.top"] = "off"
    rcParams["axes.spines.left"] = "on"
    rcParams["axes.spines.right"] = "off"
    rcParams['axes.axisbelow'] = True

    # color map
    rcParams["image.cmap"] = "Spectral_r" if color_map is None else color_map


    # frame
    frameon = True if frameon is None else frameon
    global _frameon
    _frameon = frameon