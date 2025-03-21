#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle

# global variables coming from the CLI
# CLI_INSERT_CONFIG_BEGIN
# add some placeholders meanwhile
title = ""
xlabel = ""
ylabel = ""
xlog = False
ylog = False
scales = 1
plot_panels = {
    "int": False,
    "abs": False,
    "rel_ewonoff": False,
    "abs_pdfs": False,
    "ratio_pdf": False,
    "double_ratio_pdf": False,
    "rel_pdfunc": False,
    "rel_pdfpull": False,
}
output = ""
data = {}
metadata = {}
# CLI_INSERT_CONFIG_END

# color cycler for different PDF results
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# color for the first PDF result with QCD-only predictions
colors0_qcd = "red"

# stylesheet for plot
stylesheet = {
    "axes.axisbelow": True,
    "axes.grid": True,
    "axes.labelsize": "small",
    "figure.constrained_layout.hspace": 0.0,
    "figure.constrained_layout.use": True,
    "figure.constrained_layout.wspace": 0.0,
    "font.family": "serif",
    "font.size": 14.0,
    "grid.linestyle": "dotted",
    "legend.borderpad": 0.0,
    "legend.fontsize": "x-small",
    "legend.frameon": False,
    "pdf.compression": 0,
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{siunitx}\usepackage{lmodern}\usepackage[T1]{fontenc}",
    "xtick.bottom": True,
    "xtick.top": True,
    "xtick.direction": "in",
    "xtick.major.width": 0.5,
    "xtick.minor.bottom": True,
    "xtick.minor.top": True,
    "xtick.minor.width": 0.5,
    "ytick.direction": "in",
    "ytick.left": True,
    "ytick.right": True,
    "ytick.major.width": 0.5,
    "ytick.minor.visible": True,
    "ytick.minor.width": 0.5,
}

# panel plot labels
ylabel_ratio_pdf = r"Ratio to {{central_pdf}}"
ylabel_double_ratio_pdf = r"Ratio to NLO"
ylabel_rel_ewonoff = r"NLO EW on/off [\si{\percent}]"
ylabel_rel_pdfunc = r"PDF uncertainty [\si{\percent}]"
ylabel_rel_pdfpull = r"Pull [$\sigma$]"

label_rel_ewonoff_qcd = r"NLO QCD"
label_rel_ewonoff_ew = r"NLO QCD+EW"
label_rel_ewonoff_scale_unc = f"{scales}-p. scale var."
label_rel_ewonoff_pdf_unc = r"PDF uncertainty"


# linestyle for the channel breakdown shown in the panel `plot_abs_pdfs`. If the array
# is empty, no channel breakdown will be shown, otherwise the most important channels,
# as many as linestyles are given. See also
# https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
channel_breakdown_linestyles = []


def main(active_panels):
    """Build a plot figure with various panels."""
    # Find the active panels
    panels = [
        PANEL_FNC_MAP[panel] for panel, enabled in active_panels.items() if enabled
    ]

    # prepare the figure
    mpl.rcParams.update(stylesheet)
    if len(panels) == 1:
        plt.rc("figure", figsize=(4.2, 2.6))
    else:
        plt.rc("figure", figsize=(6.4, 2.4 * len(panels)))

    # Plot all data
    for index, kwargs in enumerate(data):
        figure, axes = plt.subplots(len(panels), 1, sharex=True, squeeze=False)

        if len(kwargs["x"]) > 2 and xlog:
            axes[0, 0].set_xscale("log")

        axes[0, 0].set_title(title)
        axes[-1, 0].set_xlabel(xlabel)

        for plot, axis in zip(panels, axes[:, 0]):
            plot(axis, **kwargs)

        if len(data) == 1:
            figure.savefig(f"{output}.pdf")
        else:
            figure.savefig(f"{output}-{index}.pdf")
        plt.close(figure)


def percent_diff(a, b):
    return (a / b - 1.0) * 100.0


def set_ylim(axis, save, load, filename):
    # extract the y limits *not* considering margins
    margins = axis.margins()
    axis.margins(y=0.0)
    ymin, ymax = axis.get_ylim()
    axis.margins(y=margins[1])

    inc = 1.0

    if (ymax - ymin) > 100.0:
        ymin = -50.0
        ymax = 50.0
        inc = 25.0
    elif (ymax - ymin) > 30.5:
        inc = 10.0
    elif (ymax - ymin) > 20.5:
        inc = 5.0
    elif (ymax - ymin) > 10.5:
        inc = 2.0
    elif (ymax - ymin) < 3.0:
        inc = 0.5

    ymin = np.floor(ymin / inc) * inc
    ymax = np.ceil(ymax / inc) * inc

    if save:
        with open(filename, "wb") as f:
            pickle.dump([ymin, ymax, inc], f)

    if load:
        resave = False

        with open(filename, "rb") as f:
            [saved_ymin, saved_ymax, saved_inc] = pickle.load(f)

        if saved_ymin < ymin:
            ymin = saved_ymin
            resave = True

        if saved_ymax > ymax:
            ymax = saved_ymax
            resave = True

        if saved_inc > inc:
            inc = saved_inc
            resave = True

        if resave:
            with open(filename, "wb") as f:
                pickle.dump([ymin, ymax, inc], f)

    axis.set_yticks(np.arange(ymin, ymax + inc, inc))
    space = 0.05 * (ymax - ymin)
    axis.set_ylim((ymin - space, ymax + space))


def plot_int(axis, /, pdf_results, **_kwargs):
    xmin = np.array([])
    xmax = np.array([])
    x = np.array([])
    y = np.array([])

    for index, (label, ycentral, ymin, ymax) in enumerate(pdf_results):
        x = np.append(x, ycentral[:-1])
        xmin = np.append(xmin, ymin[:-1])
        xmax = np.append(xmax, ymax[:-1])
        y = np.append(y, label)

        # draw one- and two-sigma bands
        if label == "CENTRAL-PDF":
            axis.axvspan(
                xmin[-1], xmax[-1], alpha=0.3, color=colors[index], linewidth=0
            )
            # TODO: this is only correct for MC PDF uncertainties
            axis.axvspan(
                x[-1] - 2.0 * (x[-1] - xmin[-1]),
                x[-1] + 2.0 * (xmax[-1] - x[-1]),
                alpha=0.1,
                color=colors[index],
                linewidth=0,
            )

    axis.errorbar(
        x, y, xerr=(x - xmin, xmax - x), fmt=".", capsize=3, markersize=5, linewidth=1.5
    )
    axis.margins(x=0.1, y=0.1)


def plot_abs(axis, /, x, y, ymin, ymax, slice_label="", **_kwargs):
    axis.set_yscale("log" if ylog else "linear")
    axis.step(x, y, colors[0], linewidth=1.0, where="post", label=slice_label)
    axis.fill_between(
        x,
        ymin,
        ymax,
        alpha=0.4,
        color=colors[0],
        linewidth=0.5,
        step="post",
    )
    axis.set_ylabel(ylabel)

    if slice_label != "":
        axis.legend()


def plot_ratio_pdf(axis, /, x, pdf_results, slice_label="", **_kwargs):
    axis.set_ylabel(ylabel_ratio_pdf.format(central_pdf=pdf_results[0][0]))

    for index, i in enumerate(pdf_results):
        label, y, ymin, ymax = i
        y = y / pdf_results[0][1]
        ymin = ymin / pdf_results[0][1]
        ymax = ymax / pdf_results[0][1]

        axis.step(x, y, color=colors[index], linewidth=1.0, where="post")
        axis.fill_between(
            x,
            ymin,
            ymax,
            alpha=0.4,
            color=colors[index],
            label=label,
            linewidth=0.5,
            step="post",
        )

    axis.legend(
        bbox_to_anchor=(0, -0.24, 1, 0.2),
        loc="upper left",
        mode="expand",
        borderaxespad=0,
        ncol=min(4, len(pdf_results)),
    )

    if slice_label != "":
        t = axis.text(
            0.98,
            0.98,
            slice_label,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axis.transAxes,
            fontsize="x-small",
        )
        t.set_bbox(
            {
                "alpha": 0.7,
                "boxstyle": "square, pad=0.0",
                "edgecolor": "white",
                "facecolor": "white",
            }
        )


def plot_double_ratio_pdf(axis, /, x, pdf_results, slice_label="", **_kwargs):
    axis.set_ylabel(ylabel_double_ratio_pdf)

    for index, i in enumerate(pdf_results):
        label, y, ymin, ymax = i
        if index == 0 or index == 2:
            y = y / pdf_results[0][1]
            ymin = ymin / pdf_results[0][1]
            ymax = ymax / pdf_results[0][1]
        else:
            y = y / pdf_results[1][1]
            ymin = ymin / pdf_results[1][1]
            ymax = ymax / pdf_results[1][1]
        axis.step(x, y, color=colors[index], linewidth=1.0, where="post")
        axis.fill_between(
            x,
            ymin,
            ymax,
            alpha=0.4,
            color=colors[index],
            label=label,
            linewidth=0.5,
            step="post",
        )

    axis.legend(
        bbox_to_anchor=(0, -0.24, 1, 0.2),
        loc="upper left",
        mode="expand",
        borderaxespad=0,
        ncol=min(4, len(pdf_results)),
    )

    if slice_label != "":
        t = axis.text(
            0.98,
            0.98,
            slice_label,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axis.transAxes,
            fontsize="x-small",
        )
        t.set_bbox(
            {
                "alpha": 0.7,
                "boxstyle": "square, pad=0.0",
                "edgecolor": "white",
                "facecolor": "white",
            }
        )


def plot_abs_pdfs(axis, /, x, pdf_results, channels, slice_label="", **_kwargs):
    axis.set_yscale("log" if ylog else "linear")
    axis.set_ylabel(ylabel)

    for index, i in enumerate(pdf_results):
        label, y, ymin, ymax = i
        axis.step(x, y, color=colors[index], linewidth=1.0, where="post")
        axis.fill_between(
            x,
            ymin,
            ymax,
            alpha=0.4,
            color=colors[index],
            label=label,
            linewidth=0.5,
            step="post",
        )

    for index, ((label, y), linestyle) in enumerate(
        zip(channels, channel_breakdown_linestyles)
    ):
        axis.step(
            x,
            y,
            color=colors[0],
            label=label,
            linestyle=linestyle,
            linewidth=1.0,
            where="post",
        )

    axis.legend(
        bbox_to_anchor=(0, -0.24, 1, 0.2),
        loc="upper left",
        mode="expand",
        borderaxespad=0,
        ncol=min(4, len(pdf_results) + len(channel_breakdown_linestyles)),
    )

    if slice_label != "":
        t = axis.text(
            0.98,
            0.98,
            slice_label,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axis.transAxes,
            fontsize="x-small",
        )
        t.set_bbox(
            {
                "alpha": 0.7,
                "boxstyle": "square, pad=0.0",
                "edgecolor": "white",
                "facecolor": "white",
            }
        )


def plot_rel_ewonoff(axis, /, x, mid, y, ymin, ymax, qcd_y, pdf_results, **_kwargs):
    y = percent_diff(y, qcd_y)
    qcd_y = percent_diff(qcd_y, qcd_y)
    # qcd_ymin = percent_diff(kwargs["qcd_min"], kwargs["qcd_y"])
    # qcd_ymax = percent_diff(kwargs["qcd_max"], kwargs["qcd_y"])
    ymin = percent_diff(ymin, qcd_y)
    ymax = percent_diff(ymax, qcd_y)
    pdf_min = abs(percent_diff(pdf_results[0][2], pdf_results[0][1]))[:-1]
    pdf_max = abs(percent_diff(pdf_results[0][3], pdf_results[0][1]))[:-1]

    axis.step(
        x, qcd_y, colors0_qcd, label=label_rel_ewonoff_qcd, linewidth=1.0, where="post"
    )
    # axis.fill_between(x, qcd_ymin, qcd_ymax, alpha=0.4, color=colors0_qcd, label=label_rel_ewonoff_scale_unc, linewidth=0.5, step="post")
    axis.step(x, y, colors[0], label=label_rel_ewonoff_ew, linewidth=1.0, where="post")
    axis.fill_between(
        x,
        ymin,
        ymax,
        alpha=0.4,
        color=colors[0],
        label=label_rel_ewonoff_scale_unc,
        linewidth=0.5,
        step="post",
    )
    axis.errorbar(
        mid,
        y[:-1],
        yerr=(pdf_min, pdf_max),
        color=colors[0],
        label=label_rel_ewonoff_pdf_unc,
        fmt=".",
        capsize=1,
        markersize=0,
        linewidth=1,
    )
    axis.set_ylabel(ylabel_rel_ewonoff)
    axis.legend(
        bbox_to_anchor=(0, 1.03, 1, 0.2),
        loc="lower left",
        mode="expand",
        borderaxespad=0,
        ncol=4,
    )


def plot_rel_pdfunc(axis, /, x, pdf_results, **_kwargs):
    for index, i in enumerate(pdf_results):
        label, y, ymin, ymax = i
        ymin = percent_diff(ymin, y)
        ymax = percent_diff(ymax, y)
        axis.step(x, ymax, color=colors[index], label=label, linewidth=1, where="post")
        axis.step(x, ymin, color=colors[index], linewidth=1, where="post")

    axis.set_ylabel(ylabel_rel_pdfunc)

    set_ylim(axis, False, False, "rel_pdfunc")


def plot_rel_pdfpull(axis, /, x, y, pdf_results, **_kwargs):
    central_y = pdf_results[0][1]
    central_ymin = pdf_results[0][2]
    central_ymax = pdf_results[0][3]

    for index, i in enumerate(pdf_results):
        label, y, ymin, ymax = i
        diff = y - central_y
        yerr = np.where(diff > 0.0, y - ymin, ymax - y)
        cerr = np.where(diff > 0.0, central_ymax - central_y, central_y - central_ymin)
        pull = diff / np.sqrt(np.power(yerr, 2) + np.power(cerr, 2))

        axis.step(
            x,
            pull,
            color=colors[index],
            label=label,
            linewidth=1,
            where="post",
            zorder=2 * index + 1,
        )

    axis.legend(
        bbox_to_anchor=(0, 1.03, 1, 0.2),
        loc="lower left",
        mode="expand",
        borderaxespad=0,
        ncol=min(4, len(pdf_results)),
    )  # rel_pdfpull
    axis.set_ylabel(ylabel_rel_pdfpull)

    set_ylim(axis, False, False, "rel_pdfpull")


PANEL_FNC_MAP = {
    "int": plot_int,
    "abs": plot_abs,
    "rel_ewonoff": plot_rel_ewonoff,
    "abs_pdfs": plot_abs_pdfs,
    "ratio_pdf": plot_ratio_pdf,
    "double_ratio_pdf": plot_double_ratio_pdf,
    "rel_pdfunc": plot_rel_pdfunc,
    "rel_pdfpull": plot_rel_pdfpull,
}


# CLI data variables
# CLI_INSERT_DATA
# end CLI data variables


if __name__ == "__main__":
    main(plot_panels)
