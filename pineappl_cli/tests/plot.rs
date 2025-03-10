#![allow(missing_docs)]

use assert_cmd::Command;
use predicates::str;
use std::num::NonZeroUsize;
use std::thread;

const HELP_STR: &str = "Creates a matplotlib script plotting the contents of the grid

Usage: pineappl plot [OPTIONS] <INPUT> <CONV_FUNS>...

Arguments:
  <INPUT>         Path to the input grid
  <CONV_FUNS>...  LHAPDF id(s) or name of the PDF set(s)

Options:
      --conv-fun-uncert-from <IDX>  Choose for which convolution function the uncertainty should be calculated [default: 0]
  -s, --scales <SCALES>             Set the number of scale variations [default: 7] [possible values: 1, 3, 7, 9]
      --asymmetry                   Plot the asymmetry
      --threads <THREADS>           Number of threads to utilize [default: {}]
      --no-conv-fun-unc             Disable the (time-consuming) calculation of PDF uncertainties
  -h, --help                        Print help
";

const DEFAULT_STR: &str = r#"#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle

# global variables coming from the CLI
title = r"LHCb differential W-boson production cross section at 7 TeV"
xlabel = r"$\eta_{\bar{\ell}}$"
ylabel = r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\eta_{\bar{\ell}}}$ [\si{pb}]"
xlog = False
ylog = False
scales = 7
plot_panels = {
    "int": False,
    "abs": True,
    "rel_ewonoff": True,
    "abs_pdfs": True,
    "ratio_pdf": True,
    "double_ratio_pdf": True,
    "rel_pdfunc": True,
    "rel_pdfpull": True,
}
output = r"LHCB_WP_7TEV_opt"

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
data = [
    {
        "slice_label"    : r"",
        "x"        : np.array([2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4, 4.5]),
        "mid"      : np.array([2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.75, 4.25]),
        "pdf_results" : [
            (
                    r"NNPDF31\_nlo\_as\_0118\_luxqed",
                    np.array([7.5461655e2, 6.9027941e2, 6.0022595e2, 4.8548211e2, 3.6191001e2, 2.4582640e2, 1.1584074e2, 2.7504644e1, 2.7504644e1]),
                    np.array([7.4602445e2, 6.8228693e2, 5.9312383e2, 4.7957251e2, 3.5732461e2, 2.4250580e2, 1.1409722e2, 2.6743047e1, 2.6743047e1]),
                    np.array([7.6320866e2, 6.9827189e2, 6.0732807e2, 4.9139170e2, 3.6649540e2, 2.4914699e2, 1.1758425e2, 2.8266240e1, 2.8266240e1]),
                ),
            (
                    r"NNPDF4.0",
                    np.array([7.8845642e2, 7.2060104e2, 6.2525179e2, 5.0384928e2, 3.7399422e2, 2.5300320e2, 1.1909100e2, 2.9002507e1, 2.9002507e1]),
                    np.array([7.8450378e2, 7.1690933e2, 6.2191897e2, 5.0097415e2, 3.7162070e2, 2.5109709e2, 1.1772232e2, 2.8044400e1, 2.8044400e1]),
                    np.array([7.9240905e2, 7.2429275e2, 6.2858461e2, 5.0672441e2, 3.7636773e2, 2.5490931e2, 1.2045969e2, 2.9960613e1, 2.9960613e1]),
                ),
        ],
        "qcd_y"    : np.array([7.6246034e2, 6.9684577e2, 6.0548681e2, 4.8928139e2, 3.6454175e2, 2.4754316e2, 1.1667878e2, 2.7737493e1, 2.7737493e1]),
        "qcd_min"  : np.array([7.3365413e2, 6.7075857e2, 5.8286230e2, 4.7122879e2, 3.5126111e2, 2.3864984e2, 1.1258836e2, 2.6810991e1, 2.6810991e1]),
        "qcd_max"  : np.array([7.8308860e2, 7.1589472e2, 6.2240294e2, 5.0304405e2, 3.7486018e2, 2.5457732e2, 1.1997317e2, 2.8492180e1, 2.8492180e1]),
        "y"        : np.array([7.5459110e2, 6.9028342e2, 6.0025198e2, 4.8552235e2, 3.6195456e2, 2.4586691e2, 1.1586851e2, 2.7517266e1, 2.7517266e1]),
        "ymin"     : np.array([7.2595107e2, 6.6417441e2, 5.7747295e2, 4.6723687e2, 3.4843600e2, 2.3677625e2, 1.1166942e2, 2.6562471e1, 2.6562471e1]),
        "ymax"     : np.array([7.7529764e2, 7.0957480e2, 6.1750966e2, 4.9966237e2, 3.7261436e2, 2.5316566e2, 1.1930157e2, 2.8306905e1, 2.8306905e1]),
        "channels" : [
            (r"$\mathrm{u}\bar{\mathrm{d}} + \mathrm{c}\bar{\mathrm{s}}$", np.array([8.4002759e2, 7.7448295e2, 6.7891182e2, 5.5341626e2, 4.1562095e2, 2.8427837e2, 1.3470473e2, 3.1886258e1, 3.1886258e1])),
            (r"$\mathrm{u}\mathrm{g} + \mathrm{c}\mathrm{g}$", np.array([-6.0727462e1, -6.1109036e1, -5.7385834e1, -4.9385114e1, -3.8287410e1, -2.6578788e1, -1.2142190e1, -2.3686722e0, -2.3686722e0])),
            (r"$\mathrm{g}\bar{\mathrm{s}} + \mathrm{g}\bar{\mathrm{d}}$", np.array([-2.4969360e1, -2.3319483e1, -2.1436419e1, -1.8639887e1, -1.5462782e1, -1.1889878e1, -6.7199873e0, -2.0056686e0, -2.0056686e0])),
            (r"$\mathrm{u}\gamma + \mathrm{c}\gamma$", np.array([1.7176328e-1, 1.4518685e-1, 1.1534278e-1, 7.2943823e-2, 4.9352954e-2, 3.8564621e-2, 1.2734974e-2, 3.4154203e-3, 3.4154203e-3])),
            (r"$\gamma\bar{\mathrm{s}} + \gamma\bar{\mathrm{d}}$", np.array([8.8565923e-2, 8.3802762e-2, 4.7074109e-2, 5.8147927e-2, 3.4452663e-2, 1.8643688e-2, 1.3223117e-2, 1.9334685e-3, 1.9334685e-3]))
        ],
    },
]
metadata = {
    "arxiv": r"1505.07024",
    "description": r"LHCb differential W-boson production cross section at 7 TeV",
    "hepdata": r"10.17182/hepdata.2114.v1/t4",
    "initial_state_1": r"2212",
    "initial_state_2": r"2212",
    "lumi_id_types": r"pdg_mc_ids",
    "mg5amc_repo": r"http://bazaar.launchpad.net/~maddevelopers/mg5amcnlo/3.1.2/",
    "mg5amc_revno": r"983",
    "nnpdf_id": r"LHCBWZMU7TEV",
    "pineappl_gitversion": r"v0.4.1-114-gdce19e0",
    "runcard_gitversion": r"82de4ad",
    "x1_label": r"etal",
    "x1_label_tex": r"$\eta_{\bar{\ell}}$",
    "x1_unit": r"",
    "y_label": r"dsig/detal",
    "y_label_tex": r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\eta_{\bar{\ell}}}$",
    "y_unit": r"\pico\barn",
}
# end CLI data variables


if __name__ == "__main__":
    main(plot_panels)
"#;

const DRELL_YAN_AFB_STR: &str = r#"#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle

# global variables coming from the CLI
title = r""
xlabel = r"$\cos \theta^*$"
ylabel = r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\cos \theta^*}$ [\si{pb}]"
xlog = False
ylog = False
scales = 7
plot_panels = {
    "int": False,
    "abs": True,
    "rel_ewonoff": True,
    "abs_pdfs": False,
    "ratio_pdf": False,
    "double_ratio_pdf": False,
    "rel_pdfunc": False,
    "rel_pdfpull": False,
}
output = r"CMS_DY_14TEV_MLL_6000_COSTH"

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
data = [
    {
        "slice_label"    : r"",
        "x"        : np.array([0, 0.040000000000000036, 0.08000000000000007, 0.1200000000000001, 0.15999999999999992, 0.19999999999999996, 0.24, 0.28, 0.32000000000000006, 0.3600000000000001, 0.40000000000000013, 0.43999999999999995, 0.48, 0.52, 0.56, 0.6000000000000001, 0.6400000000000001, 0.6799999999999999, 0.72, 0.76, 0.8, 0.8400000000000001, 0.8800000000000001, 0.9199999999999999, 0.96, 1]),
        "mid"      : np.array([0.020000000000000018, 0.06000000000000005, 0.10000000000000009, 0.14, 0.17999999999999994, 0.21999999999999997, 0.26, 0.30000000000000004, 0.3400000000000001, 0.3800000000000001, 0.42000000000000004, 0.45999999999999996, 0.5, 0.54, 0.5800000000000001, 0.6200000000000001, 0.66, 0.7, 0.74, 0.78, 0.8200000000000001, 0.8600000000000001, 0.9, 0.94, 0.98]),
        "pdf_results" : [
            (
                    r"NNPDF40\_nnlo\_as\_01180",
                    np.array([-8.9073097e-3, -2.5258011e-2, -4.2074993e-2, -5.7119354e-2, -7.2072246e-2, -8.8769233e-2, -1.0168073e-1, -1.1291620e-1, -1.2694004e-1, -1.3981641e-1, -1.4822989e-1, -1.5828910e-1, -1.6656074e-1, -1.7672691e-1, -1.8106409e-1, -1.8604958e-1, -1.9264779e-1, -1.9711216e-1, -1.9947876e-1, -2.0333702e-1, -2.0602952e-1, -2.0751500e-1, -2.0920588e-1, -2.0522087e-1, -1.5755961e-1, -1.5755961e-1]),
                    np.array([-1.3416731e-2, -3.8485931e-2, -6.4022917e-2, -8.7525834e-2, -1.1059996e-1, -1.3506164e-1, -1.5518895e-1, -1.7370393e-1, -1.9408882e-1, -2.1338899e-1, -2.2702283e-1, -2.4212719e-1, -2.5490585e-1, -2.6928706e-1, -2.7696185e-1, -2.8511886e-1, -2.9435668e-1, -3.0107075e-1, -3.0528544e-1, -3.1066117e-1, -3.1459806e-1, -3.1687573e-1, -3.1925830e-1, -3.1465864e-1, -2.4456638e-1, -2.4456638e-1]),
                    np.array([-4.3978881e-3, -1.2030092e-2, -2.0127068e-2, -2.6712874e-2, -3.3544535e-2, -4.2476821e-2, -4.8172510e-2, -5.2128475e-2, -5.9791257e-2, -6.6243833e-2, -6.9436948e-2, -7.4451017e-2, -7.8215644e-2, -8.4166747e-2, -8.5166318e-2, -8.6980291e-2, -9.0938899e-2, -9.3153568e-2, -9.3672084e-2, -9.6012876e-2, -9.7460983e-2, -9.8154272e-2, -9.9153450e-2, -9.5783090e-2, -7.0552833e-2, -7.0552833e-2]),
                ),
        ],
        "qcd_y"    : np.array([-9.7774455e-3, -2.7803070e-2, -4.6322000e-2, -6.2975151e-2, -7.9511210e-2, -9.7716412e-2, -1.1201228e-1, -1.2465352e-1, -1.3992382e-1, -1.5403592e-1, -1.6342388e-1, -1.7449431e-1, -1.8365024e-1, -1.9465371e-1, -1.9962541e-1, -2.0517482e-1, -2.1229826e-1, -2.1722257e-1, -2.1993325e-1, -2.2407089e-1, -2.2699239e-1, -2.2865013e-1, -2.3047026e-1, -2.2613492e-1, -1.7065560e-1, -1.7065560e-1]),
        "qcd_min"  : np.array([-9.9495415e-3, -2.8039568e-2, -4.6641585e-2, -6.3405979e-2, -8.0058809e-2, -9.8477267e-2, -1.1292865e-1, -1.2565850e-1, -1.4105077e-1, -1.5520442e-1, -1.6480071e-1, -1.7585102e-1, -1.8515335e-1, -1.9618737e-1, -2.0119655e-1, -2.0679841e-1, -2.1382807e-1, -2.1888069e-1, -2.2163057e-1, -2.2580492e-1, -2.2871629e-1, -2.3041863e-1, -2.3219166e-1, -2.2801534e-1, -1.7200834e-1, -1.7200834e-1]),
        "qcd_max"  : np.array([-9.6561309e-3, -2.7599708e-2, -4.6069879e-2, -6.2637658e-2, -7.9083183e-2, -9.7073087e-2, -1.1122072e-1, -1.2380361e-1, -1.3896519e-1, -1.5306307e-1, -1.6222176e-1, -1.7335756e-1, -1.8235507e-1, -1.9334313e-1, -1.9829474e-1, -2.0380214e-1, -2.1106506e-1, -2.1583806e-1, -2.1851565e-1, -2.2261390e-1, -2.2555043e-1, -2.2716600e-1, -2.2905247e-1, -2.2452090e-1, -1.6964796e-1, -1.6964796e-1]),
        "y"        : np.array([-9.7774455e-3, -2.7803070e-2, -4.6322000e-2, -6.2975151e-2, -7.9511210e-2, -9.7716412e-2, -1.1201228e-1, -1.2465352e-1, -1.3992382e-1, -1.5403592e-1, -1.6342388e-1, -1.7449431e-1, -1.8365024e-1, -1.9465371e-1, -1.9962541e-1, -2.0517482e-1, -2.1229826e-1, -2.1722257e-1, -2.1993325e-1, -2.2407089e-1, -2.2699239e-1, -2.2865013e-1, -2.3047026e-1, -2.2613492e-1, -1.7065560e-1, -1.7065560e-1]),
        "ymin"     : np.array([-9.9495415e-3, -2.8039568e-2, -4.6641585e-2, -6.3405979e-2, -8.0058809e-2, -9.8477267e-2, -1.1292865e-1, -1.2565850e-1, -1.4105077e-1, -1.5520442e-1, -1.6480071e-1, -1.7585102e-1, -1.8515335e-1, -1.9618737e-1, -2.0119655e-1, -2.0679841e-1, -2.1382807e-1, -2.1888069e-1, -2.2163057e-1, -2.2580492e-1, -2.2871629e-1, -2.3041863e-1, -2.3219166e-1, -2.2801534e-1, -1.7200834e-1, -1.7200834e-1]),
        "ymax"     : np.array([-9.6561309e-3, -2.7599708e-2, -4.6069879e-2, -6.2637658e-2, -7.9083183e-2, -9.7073087e-2, -1.1122072e-1, -1.2380361e-1, -1.3896519e-1, -1.5306307e-1, -1.6222176e-1, -1.7335756e-1, -1.8235507e-1, -1.9334313e-1, -1.9829474e-1, -2.0380214e-1, -2.1106506e-1, -2.1583806e-1, -2.1851565e-1, -2.2261390e-1, -2.2555043e-1, -2.2716600e-1, -2.2905247e-1, -2.2452090e-1, -1.6964796e-1, -1.6964796e-1]),
        "channels" : [

        ],
    },
]
metadata = {
    "arxiv": r"",
    "description": r"",
    "hepdata": r"",
    "initial_state_1": r"2212",
    "initial_state_2": r"2212",
    "lumi_id_types": r"pdg_mc_ids",
    "mg5amc_repo": r"http://bazaar.launchpad.net/~maddevelopers/mg5amcnlo/3.3.1/",
    "mg5amc_revno": r"981",
    "patches": r"",
    "pineappl_gitversion": r"v0.5.4-49-g04650d9",
    "results_pdf": r"MSHT20nnlo_as118",
    "runcard_gitversion": r"7b5180d",
    "tau_min": r"",
    "user_cuts": r"",
    "x1_label": r"costh",
    "x1_label_tex": r"$\cos \theta^*$",
    "x1_unit": r"",
    "y_label": r"dsig/dcosth",
    "y_label_tex": r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\cos \theta^*}$",
    "y_unit": r"\pico\barn",
}
# end CLI data variables


if __name__ == "__main__":
    main(plot_panels)
"#;

const DRELL_YAN_MASS_SLICES_STR: &str = r#"#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pickle

# global variables coming from the CLI
title = r""
xlabel = r"$\cos \theta^*$"
ylabel = r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\cos \theta^*}$ [\si{pb}]"
xlog = False
ylog = False
scales = 7
plot_panels = {
    "int": False,
    "abs": True,
    "rel_ewonoff": True,
    "abs_pdfs": False,
    "ratio_pdf": False,
    "double_ratio_pdf": False,
    "rel_pdfunc": False,
    "rel_pdfpull": False,
}
output = r"NNPDF_DY_14TEV_BSM_AFB"

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
data = [
    {
        "slice_label"    : r"$\SI{200}{GeV} < M_{\ell\bar{\ell}} < \SI{500}{GeV}$",
        "x"        : np.array([-1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1]),
        "mid"      : np.array([-0.98, -0.94, -0.9, -0.86, -0.8200000000000001, -0.78, -0.74, -0.7, -0.66, -0.62, -0.5800000000000001, -0.54, -0.5, -0.45999999999999996, -0.42000000000000004, -0.38, -0.33999999999999997, -0.30000000000000004, -0.26, -0.22, -0.18, -0.14, -0.1, -0.06, -0.02, 0.02, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30000000000000004, 0.33999999999999997, 0.38, 0.42000000000000004, 0.45999999999999996, 0.5, 0.54, 0.5800000000000001, 0.62, 0.66, 0.7, 0.74, 0.78, 0.8200000000000001, 0.86, 0.9, 0.94, 0.98]),
        "pdf_results" : [
            (
                    r"NNPDF40\_nnlo\_as\_01180",
                    np.array([1.2030676e0, 1.1565560e0, 1.1130457e0, 1.0725099e0, 1.0387691e0, 1.0063348e0, 9.7581283e-1, 9.5075401e-1, 9.2789806e-1, 9.0985382e-1, 8.9310798e-1, 8.8180179e-1, 8.7261653e-1, 8.6606981e-1, 8.6360260e-1, 8.6373393e-1, 8.6791957e-1, 8.7409325e-1, 8.8588921e-1, 8.9893349e-1, 9.1509988e-1, 9.3716398e-1, 9.6097839e-1, 9.8834317e-1, 1.0180896e0, 1.0535762e0, 1.0895800e0, 1.1288776e0, 1.1737558e0, 1.2197631e0, 1.2689601e0, 1.3247334e0, 1.3811064e0, 1.4397822e0, 1.5032413e0, 1.5720274e0, 1.6419769e0, 1.7134034e0, 1.7913589e0, 1.8745908e0, 1.9576793e0, 2.0405341e0, 2.1340674e0, 2.2269273e0, 2.3239393e0, 2.4234584e0, 2.5266015e0, 2.6334832e0, 2.7449504e0, 2.8548557e0, 2.8548557e0]),
                    np.array([1.2030676e0, 1.1565560e0, 1.1130457e0, 1.0725099e0, 1.0387691e0, 1.0063348e0, 9.7581283e-1, 9.5075401e-1, 9.2789806e-1, 9.0985382e-1, 8.9310798e-1, 8.8180179e-1, 8.7261653e-1, 8.6606981e-1, 8.6360260e-1, 8.6373393e-1, 8.6791957e-1, 8.7409325e-1, 8.8588921e-1, 8.9893349e-1, 9.1509988e-1, 9.3716398e-1, 9.6097839e-1, 9.8834317e-1, 1.0180896e0, 1.0535762e0, 1.0895800e0, 1.1288776e0, 1.1737558e0, 1.2197631e0, 1.2689601e0, 1.3247334e0, 1.3811064e0, 1.4397822e0, 1.5032413e0, 1.5720274e0, 1.6419769e0, 1.7134034e0, 1.7913589e0, 1.8745908e0, 1.9576793e0, 2.0405341e0, 2.1340674e0, 2.2269273e0, 2.3239393e0, 2.4234584e0, 2.5266015e0, 2.6334832e0, 2.7449504e0, 2.8548557e0, 2.8548557e0]),
                    np.array([1.2030676e0, 1.1565560e0, 1.1130457e0, 1.0725099e0, 1.0387691e0, 1.0063348e0, 9.7581283e-1, 9.5075401e-1, 9.2789806e-1, 9.0985382e-1, 8.9310798e-1, 8.8180179e-1, 8.7261653e-1, 8.6606981e-1, 8.6360260e-1, 8.6373393e-1, 8.6791957e-1, 8.7409325e-1, 8.8588921e-1, 8.9893349e-1, 9.1509988e-1, 9.3716398e-1, 9.6097839e-1, 9.8834317e-1, 1.0180896e0, 1.0535762e0, 1.0895800e0, 1.1288776e0, 1.1737558e0, 1.2197631e0, 1.2689601e0, 1.3247334e0, 1.3811064e0, 1.4397822e0, 1.5032413e0, 1.5720274e0, 1.6419769e0, 1.7134034e0, 1.7913589e0, 1.8745908e0, 1.9576793e0, 2.0405341e0, 2.1340674e0, 2.2269273e0, 2.3239393e0, 2.4234584e0, 2.5266015e0, 2.6334832e0, 2.7449504e0, 2.8548557e0, 2.8548557e0]),
                ),
        ],
        "qcd_y"    : np.array([1.2030676e0, 1.1565560e0, 1.1130457e0, 1.0725099e0, 1.0387691e0, 1.0063348e0, 9.7581283e-1, 9.5075401e-1, 9.2789806e-1, 9.0985382e-1, 8.9310798e-1, 8.8180179e-1, 8.7261653e-1, 8.6606981e-1, 8.6360260e-1, 8.6373393e-1, 8.6791957e-1, 8.7409325e-1, 8.8588921e-1, 8.9893349e-1, 9.1509988e-1, 9.3716398e-1, 9.6097839e-1, 9.8834317e-1, 1.0180896e0, 1.0535762e0, 1.0895800e0, 1.1288776e0, 1.1737558e0, 1.2197631e0, 1.2689601e0, 1.3247334e0, 1.3811064e0, 1.4397822e0, 1.5032413e0, 1.5720274e0, 1.6419769e0, 1.7134034e0, 1.7913589e0, 1.8745908e0, 1.9576793e0, 2.0405341e0, 2.1340674e0, 2.2269273e0, 2.3239393e0, 2.4234584e0, 2.5266015e0, 2.6334832e0, 2.7449504e0, 2.8548557e0, 2.8548557e0]),
        "qcd_min"  : np.array([1.1476341e0, 1.1032579e0, 1.0617385e0, 1.0231508e0, 9.9097579e-1, 9.6014207e-1, 9.3113469e-1, 9.0728416e-1, 8.8560674e-1, 8.6851930e-1, 8.5267109e-1, 8.4203435e-1, 8.3342689e-1, 8.2736426e-1, 8.2520004e-1, 8.2555967e-1, 8.2975676e-1, 8.3587859e-1, 8.4733087e-1, 8.6002996e-1, 8.7572987e-1, 8.9705641e-1, 9.2007984e-1, 9.4649062e-1, 9.7517487e-1, 1.0093523e0, 1.0440650e0, 1.0819022e0, 1.1251039e0, 1.1693745e0, 1.2166754e0, 1.2702912e0, 1.3245299e0, 1.3809292e0, 1.4419071e0, 1.5079840e0, 1.5752031e0, 1.6438283e0, 1.7186871e0, 1.7987148e0, 1.8783985e0, 1.9580267e0, 2.0477772e0, 2.1369454e0, 2.2300739e0, 2.3256051e0, 2.4246078e0, 2.5271474e0, 2.6341671e0, 2.7396171e0, 2.7396171e0]),
        "qcd_max"  : np.array([1.2482967e0, 1.2000452e0, 1.1549126e0, 1.1127812e0, 1.0777638e0, 1.0440155e0, 1.0122494e0, 9.8620414e-1, 9.6237779e-1, 9.4354507e-1, 9.2605954e-1, 9.1419910e-1, 9.0453198e-1, 8.9757973e-1, 8.9485425e-1, 8.9478263e-1, 8.9894680e-1, 9.0514901e-1, 9.1721599e-1, 9.3052688e-1, 9.4705558e-1, 9.6970487e-1, 9.9414526e-1, 1.0222704e0, 1.0528665e0, 1.0894022e0, 1.1264355e0, 1.1669048e0, 1.2131264e0, 1.2605276e0, 1.3112501e0, 1.3687598e0, 1.4268454e0, 1.4873516e0, 1.5528115e0, 1.6237809e0, 1.6959271e0, 1.7696100e0, 1.8500634e0, 1.9358676e0, 2.0217144e0, 2.1071620e0, 2.2037519e0, 2.2995896e0, 2.3997365e0, 2.5024724e0, 2.6089569e0, 2.7193485e0, 2.8344013e0, 2.9479075e0, 2.9479075e0]),
        "y"        : np.array([1.2030676e0, 1.1565560e0, 1.1130457e0, 1.0725099e0, 1.0387691e0, 1.0063348e0, 9.7581283e-1, 9.5075401e-1, 9.2789806e-1, 9.0985382e-1, 8.9310798e-1, 8.8180179e-1, 8.7261653e-1, 8.6606981e-1, 8.6360260e-1, 8.6373393e-1, 8.6791957e-1, 8.7409325e-1, 8.8588921e-1, 8.9893349e-1, 9.1509988e-1, 9.3716398e-1, 9.6097839e-1, 9.8834317e-1, 1.0180896e0, 1.0535762e0, 1.0895800e0, 1.1288776e0, 1.1737558e0, 1.2197631e0, 1.2689601e0, 1.3247334e0, 1.3811064e0, 1.4397822e0, 1.5032413e0, 1.5720274e0, 1.6419769e0, 1.7134034e0, 1.7913589e0, 1.8745908e0, 1.9576793e0, 2.0405341e0, 2.1340674e0, 2.2269273e0, 2.3239393e0, 2.4234584e0, 2.5266015e0, 2.6334832e0, 2.7449504e0, 2.8548557e0, 2.8548557e0]),
        "ymin"     : np.array([1.1476341e0, 1.1032579e0, 1.0617385e0, 1.0231508e0, 9.9097579e-1, 9.6014207e-1, 9.3113469e-1, 9.0728416e-1, 8.8560674e-1, 8.6851930e-1, 8.5267109e-1, 8.4203435e-1, 8.3342689e-1, 8.2736426e-1, 8.2520004e-1, 8.2555967e-1, 8.2975676e-1, 8.3587859e-1, 8.4733087e-1, 8.6002996e-1, 8.7572987e-1, 8.9705641e-1, 9.2007984e-1, 9.4649062e-1, 9.7517487e-1, 1.0093523e0, 1.0440650e0, 1.0819022e0, 1.1251039e0, 1.1693745e0, 1.2166754e0, 1.2702912e0, 1.3245299e0, 1.3809292e0, 1.4419071e0, 1.5079840e0, 1.5752031e0, 1.6438283e0, 1.7186871e0, 1.7987148e0, 1.8783985e0, 1.9580267e0, 2.0477772e0, 2.1369454e0, 2.2300739e0, 2.3256051e0, 2.4246078e0, 2.5271474e0, 2.6341671e0, 2.7396171e0, 2.7396171e0]),
        "ymax"     : np.array([1.2482967e0, 1.2000452e0, 1.1549126e0, 1.1127812e0, 1.0777638e0, 1.0440155e0, 1.0122494e0, 9.8620414e-1, 9.6237779e-1, 9.4354507e-1, 9.2605954e-1, 9.1419910e-1, 9.0453198e-1, 8.9757973e-1, 8.9485425e-1, 8.9478263e-1, 8.9894680e-1, 9.0514901e-1, 9.1721599e-1, 9.3052688e-1, 9.4705558e-1, 9.6970487e-1, 9.9414526e-1, 1.0222704e0, 1.0528665e0, 1.0894022e0, 1.1264355e0, 1.1669048e0, 1.2131264e0, 1.2605276e0, 1.3112501e0, 1.3687598e0, 1.4268454e0, 1.4873516e0, 1.5528115e0, 1.6237809e0, 1.6959271e0, 1.7696100e0, 1.8500634e0, 1.9358676e0, 2.0217144e0, 2.1071620e0, 2.2037519e0, 2.2995896e0, 2.3997365e0, 2.5024724e0, 2.6089569e0, 2.7193485e0, 2.8344013e0, 2.9479075e0, 2.9479075e0]),
        "channels" : [
            (r"$\mathrm{u}\bar{\mathrm{u}} + \mathrm{c}\bar{\mathrm{c}}$", np.array([6.8792943e-1, 6.6125562e-1, 6.3777395e-1, 6.1486127e-1, 5.9741951e-1, 5.7963946e-1, 5.6412672e-1, 5.5280097e-1, 5.4155222e-1, 5.3440638e-1, 5.2772453e-1, 5.2540250e-1, 5.2312552e-1, 5.2336775e-1, 5.2636757e-1, 5.3147731e-1, 5.3918965e-1, 5.4717287e-1, 5.6056341e-1, 5.7337587e-1, 5.8836121e-1, 6.0791258e-1, 6.2862964e-1, 6.4938700e-1, 6.7588930e-1, 7.0371242e-1, 7.3209233e-1, 7.6402657e-1, 7.9766655e-1, 8.3289425e-1, 8.7035507e-1, 9.1152732e-1, 9.5376500e-1, 9.9877556e-1, 1.0445475e0, 1.0949400e0, 1.1456393e0, 1.1986085e0, 1.2540414e0, 1.3152467e0, 1.3747164e0, 1.4336967e0, 1.5021011e0, 1.5683542e0, 1.6368719e0, 1.7086013e0, 1.7810052e0, 1.8565613e0, 1.9354914e0, 2.0139646e0, 2.0139646e0])),
            (r"$\mathrm{d}\bar{\mathrm{d}} + \mathrm{s}\bar{\mathrm{s}} + \mathrm{b}\bar{\mathrm{b}}$", np.array([5.1513820e-1, 4.9530042e-1, 4.7527174e-1, 4.5764866e-1, 4.4134961e-1, 4.2669534e-1, 4.1168610e-1, 3.9795304e-1, 3.8634584e-1, 3.7544744e-1, 3.6538344e-1, 3.5639929e-1, 3.4949102e-1, 3.4270206e-1, 3.3723503e-1, 3.3225662e-1, 3.2872992e-1, 3.2692037e-1, 3.2532580e-1, 3.2555763e-1, 3.2673866e-1, 3.2925139e-1, 3.3234876e-1, 3.3895617e-1, 3.4220033e-1, 3.4986376e-1, 3.5748772e-1, 3.6485107e-1, 3.7608929e-1, 3.8686886e-1, 3.9860504e-1, 4.1320603e-1, 4.2734137e-1, 4.4100662e-1, 4.5869384e-1, 4.7708743e-1, 4.9633757e-1, 5.1479498e-1, 5.3731751e-1, 5.5934418e-1, 5.8296294e-1, 6.0683737e-1, 6.3196624e-1, 6.5857309e-1, 6.8706738e-1, 7.1485710e-1, 7.4559634e-1, 7.7692189e-1, 8.0945902e-1, 8.4089109e-1, 8.4089109e-1]))
        ],
    },
    {
        "slice_label"    : r"$\SI{500}{GeV} < M_{\ell\bar{\ell}} < \SI{1000}{GeV}$",
        "x"        : np.array([-1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1]),
        "mid"      : np.array([-0.98, -0.94, -0.9, -0.86, -0.8200000000000001, -0.78, -0.74, -0.7, -0.66, -0.62, -0.5800000000000001, -0.54, -0.5, -0.45999999999999996, -0.42000000000000004, -0.38, -0.33999999999999997, -0.30000000000000004, -0.26, -0.22, -0.18, -0.14, -0.1, -0.06, -0.02, 0.02, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30000000000000004, 0.33999999999999997, 0.38, 0.42000000000000004, 0.45999999999999996, 0.5, 0.54, 0.5800000000000001, 0.62, 0.66, 0.7, 0.74, 0.78, 0.8200000000000001, 0.86, 0.9, 0.94, 0.98]),
        "pdf_results" : [
            (
                    r"NNPDF40\_nnlo\_as\_01180",
                    np.array([4.0594871e-2, 3.8970516e-2, 3.7589145e-2, 3.6287278e-2, 3.5138659e-2, 3.4083575e-2, 3.3138203e-2, 3.2369514e-2, 3.1687339e-2, 3.1130680e-2, 3.0694113e-2, 3.0348044e-2, 3.0193553e-2, 3.0118871e-2, 3.0174121e-2, 3.0339279e-2, 3.0573898e-2, 3.1016543e-2, 3.1542002e-2, 3.2184767e-2, 3.2964449e-2, 3.3870335e-2, 3.4866769e-2, 3.6057960e-2, 3.7346103e-2, 3.8648244e-2, 4.0227527e-2, 4.1803650e-2, 4.3536855e-2, 4.5451532e-2, 4.7387640e-2, 4.9457643e-2, 5.1733512e-2, 5.4105127e-2, 5.6588818e-2, 5.9238183e-2, 6.1866112e-2, 6.4704350e-2, 6.7748017e-2, 7.0789065e-2, 7.4060501e-2, 7.7264134e-2, 8.0877170e-2, 8.4406261e-2, 8.8038767e-2, 9.1920534e-2, 9.5741259e-2, 9.9856666e-2, 1.0403280e-1, 1.0829409e-1, 1.0829409e-1]),
                    np.array([4.0594871e-2, 3.8970516e-2, 3.7589145e-2, 3.6287278e-2, 3.5138659e-2, 3.4083575e-2, 3.3138203e-2, 3.2369514e-2, 3.1687339e-2, 3.1130680e-2, 3.0694113e-2, 3.0348044e-2, 3.0193553e-2, 3.0118871e-2, 3.0174121e-2, 3.0339279e-2, 3.0573898e-2, 3.1016543e-2, 3.1542002e-2, 3.2184767e-2, 3.2964449e-2, 3.3870335e-2, 3.4866769e-2, 3.6057960e-2, 3.7346103e-2, 3.8648244e-2, 4.0227527e-2, 4.1803650e-2, 4.3536855e-2, 4.5451532e-2, 4.7387640e-2, 4.9457643e-2, 5.1733512e-2, 5.4105127e-2, 5.6588818e-2, 5.9238183e-2, 6.1866112e-2, 6.4704350e-2, 6.7748017e-2, 7.0789065e-2, 7.4060501e-2, 7.7264134e-2, 8.0877170e-2, 8.4406261e-2, 8.8038767e-2, 9.1920534e-2, 9.5741259e-2, 9.9856666e-2, 1.0403280e-1, 1.0829409e-1, 1.0829409e-1]),
                    np.array([4.0594871e-2, 3.8970516e-2, 3.7589145e-2, 3.6287278e-2, 3.5138659e-2, 3.4083575e-2, 3.3138203e-2, 3.2369514e-2, 3.1687339e-2, 3.1130680e-2, 3.0694113e-2, 3.0348044e-2, 3.0193553e-2, 3.0118871e-2, 3.0174121e-2, 3.0339279e-2, 3.0573898e-2, 3.1016543e-2, 3.1542002e-2, 3.2184767e-2, 3.2964449e-2, 3.3870335e-2, 3.4866769e-2, 3.6057960e-2, 3.7346103e-2, 3.8648244e-2, 4.0227527e-2, 4.1803650e-2, 4.3536855e-2, 4.5451532e-2, 4.7387640e-2, 4.9457643e-2, 5.1733512e-2, 5.4105127e-2, 5.6588818e-2, 5.9238183e-2, 6.1866112e-2, 6.4704350e-2, 6.7748017e-2, 7.0789065e-2, 7.4060501e-2, 7.7264134e-2, 8.0877170e-2, 8.4406261e-2, 8.8038767e-2, 9.1920534e-2, 9.5741259e-2, 9.9856666e-2, 1.0403280e-1, 1.0829409e-1, 1.0829409e-1]),
                ),
        ],
        "qcd_y"    : np.array([4.0594871e-2, 3.8970516e-2, 3.7589145e-2, 3.6287278e-2, 3.5138659e-2, 3.4083575e-2, 3.3138203e-2, 3.2369514e-2, 3.1687339e-2, 3.1130680e-2, 3.0694113e-2, 3.0348044e-2, 3.0193553e-2, 3.0118871e-2, 3.0174121e-2, 3.0339279e-2, 3.0573898e-2, 3.1016543e-2, 3.1542002e-2, 3.2184767e-2, 3.2964449e-2, 3.3870335e-2, 3.4866769e-2, 3.6057960e-2, 3.7346103e-2, 3.8648244e-2, 4.0227527e-2, 4.1803650e-2, 4.3536855e-2, 4.5451532e-2, 4.7387640e-2, 4.9457643e-2, 5.1733512e-2, 5.4105127e-2, 5.6588818e-2, 5.9238183e-2, 6.1866112e-2, 6.4704350e-2, 6.7748017e-2, 7.0789065e-2, 7.4060501e-2, 7.7264134e-2, 8.0877170e-2, 8.4406261e-2, 8.8038767e-2, 9.1920534e-2, 9.5741259e-2, 9.9856666e-2, 1.0403280e-1, 1.0829409e-1, 1.0829409e-1]),
        "qcd_min"  : np.array([4.0037777e-2, 3.8435726e-2, 3.7073393e-2, 3.5789220e-2, 3.4655883e-2, 3.3613703e-2, 3.2680468e-2, 3.1921166e-2, 3.1247155e-2, 3.0695910e-2, 3.0263715e-2, 2.9921326e-2, 2.9765787e-2, 2.9690226e-2, 2.9743185e-2, 2.9903344e-2, 3.0131781e-2, 3.0565874e-2, 3.1081207e-2, 3.1711632e-2, 3.2477115e-2, 3.3366364e-2, 3.4346364e-2, 3.5517091e-2, 3.6783757e-2, 3.8063466e-2, 3.9616935e-2, 4.1166189e-2, 4.2871525e-2, 4.4755097e-2, 4.6659767e-2, 4.8697082e-2, 5.0934848e-2, 5.3268112e-2, 5.5712599e-2, 5.8320365e-2, 6.0905504e-2, 6.3698651e-2, 6.6694216e-2, 6.9687107e-2, 7.2907813e-2, 7.6059137e-2, 7.9616205e-2, 8.3090171e-2, 8.6665976e-2, 9.0484669e-2, 9.4247104e-2, 9.8297870e-2, 1.0240852e-1, 1.0660338e-1, 1.0660338e-1]),
        "qcd_max"  : np.array([4.1061145e-2, 3.9418149e-2, 3.8020807e-2, 3.6704197e-2, 3.5542954e-2, 3.4477622e-2, 3.3522380e-2, 3.2746246e-2, 3.2057661e-2, 3.1497253e-2, 3.1057591e-2, 3.0708825e-2, 3.0556281e-2, 3.0483016e-2, 3.0540730e-2, 3.0711032e-2, 3.0951879e-2, 3.1402556e-2, 3.1937511e-2, 3.2591833e-2, 3.3384624e-2, 3.4305895e-2, 3.5317075e-2, 3.6526854e-2, 3.7834301e-2, 3.9156804e-2, 4.0759176e-2, 4.2359581e-2, 4.4117603e-2, 4.6059981e-2, 4.8024100e-2, 5.0123011e-2, 5.2433124e-2, 5.4838902e-2, 5.7357234e-2, 6.0043264e-2, 6.2709323e-2, 6.5587507e-2, 6.8673601e-2, 7.1757263e-2, 7.5073207e-2, 7.8323545e-2, 8.1985664e-2, 8.5563285e-2, 8.9245686e-2, 9.3183548e-2, 9.7055206e-2, 1.0122759e-1, 1.0546145e-1, 1.0978108e-1, 1.0978108e-1]),
        "y"        : np.array([4.0594871e-2, 3.8970516e-2, 3.7589145e-2, 3.6287278e-2, 3.5138659e-2, 3.4083575e-2, 3.3138203e-2, 3.2369514e-2, 3.1687339e-2, 3.1130680e-2, 3.0694113e-2, 3.0348044e-2, 3.0193553e-2, 3.0118871e-2, 3.0174121e-2, 3.0339279e-2, 3.0573898e-2, 3.1016543e-2, 3.1542002e-2, 3.2184767e-2, 3.2964449e-2, 3.3870335e-2, 3.4866769e-2, 3.6057960e-2, 3.7346103e-2, 3.8648244e-2, 4.0227527e-2, 4.1803650e-2, 4.3536855e-2, 4.5451532e-2, 4.7387640e-2, 4.9457643e-2, 5.1733512e-2, 5.4105127e-2, 5.6588818e-2, 5.9238183e-2, 6.1866112e-2, 6.4704350e-2, 6.7748017e-2, 7.0789065e-2, 7.4060501e-2, 7.7264134e-2, 8.0877170e-2, 8.4406261e-2, 8.8038767e-2, 9.1920534e-2, 9.5741259e-2, 9.9856666e-2, 1.0403280e-1, 1.0829409e-1, 1.0829409e-1]),
        "ymin"     : np.array([4.0037777e-2, 3.8435726e-2, 3.7073393e-2, 3.5789220e-2, 3.4655883e-2, 3.3613703e-2, 3.2680468e-2, 3.1921166e-2, 3.1247155e-2, 3.0695910e-2, 3.0263715e-2, 2.9921326e-2, 2.9765787e-2, 2.9690226e-2, 2.9743185e-2, 2.9903344e-2, 3.0131781e-2, 3.0565874e-2, 3.1081207e-2, 3.1711632e-2, 3.2477115e-2, 3.3366364e-2, 3.4346364e-2, 3.5517091e-2, 3.6783757e-2, 3.8063466e-2, 3.9616935e-2, 4.1166189e-2, 4.2871525e-2, 4.4755097e-2, 4.6659767e-2, 4.8697082e-2, 5.0934848e-2, 5.3268112e-2, 5.5712599e-2, 5.8320365e-2, 6.0905504e-2, 6.3698651e-2, 6.6694216e-2, 6.9687107e-2, 7.2907813e-2, 7.6059137e-2, 7.9616205e-2, 8.3090171e-2, 8.6665976e-2, 9.0484669e-2, 9.4247104e-2, 9.8297870e-2, 1.0240852e-1, 1.0660338e-1, 1.0660338e-1]),
        "ymax"     : np.array([4.1061145e-2, 3.9418149e-2, 3.8020807e-2, 3.6704197e-2, 3.5542954e-2, 3.4477622e-2, 3.3522380e-2, 3.2746246e-2, 3.2057661e-2, 3.1497253e-2, 3.1057591e-2, 3.0708825e-2, 3.0556281e-2, 3.0483016e-2, 3.0540730e-2, 3.0711032e-2, 3.0951879e-2, 3.1402556e-2, 3.1937511e-2, 3.2591833e-2, 3.3384624e-2, 3.4305895e-2, 3.5317075e-2, 3.6526854e-2, 3.7834301e-2, 3.9156804e-2, 4.0759176e-2, 4.2359581e-2, 4.4117603e-2, 4.6059981e-2, 4.8024100e-2, 5.0123011e-2, 5.2433124e-2, 5.4838902e-2, 5.7357234e-2, 6.0043264e-2, 6.2709323e-2, 6.5587507e-2, 6.8673601e-2, 7.1757263e-2, 7.5073207e-2, 7.8323545e-2, 8.1985664e-2, 8.5563285e-2, 8.9245686e-2, 9.3183548e-2, 9.7055206e-2, 1.0122759e-1, 1.0546145e-1, 1.0978108e-1, 1.0978108e-1]),
        "channels" : [
            (r"$\mathrm{u}\bar{\mathrm{u}} + \mathrm{c}\bar{\mathrm{c}}$", np.array([2.6199757e-2, 2.5141774e-2, 2.4301385e-2, 2.3475409e-2, 2.2742818e-2, 2.2111053e-2, 2.1562477e-2, 2.1150586e-2, 2.0753029e-2, 2.0527176e-2, 2.0286675e-2, 2.0198805e-2, 2.0196503e-2, 2.0277685e-2, 2.0455607e-2, 2.0654565e-2, 2.0963610e-2, 2.1393063e-2, 2.1909943e-2, 2.2517162e-2, 2.3181602e-2, 2.3967128e-2, 2.4795910e-2, 2.5775243e-2, 2.6848306e-2, 2.7911025e-2, 2.9189817e-2, 3.0427277e-2, 3.1802230e-2, 3.3315181e-2, 3.4815501e-2, 3.6458797e-2, 3.8188509e-2, 4.0036615e-2, 4.1945442e-2, 4.3957344e-2, 4.6030617e-2, 4.8180061e-2, 5.0492953e-2, 5.2808985e-2, 5.5264729e-2, 5.7718158e-2, 6.0458116e-2, 6.3119505e-2, 6.5844881e-2, 6.8772294e-2, 7.1647003e-2, 7.4759927e-2, 7.7907487e-2, 8.1060130e-2, 8.1060130e-2])),
            (r"$\mathrm{d}\bar{\mathrm{d}} + \mathrm{s}\bar{\mathrm{s}} + \mathrm{b}\bar{\mathrm{b}}$", np.array([1.4395113e-2, 1.3828743e-2, 1.3287760e-2, 1.2811869e-2, 1.2395841e-2, 1.1972522e-2, 1.1575727e-2, 1.1218927e-2, 1.0934310e-2, 1.0603503e-2, 1.0407438e-2, 1.0149239e-2, 9.9970507e-3, 9.8411857e-3, 9.7185140e-3, 9.6847141e-3, 9.6102880e-3, 9.6234798e-3, 9.6320596e-3, 9.6676052e-3, 9.7828473e-3, 9.9032066e-3, 1.0070859e-2, 1.0282716e-2, 1.0497797e-2, 1.0737219e-2, 1.1037711e-2, 1.1376372e-2, 1.1734625e-2, 1.2136350e-2, 1.2572139e-2, 1.2998846e-2, 1.3545003e-2, 1.4068512e-2, 1.4643376e-2, 1.5280839e-2, 1.5835494e-2, 1.6524288e-2, 1.7255064e-2, 1.7980080e-2, 1.8795773e-2, 1.9545976e-2, 2.0419054e-2, 2.1286756e-2, 2.2193887e-2, 2.3148239e-2, 2.4094256e-2, 2.5096739e-2, 2.6125312e-2, 2.7233960e-2, 2.7233960e-2]))
        ],
    },
    {
        "slice_label"    : r"$\SI{1000}{GeV} < M_{\ell\bar{\ell}} < \SI{3000}{GeV}$",
        "x"        : np.array([-1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1]),
        "mid"      : np.array([-0.98, -0.94, -0.9, -0.86, -0.8200000000000001, -0.78, -0.74, -0.7, -0.66, -0.62, -0.5800000000000001, -0.54, -0.5, -0.45999999999999996, -0.42000000000000004, -0.38, -0.33999999999999997, -0.30000000000000004, -0.26, -0.22, -0.18, -0.14, -0.1, -0.06, -0.02, 0.02, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30000000000000004, 0.33999999999999997, 0.38, 0.42000000000000004, 0.45999999999999996, 0.5, 0.54, 0.5800000000000001, 0.62, 0.66, 0.7, 0.74, 0.78, 0.8200000000000001, 0.86, 0.9, 0.94, 0.98]),
        "pdf_results" : [
            (
                    r"NNPDF40\_nnlo\_as\_01180",
                    np.array([2.5617082e-3, 2.4598979e-3, 2.3739593e-3, 2.2904004e-3, 2.2169883e-3, 2.1545803e-3, 2.0954206e-3, 2.0498675e-3, 2.0070737e-3, 1.9754803e-3, 1.9463467e-3, 1.9278404e-3, 1.9236438e-3, 1.9198349e-3, 1.9306588e-3, 1.9387108e-3, 1.9583331e-3, 1.9907351e-3, 2.0286806e-3, 2.0733607e-3, 2.1270186e-3, 2.1863925e-3, 2.2577539e-3, 2.3390006e-3, 2.4212326e-3, 2.5112892e-3, 2.6161915e-3, 2.7238882e-3, 2.8385163e-3, 2.9647013e-3, 3.0959433e-3, 3.2331088e-3, 3.3817486e-3, 3.5385707e-3, 3.7031757e-3, 3.8780884e-3, 4.0550706e-3, 4.2388498e-3, 4.4433747e-3, 4.6413012e-3, 4.8539670e-3, 5.0680268e-3, 5.3044183e-3, 5.5388017e-3, 5.7755439e-3, 6.0326566e-3, 6.2848452e-3, 6.5558742e-3, 6.8328348e-3, 7.1072149e-3, 7.1072149e-3]),
                    np.array([2.5617082e-3, 2.4598979e-3, 2.3739593e-3, 2.2904004e-3, 2.2169883e-3, 2.1545803e-3, 2.0954206e-3, 2.0498675e-3, 2.0070737e-3, 1.9754803e-3, 1.9463467e-3, 1.9278404e-3, 1.9236438e-3, 1.9198349e-3, 1.9306588e-3, 1.9387108e-3, 1.9583331e-3, 1.9907351e-3, 2.0286806e-3, 2.0733607e-3, 2.1270186e-3, 2.1863925e-3, 2.2577539e-3, 2.3390006e-3, 2.4212326e-3, 2.5112892e-3, 2.6161915e-3, 2.7238882e-3, 2.8385163e-3, 2.9647013e-3, 3.0959433e-3, 3.2331088e-3, 3.3817486e-3, 3.5385707e-3, 3.7031757e-3, 3.8780884e-3, 4.0550706e-3, 4.2388498e-3, 4.4433747e-3, 4.6413012e-3, 4.8539670e-3, 5.0680268e-3, 5.3044183e-3, 5.5388017e-3, 5.7755439e-3, 6.0326566e-3, 6.2848452e-3, 6.5558742e-3, 6.8328348e-3, 7.1072149e-3, 7.1072149e-3]),
                    np.array([2.5617082e-3, 2.4598979e-3, 2.3739593e-3, 2.2904004e-3, 2.2169883e-3, 2.1545803e-3, 2.0954206e-3, 2.0498675e-3, 2.0070737e-3, 1.9754803e-3, 1.9463467e-3, 1.9278404e-3, 1.9236438e-3, 1.9198349e-3, 1.9306588e-3, 1.9387108e-3, 1.9583331e-3, 1.9907351e-3, 2.0286806e-3, 2.0733607e-3, 2.1270186e-3, 2.1863925e-3, 2.2577539e-3, 2.3390006e-3, 2.4212326e-3, 2.5112892e-3, 2.6161915e-3, 2.7238882e-3, 2.8385163e-3, 2.9647013e-3, 3.0959433e-3, 3.2331088e-3, 3.3817486e-3, 3.5385707e-3, 3.7031757e-3, 3.8780884e-3, 4.0550706e-3, 4.2388498e-3, 4.4433747e-3, 4.6413012e-3, 4.8539670e-3, 5.0680268e-3, 5.3044183e-3, 5.5388017e-3, 5.7755439e-3, 6.0326566e-3, 6.2848452e-3, 6.5558742e-3, 6.8328348e-3, 7.1072149e-3, 7.1072149e-3]),
                ),
        ],
        "qcd_y"    : np.array([2.5617082e-3, 2.4598979e-3, 2.3739593e-3, 2.2904004e-3, 2.2169883e-3, 2.1545803e-3, 2.0954206e-3, 2.0498675e-3, 2.0070737e-3, 1.9754803e-3, 1.9463467e-3, 1.9278404e-3, 1.9236438e-3, 1.9198349e-3, 1.9306588e-3, 1.9387108e-3, 1.9583331e-3, 1.9907351e-3, 2.0286806e-3, 2.0733607e-3, 2.1270186e-3, 2.1863925e-3, 2.2577539e-3, 2.3390006e-3, 2.4212326e-3, 2.5112892e-3, 2.6161915e-3, 2.7238882e-3, 2.8385163e-3, 2.9647013e-3, 3.0959433e-3, 3.2331088e-3, 3.3817486e-3, 3.5385707e-3, 3.7031757e-3, 3.8780884e-3, 4.0550706e-3, 4.2388498e-3, 4.4433747e-3, 4.6413012e-3, 4.8539670e-3, 5.0680268e-3, 5.3044183e-3, 5.5388017e-3, 5.7755439e-3, 6.0326566e-3, 6.2848452e-3, 6.5558742e-3, 6.8328348e-3, 7.1072149e-3, 7.1072149e-3]),
        "qcd_min"  : np.array([2.4433990e-3, 2.3462555e-3, 2.2643201e-3, 2.1845478e-3, 2.1145236e-3, 2.0548881e-3, 1.9983845e-3, 1.9548476e-3, 1.9139322e-3, 1.8836686e-3, 1.8557506e-3, 1.8379680e-3, 1.8337991e-3, 1.8299901e-3, 1.8401805e-3, 1.8476284e-3, 1.8662068e-3, 1.8968667e-3, 1.9328299e-3, 1.9751895e-3, 2.0260910e-3, 2.0824733e-3, 2.1502855e-3, 2.2275197e-3, 2.3055794e-3, 2.3911361e-3, 2.4909493e-3, 2.5933063e-3, 2.7022640e-3, 2.8222590e-3, 2.9470724e-3, 3.0775451e-3, 3.2189071e-3, 3.3680276e-3, 3.5246274e-3, 3.6910487e-3, 3.8593518e-3, 4.0342136e-3, 4.2287372e-3, 4.4170780e-3, 4.6194516e-3, 4.8230030e-3, 5.0480045e-3, 5.2710327e-3, 5.4962622e-3, 5.7409237e-3, 5.9809253e-3, 6.2388890e-3, 6.5023785e-3, 6.7635515e-3, 6.7635515e-3]),
        "qcd_max"  : np.array([2.6900604e-3, 2.5831920e-3, 2.4929066e-3, 2.4052452e-3, 2.3281584e-3, 2.2627536e-3, 2.2007167e-3, 2.1529838e-3, 2.1081610e-3, 2.0751367e-3, 2.0446966e-3, 2.0254163e-3, 2.0212045e-3, 2.0174111e-3, 2.0289322e-3, 2.0376618e-3, 2.0584262e-3, 2.0927416e-3, 2.1328564e-3, 2.1800780e-3, 2.2367509e-3, 2.2993910e-3, 2.3746249e-3, 2.4602457e-3, 2.5470394e-3, 2.6420082e-3, 2.7524493e-3, 2.8659720e-3, 2.9867879e-3, 3.1197183e-3, 3.2579655e-3, 3.4024228e-3, 3.5589936e-3, 3.7242136e-3, 3.8975378e-3, 4.0816963e-3, 4.2681397e-3, 4.4616317e-3, 4.6770565e-3, 4.8854207e-3, 5.1092919e-3, 5.3348103e-3, 5.5835965e-3, 5.8303478e-3, 6.0796332e-3, 6.3502993e-3, 6.6157573e-3, 6.9010042e-3, 7.1926442e-3, 7.4813962e-3, 7.4813962e-3]),
        "y"        : np.array([2.5617082e-3, 2.4598979e-3, 2.3739593e-3, 2.2904004e-3, 2.2169883e-3, 2.1545803e-3, 2.0954206e-3, 2.0498675e-3, 2.0070737e-3, 1.9754803e-3, 1.9463467e-3, 1.9278404e-3, 1.9236438e-3, 1.9198349e-3, 1.9306588e-3, 1.9387108e-3, 1.9583331e-3, 1.9907351e-3, 2.0286806e-3, 2.0733607e-3, 2.1270186e-3, 2.1863925e-3, 2.2577539e-3, 2.3390006e-3, 2.4212326e-3, 2.5112892e-3, 2.6161915e-3, 2.7238882e-3, 2.8385163e-3, 2.9647013e-3, 3.0959433e-3, 3.2331088e-3, 3.3817486e-3, 3.5385707e-3, 3.7031757e-3, 3.8780884e-3, 4.0550706e-3, 4.2388498e-3, 4.4433747e-3, 4.6413012e-3, 4.8539670e-3, 5.0680268e-3, 5.3044183e-3, 5.5388017e-3, 5.7755439e-3, 6.0326566e-3, 6.2848452e-3, 6.5558742e-3, 6.8328348e-3, 7.1072149e-3, 7.1072149e-3]),
        "ymin"     : np.array([2.4433990e-3, 2.3462555e-3, 2.2643201e-3, 2.1845478e-3, 2.1145236e-3, 2.0548881e-3, 1.9983845e-3, 1.9548476e-3, 1.9139322e-3, 1.8836686e-3, 1.8557506e-3, 1.8379680e-3, 1.8337991e-3, 1.8299901e-3, 1.8401805e-3, 1.8476284e-3, 1.8662068e-3, 1.8968667e-3, 1.9328299e-3, 1.9751895e-3, 2.0260910e-3, 2.0824733e-3, 2.1502855e-3, 2.2275197e-3, 2.3055794e-3, 2.3911361e-3, 2.4909493e-3, 2.5933063e-3, 2.7022640e-3, 2.8222590e-3, 2.9470724e-3, 3.0775451e-3, 3.2189071e-3, 3.3680276e-3, 3.5246274e-3, 3.6910487e-3, 3.8593518e-3, 4.0342136e-3, 4.2287372e-3, 4.4170780e-3, 4.6194516e-3, 4.8230030e-3, 5.0480045e-3, 5.2710327e-3, 5.4962622e-3, 5.7409237e-3, 5.9809253e-3, 6.2388890e-3, 6.5023785e-3, 6.7635515e-3, 6.7635515e-3]),
        "ymax"     : np.array([2.6900604e-3, 2.5831920e-3, 2.4929066e-3, 2.4052452e-3, 2.3281584e-3, 2.2627536e-3, 2.2007167e-3, 2.1529838e-3, 2.1081610e-3, 2.0751367e-3, 2.0446966e-3, 2.0254163e-3, 2.0212045e-3, 2.0174111e-3, 2.0289322e-3, 2.0376618e-3, 2.0584262e-3, 2.0927416e-3, 2.1328564e-3, 2.1800780e-3, 2.2367509e-3, 2.2993910e-3, 2.3746249e-3, 2.4602457e-3, 2.5470394e-3, 2.6420082e-3, 2.7524493e-3, 2.8659720e-3, 2.9867879e-3, 3.1197183e-3, 3.2579655e-3, 3.4024228e-3, 3.5589936e-3, 3.7242136e-3, 3.8975378e-3, 4.0816963e-3, 4.2681397e-3, 4.4616317e-3, 4.6770565e-3, 4.8854207e-3, 5.1092919e-3, 5.3348103e-3, 5.5835965e-3, 5.8303478e-3, 6.0796332e-3, 6.3502993e-3, 6.6157573e-3, 6.9010042e-3, 7.1926442e-3, 7.4813962e-3, 7.4813962e-3]),
        "channels" : [
            (r"$\mathrm{u}\bar{\mathrm{u}} + \mathrm{c}\bar{\mathrm{c}}$", np.array([1.7616291e-3, 1.6911649e-3, 1.6342583e-3, 1.5774106e-3, 1.5296227e-3, 1.4870592e-3, 1.4501547e-3, 1.4231677e-3, 1.3961418e-3, 1.3812502e-3, 1.3643810e-3, 1.3567734e-3, 1.3608923e-3, 1.3632562e-3, 1.3766557e-3, 1.3895131e-3, 1.4112073e-3, 1.4417942e-3, 1.4765595e-3, 1.5173010e-3, 1.5614209e-3, 1.6126977e-3, 1.6721940e-3, 1.7400618e-3, 1.8065100e-3, 1.8812026e-3, 1.9661448e-3, 2.0527345e-3, 2.1433356e-3, 2.2450972e-3, 2.3489923e-3, 2.4583112e-3, 2.5749841e-3, 2.6984465e-3, 2.8287488e-3, 2.9650118e-3, 3.1027853e-3, 3.2480645e-3, 3.4074955e-3, 3.5603053e-3, 3.7263245e-3, 3.8911284e-3, 4.0776944e-3, 4.2584114e-3, 4.4403143e-3, 4.6397117e-3, 4.8333550e-3, 5.0415868e-3, 5.2559809e-3, 5.4675131e-3, 5.4675131e-3])),
            (r"$\mathrm{d}\bar{\mathrm{d}} + \mathrm{s}\bar{\mathrm{s}} + \mathrm{b}\bar{\mathrm{b}}$", np.array([8.0007914e-4, 7.6873302e-4, 7.3970097e-4, 7.1298983e-4, 6.8736553e-4, 6.6752119e-4, 6.4526588e-4, 6.2669976e-4, 6.1093192e-4, 5.9423014e-4, 5.8196564e-4, 5.7106698e-4, 5.6275141e-4, 5.5657872e-4, 5.5400311e-4, 5.4919766e-4, 5.4712581e-4, 5.4894089e-4, 5.5212109e-4, 5.5605969e-4, 5.6559775e-4, 5.7369482e-4, 5.8555990e-4, 5.9893885e-4, 6.1472260e-4, 6.3008659e-4, 6.5004668e-4, 6.7115370e-4, 6.9518073e-4, 7.1960411e-4, 7.4695103e-4, 7.7479758e-4, 8.0676452e-4, 8.4012417e-4, 8.7442688e-4, 9.1307666e-4, 9.5228530e-4, 9.9078535e-4, 1.0358792e-3, 1.0809959e-3, 1.1276426e-3, 1.1768984e-3, 1.2267239e-3, 1.2803903e-3, 1.3352297e-3, 1.3929449e-3, 1.4514902e-3, 1.5142873e-3, 1.5768539e-3, 1.6397017e-3, 1.6397017e-3]))
        ],
    },
    {
        "slice_label"    : r"$\SI{3000}{GeV} < M_{\ell\bar{\ell}} < \SI{5000}{GeV}$",
        "x"        : np.array([-1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1]),
        "mid"      : np.array([-0.98, -0.94, -0.9, -0.86, -0.8200000000000001, -0.78, -0.74, -0.7, -0.66, -0.62, -0.5800000000000001, -0.54, -0.5, -0.45999999999999996, -0.42000000000000004, -0.38, -0.33999999999999997, -0.30000000000000004, -0.26, -0.22, -0.18, -0.14, -0.1, -0.06, -0.02, 0.02, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30000000000000004, 0.33999999999999997, 0.38, 0.42000000000000004, 0.45999999999999996, 0.5, 0.54, 0.5800000000000001, 0.62, 0.66, 0.7, 0.74, 0.78, 0.8200000000000001, 0.86, 0.9, 0.94, 0.98]),
        "pdf_results" : [
            (
                    r"NNPDF40\_nnlo\_as\_01180",
                    np.array([6.2465996e-6, 6.0034213e-6, 5.7797522e-6, 5.5687657e-6, 5.3773216e-6, 5.2238184e-6, 5.0756828e-6, 4.9383156e-6, 4.8214512e-6, 4.7223650e-6, 4.6294871e-6, 4.5718841e-6, 4.5342638e-6, 4.5007520e-6, 4.4839047e-6, 4.4869367e-6, 4.4982789e-6, 4.5450330e-6, 4.5977541e-6, 4.6718137e-6, 4.7661079e-6, 4.8701525e-6, 4.9939852e-6, 5.1354655e-6, 5.2846376e-6, 5.4684465e-6, 5.6623492e-6, 5.8730989e-6, 6.0903744e-6, 6.3420036e-6, 6.5968891e-6, 6.8710820e-6, 7.1723222e-6, 7.4866891e-6, 7.8198286e-6, 8.1658921e-6, 8.5241509e-6, 8.9034082e-6, 9.3249988e-6, 9.7246164e-6, 1.0161853e-5, 1.0601092e-5, 1.1090916e-5, 1.1568922e-5, 1.2067291e-5, 1.2591296e-5, 1.3116574e-5, 1.3680615e-5, 1.4255445e-5, 1.4833252e-5, 1.4833252e-5]),
                    np.array([6.2465996e-6, 6.0034213e-6, 5.7797522e-6, 5.5687657e-6, 5.3773216e-6, 5.2238184e-6, 5.0756828e-6, 4.9383156e-6, 4.8214512e-6, 4.7223650e-6, 4.6294871e-6, 4.5718841e-6, 4.5342638e-6, 4.5007520e-6, 4.4839047e-6, 4.4869367e-6, 4.4982789e-6, 4.5450330e-6, 4.5977541e-6, 4.6718137e-6, 4.7661079e-6, 4.8701525e-6, 4.9939852e-6, 5.1354655e-6, 5.2846376e-6, 5.4684465e-6, 5.6623492e-6, 5.8730989e-6, 6.0903744e-6, 6.3420036e-6, 6.5968891e-6, 6.8710820e-6, 7.1723222e-6, 7.4866891e-6, 7.8198286e-6, 8.1658921e-6, 8.5241509e-6, 8.9034082e-6, 9.3249988e-6, 9.7246164e-6, 1.0161853e-5, 1.0601092e-5, 1.1090916e-5, 1.1568922e-5, 1.2067291e-5, 1.2591296e-5, 1.3116574e-5, 1.3680615e-5, 1.4255445e-5, 1.4833252e-5, 1.4833252e-5]),
                    np.array([6.2465996e-6, 6.0034213e-6, 5.7797522e-6, 5.5687657e-6, 5.3773216e-6, 5.2238184e-6, 5.0756828e-6, 4.9383156e-6, 4.8214512e-6, 4.7223650e-6, 4.6294871e-6, 4.5718841e-6, 4.5342638e-6, 4.5007520e-6, 4.4839047e-6, 4.4869367e-6, 4.4982789e-6, 4.5450330e-6, 4.5977541e-6, 4.6718137e-6, 4.7661079e-6, 4.8701525e-6, 4.9939852e-6, 5.1354655e-6, 5.2846376e-6, 5.4684465e-6, 5.6623492e-6, 5.8730989e-6, 6.0903744e-6, 6.3420036e-6, 6.5968891e-6, 6.8710820e-6, 7.1723222e-6, 7.4866891e-6, 7.8198286e-6, 8.1658921e-6, 8.5241509e-6, 8.9034082e-6, 9.3249988e-6, 9.7246164e-6, 1.0161853e-5, 1.0601092e-5, 1.1090916e-5, 1.1568922e-5, 1.2067291e-5, 1.2591296e-5, 1.3116574e-5, 1.3680615e-5, 1.4255445e-5, 1.4833252e-5, 1.4833252e-5]),
                ),
        ],
        "qcd_y"    : np.array([6.2465996e-6, 6.0034213e-6, 5.7797522e-6, 5.5687657e-6, 5.3773216e-6, 5.2238184e-6, 5.0756828e-6, 4.9383156e-6, 4.8214512e-6, 4.7223650e-6, 4.6294871e-6, 4.5718841e-6, 4.5342638e-6, 4.5007520e-6, 4.4839047e-6, 4.4869367e-6, 4.4982789e-6, 4.5450330e-6, 4.5977541e-6, 4.6718137e-6, 4.7661079e-6, 4.8701525e-6, 4.9939852e-6, 5.1354655e-6, 5.2846376e-6, 5.4684465e-6, 5.6623492e-6, 5.8730989e-6, 6.0903744e-6, 6.3420036e-6, 6.5968891e-6, 6.8710820e-6, 7.1723222e-6, 7.4866891e-6, 7.8198286e-6, 8.1658921e-6, 8.5241509e-6, 8.9034082e-6, 9.3249988e-6, 9.7246164e-6, 1.0161853e-5, 1.0601092e-5, 1.1090916e-5, 1.1568922e-5, 1.2067291e-5, 1.2591296e-5, 1.3116574e-5, 1.3680615e-5, 1.4255445e-5, 1.4833252e-5, 1.4833252e-5]),
        "qcd_min"  : np.array([5.6866684e-6, 5.4651977e-6, 5.2613729e-6, 5.0689158e-6, 4.8941565e-6, 4.7537750e-6, 4.6182119e-6, 4.4922956e-6, 4.3849538e-6, 4.2937032e-6, 4.2079310e-6, 4.1541756e-6, 4.1184851e-6, 4.0863547e-6, 4.0694365e-6, 4.0703272e-6, 4.0787102e-6, 4.1192025e-6, 4.1650808e-6, 4.2302802e-6, 4.3136971e-6, 4.4057934e-6, 4.5160067e-6, 4.6419912e-6, 4.7748821e-6, 4.9391744e-6, 5.1125904e-6, 5.3012260e-6, 5.4957407e-6, 5.7213184e-6, 5.9497430e-6, 6.1957859e-6, 6.4659428e-6, 6.7482393e-6, 7.0474515e-6, 7.3582580e-6, 7.6802270e-6, 8.0209638e-6, 8.3999467e-6, 8.7593697e-6, 9.1525508e-6, 9.5474088e-6, 9.9882677e-6, 1.0418338e-5, 1.0866860e-5, 1.1338291e-5, 1.1811156e-5, 1.2318936e-5, 1.2836289e-5, 1.3356626e-5, 1.3356626e-5]),
        "qcd_max"  : np.array([6.9047938e-6, 6.6361078e-6, 6.3891395e-6, 6.1564219e-6, 5.9454311e-6, 5.7765647e-6, 5.6137220e-6, 5.4629977e-6, 5.3350625e-6, 5.2268797e-6, 5.1258012e-6, 5.0638247e-6, 5.0241005e-6, 4.9891562e-6, 4.9725786e-6, 4.9783504e-6, 4.9933967e-6, 5.0477573e-6, 5.1087655e-6, 5.1935058e-6, 5.3008631e-6, 5.4192638e-6, 5.5593981e-6, 5.7194262e-6, 5.8880825e-6, 6.0951841e-6, 6.3135340e-6, 6.5506543e-6, 6.7950709e-6, 7.0777272e-6, 7.3641504e-6, 7.6718472e-6, 8.0101063e-6, 8.3626263e-6, 8.7361286e-6, 9.1241232e-6, 9.5255295e-6, 9.9506079e-6, 1.0422845e-5, 1.0870221e-5, 1.1359812e-5, 1.1851819e-5, 1.2399784e-5, 1.2934737e-5, 1.3492307e-5, 1.4078770e-5, 1.4666284e-5, 1.5297122e-5, 1.5940219e-5, 1.6586246e-5, 1.6586246e-5]),
        "y"        : np.array([6.2465996e-6, 6.0034213e-6, 5.7797522e-6, 5.5687657e-6, 5.3773216e-6, 5.2238184e-6, 5.0756828e-6, 4.9383156e-6, 4.8214512e-6, 4.7223650e-6, 4.6294871e-6, 4.5718841e-6, 4.5342638e-6, 4.5007520e-6, 4.4839047e-6, 4.4869367e-6, 4.4982789e-6, 4.5450330e-6, 4.5977541e-6, 4.6718137e-6, 4.7661079e-6, 4.8701525e-6, 4.9939852e-6, 5.1354655e-6, 5.2846376e-6, 5.4684465e-6, 5.6623492e-6, 5.8730989e-6, 6.0903744e-6, 6.3420036e-6, 6.5968891e-6, 6.8710820e-6, 7.1723222e-6, 7.4866891e-6, 7.8198286e-6, 8.1658921e-6, 8.5241509e-6, 8.9034082e-6, 9.3249988e-6, 9.7246164e-6, 1.0161853e-5, 1.0601092e-5, 1.1090916e-5, 1.1568922e-5, 1.2067291e-5, 1.2591296e-5, 1.3116574e-5, 1.3680615e-5, 1.4255445e-5, 1.4833252e-5, 1.4833252e-5]),
        "ymin"     : np.array([5.6866684e-6, 5.4651977e-6, 5.2613729e-6, 5.0689158e-6, 4.8941565e-6, 4.7537750e-6, 4.6182119e-6, 4.4922956e-6, 4.3849538e-6, 4.2937032e-6, 4.2079310e-6, 4.1541756e-6, 4.1184851e-6, 4.0863547e-6, 4.0694365e-6, 4.0703272e-6, 4.0787102e-6, 4.1192025e-6, 4.1650808e-6, 4.2302802e-6, 4.3136971e-6, 4.4057934e-6, 4.5160067e-6, 4.6419912e-6, 4.7748821e-6, 4.9391744e-6, 5.1125904e-6, 5.3012260e-6, 5.4957407e-6, 5.7213184e-6, 5.9497430e-6, 6.1957859e-6, 6.4659428e-6, 6.7482393e-6, 7.0474515e-6, 7.3582580e-6, 7.6802270e-6, 8.0209638e-6, 8.3999467e-6, 8.7593697e-6, 9.1525508e-6, 9.5474088e-6, 9.9882677e-6, 1.0418338e-5, 1.0866860e-5, 1.1338291e-5, 1.1811156e-5, 1.2318936e-5, 1.2836289e-5, 1.3356626e-5, 1.3356626e-5]),
        "ymax"     : np.array([6.9047938e-6, 6.6361078e-6, 6.3891395e-6, 6.1564219e-6, 5.9454311e-6, 5.7765647e-6, 5.6137220e-6, 5.4629977e-6, 5.3350625e-6, 5.2268797e-6, 5.1258012e-6, 5.0638247e-6, 5.0241005e-6, 4.9891562e-6, 4.9725786e-6, 4.9783504e-6, 4.9933967e-6, 5.0477573e-6, 5.1087655e-6, 5.1935058e-6, 5.3008631e-6, 5.4192638e-6, 5.5593981e-6, 5.7194262e-6, 5.8880825e-6, 6.0951841e-6, 6.3135340e-6, 6.5506543e-6, 6.7950709e-6, 7.0777272e-6, 7.3641504e-6, 7.6718472e-6, 8.0101063e-6, 8.3626263e-6, 8.7361286e-6, 9.1241232e-6, 9.5255295e-6, 9.9506079e-6, 1.0422845e-5, 1.0870221e-5, 1.1359812e-5, 1.1851819e-5, 1.2399784e-5, 1.2934737e-5, 1.3492307e-5, 1.4078770e-5, 1.4666284e-5, 1.5297122e-5, 1.5940219e-5, 1.6586246e-5, 1.6586246e-5]),
        "channels" : [
            (r"$\mathrm{u}\bar{\mathrm{u}} + \mathrm{c}\bar{\mathrm{c}}$", np.array([4.6777566e-6, 4.4963962e-6, 4.3306007e-6, 4.1754689e-6, 4.0363593e-6, 3.9179347e-6, 3.8165473e-6, 3.7197417e-6, 3.6346919e-6, 3.5751766e-6, 3.5122544e-6, 3.4767729e-6, 3.4554987e-6, 3.4422327e-6, 3.4417328e-6, 3.4528854e-6, 3.4709891e-6, 3.5177672e-6, 3.5749797e-6, 3.6467432e-6, 3.7308688e-6, 3.8237671e-6, 3.9340001e-6, 4.0621100e-6, 4.1911104e-6, 4.3449958e-6, 4.5162833e-6, 4.6918556e-6, 4.8735482e-6, 5.0854161e-6, 5.2989079e-6, 5.5277044e-6, 5.7758933e-6, 6.0409618e-6, 6.3160095e-6, 6.5984170e-6, 6.8992288e-6, 7.2113545e-6, 7.5534076e-6, 7.8855627e-6, 8.2445573e-6, 8.6012535e-6, 9.0054786e-6, 9.3922951e-6, 9.8034487e-6, 1.0231410e-5, 1.0657907e-5, 1.1117212e-5, 1.1580680e-5, 1.2054302e-5, 1.2054302e-5])),
            (r"$\mathrm{d}\bar{\mathrm{d}} + \mathrm{s}\bar{\mathrm{s}} + \mathrm{b}\bar{\mathrm{b}}$", np.array([1.5688431e-6, 1.5070251e-6, 1.4491516e-6, 1.3932968e-6, 1.3409623e-6, 1.3058837e-6, 1.2591355e-6, 1.2185739e-6, 1.1867593e-6, 1.1471884e-6, 1.1172327e-6, 1.0951112e-6, 1.0787651e-6, 1.0585193e-6, 1.0421720e-6, 1.0340513e-6, 1.0272899e-6, 1.0272658e-6, 1.0227744e-6, 1.0250705e-6, 1.0352391e-6, 1.0463853e-6, 1.0599851e-6, 1.0733554e-6, 1.0935272e-6, 1.1234508e-6, 1.1460659e-6, 1.1812434e-6, 1.2168262e-6, 1.2565876e-6, 1.2979811e-6, 1.3433775e-6, 1.3964289e-6, 1.4457273e-6, 1.5038190e-6, 1.5674750e-6, 1.6249222e-6, 1.6920536e-6, 1.7715912e-6, 1.8390536e-6, 1.9172958e-6, 1.9998384e-6, 2.0854378e-6, 2.1766274e-6, 2.2638427e-6, 2.3598858e-6, 2.4586667e-6, 2.5634024e-6, 2.6747642e-6, 2.7789500e-6, 2.7789500e-6]))
        ],
    },
    {
        "slice_label"    : r"$\SI{5000}{GeV} < M_{\ell\bar{\ell}} < \SI{8000}{GeV}$",
        "x"        : np.array([-1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1]),
        "mid"      : np.array([-0.98, -0.94, -0.9, -0.86, -0.8200000000000001, -0.78, -0.74, -0.7, -0.66, -0.62, -0.5800000000000001, -0.54, -0.5, -0.45999999999999996, -0.42000000000000004, -0.38, -0.33999999999999997, -0.30000000000000004, -0.26, -0.22, -0.18, -0.14, -0.1, -0.06, -0.02, 0.02, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30000000000000004, 0.33999999999999997, 0.38, 0.42000000000000004, 0.45999999999999996, 0.5, 0.54, 0.5800000000000001, 0.62, 0.66, 0.7, 0.74, 0.78, 0.8200000000000001, 0.86, 0.9, 0.94, 0.98]),
        "pdf_results" : [
            (
                    r"NNPDF40\_nnlo\_as\_01180",
                    np.array([1.3268262e-7, 1.2786453e-7, 1.2307005e-7, 1.1780526e-7, 1.1356356e-7, 1.0892100e-7, 1.0505353e-7, 1.0109804e-7, 9.7278368e-8, 9.3859898e-8, 9.0732778e-8, 8.7265819e-8, 8.4585006e-8, 8.2016912e-8, 7.9708423e-8, 7.7419509e-8, 7.5101943e-8, 7.3713994e-8, 7.2071670e-8, 7.0690848e-8, 6.9458697e-8, 6.8574148e-8, 6.7759607e-8, 6.7352353e-8, 6.7120091e-8, 6.7021996e-8, 6.7173994e-8, 6.7534167e-8, 6.7992384e-8, 6.8908159e-8, 6.9745535e-8, 7.1108396e-8, 7.2621295e-8, 7.4232542e-8, 7.6019061e-8, 7.8080590e-8, 8.0376525e-8, 8.2963179e-8, 8.5596531e-8, 8.8529203e-8, 9.1685084e-8, 9.5096261e-8, 9.8559654e-8, 1.0240826e-7, 1.0639069e-7, 1.1076525e-7, 1.1520479e-7, 1.1965034e-7, 1.2456911e-7, 1.2975301e-7, 1.2975301e-7]),
                    np.array([1.3268262e-7, 1.2786453e-7, 1.2307005e-7, 1.1780526e-7, 1.1356356e-7, 1.0892100e-7, 1.0505353e-7, 1.0109804e-7, 9.7278368e-8, 9.3859898e-8, 9.0732778e-8, 8.7265819e-8, 8.4585006e-8, 8.2016912e-8, 7.9708423e-8, 7.7419509e-8, 7.5101943e-8, 7.3713994e-8, 7.2071670e-8, 7.0690848e-8, 6.9458697e-8, 6.8574148e-8, 6.7759607e-8, 6.7352353e-8, 6.7120091e-8, 6.7021996e-8, 6.7173994e-8, 6.7534167e-8, 6.7992384e-8, 6.8908159e-8, 6.9745535e-8, 7.1108396e-8, 7.2621295e-8, 7.4232542e-8, 7.6019061e-8, 7.8080590e-8, 8.0376525e-8, 8.2963179e-8, 8.5596531e-8, 8.8529203e-8, 9.1685084e-8, 9.5096261e-8, 9.8559654e-8, 1.0240826e-7, 1.0639069e-7, 1.1076525e-7, 1.1520479e-7, 1.1965034e-7, 1.2456911e-7, 1.2975301e-7, 1.2975301e-7]),
                    np.array([1.3268262e-7, 1.2786453e-7, 1.2307005e-7, 1.1780526e-7, 1.1356356e-7, 1.0892100e-7, 1.0505353e-7, 1.0109804e-7, 9.7278368e-8, 9.3859898e-8, 9.0732778e-8, 8.7265819e-8, 8.4585006e-8, 8.2016912e-8, 7.9708423e-8, 7.7419509e-8, 7.5101943e-8, 7.3713994e-8, 7.2071670e-8, 7.0690848e-8, 6.9458697e-8, 6.8574148e-8, 6.7759607e-8, 6.7352353e-8, 6.7120091e-8, 6.7021996e-8, 6.7173994e-8, 6.7534167e-8, 6.7992384e-8, 6.8908159e-8, 6.9745535e-8, 7.1108396e-8, 7.2621295e-8, 7.4232542e-8, 7.6019061e-8, 7.8080590e-8, 8.0376525e-8, 8.2963179e-8, 8.5596531e-8, 8.8529203e-8, 9.1685084e-8, 9.5096261e-8, 9.8559654e-8, 1.0240826e-7, 1.0639069e-7, 1.1076525e-7, 1.1520479e-7, 1.1965034e-7, 1.2456911e-7, 1.2975301e-7, 1.2975301e-7]),
                ),
        ],
        "qcd_y"    : np.array([1.3268262e-7, 1.2786453e-7, 1.2307005e-7, 1.1780526e-7, 1.1356356e-7, 1.0892100e-7, 1.0505353e-7, 1.0109804e-7, 9.7278368e-8, 9.3859898e-8, 9.0732778e-8, 8.7265819e-8, 8.4585006e-8, 8.2016912e-8, 7.9708423e-8, 7.7419509e-8, 7.5101943e-8, 7.3713994e-8, 7.2071670e-8, 7.0690848e-8, 6.9458697e-8, 6.8574148e-8, 6.7759607e-8, 6.7352353e-8, 6.7120091e-8, 6.7021996e-8, 6.7173994e-8, 6.7534167e-8, 6.7992384e-8, 6.8908159e-8, 6.9745535e-8, 7.1108396e-8, 7.2621295e-8, 7.4232542e-8, 7.6019061e-8, 7.8080590e-8, 8.0376525e-8, 8.2963179e-8, 8.5596531e-8, 8.8529203e-8, 9.1685084e-8, 9.5096261e-8, 9.8559654e-8, 1.0240826e-7, 1.0639069e-7, 1.1076525e-7, 1.1520479e-7, 1.1965034e-7, 1.2456911e-7, 1.2975301e-7, 1.2975301e-7]),
        "qcd_min"  : np.array([1.1796105e-7, 1.1367206e-7, 1.0940640e-7, 1.0473022e-7, 1.0094995e-7, 9.6817367e-8, 9.3368040e-8, 8.9841440e-8, 8.6438729e-8, 8.3389429e-8, 8.0592033e-8, 7.7501048e-8, 7.5100969e-8, 7.2799905e-8, 7.0734091e-8, 6.8678517e-8, 6.6600127e-8, 6.5343484e-8, 6.3861008e-8, 6.2611336e-8, 6.1492668e-8, 6.0680109e-8, 5.9926419e-8, 5.9534539e-8, 5.9300399e-8, 5.9178721e-8, 5.9281434e-8, 5.9570740e-8, 5.9942127e-8, 6.0724158e-8, 6.1428979e-8, 6.2599969e-8, 6.3904432e-8, 6.5297162e-8, 6.6843778e-8, 6.8632843e-8, 7.0630564e-8, 7.2883947e-8, 7.5177425e-8, 7.7737465e-8, 8.0492081e-8, 8.3474223e-8, 8.6500375e-8, 8.9868353e-8, 9.3354036e-8, 9.7186902e-8, 1.0107609e-7, 1.0497339e-7, 1.0928445e-7, 1.1383240e-7, 1.1383240e-7]),
        "qcd_max"  : np.array([1.5031713e-7, 1.4486644e-7, 1.3943889e-7, 1.3346815e-7, 1.2867540e-7, 1.2342279e-7, 1.1905669e-7, 1.1458942e-7, 1.1027138e-7, 1.0641247e-7, 1.0289371e-7, 9.8977949e-8, 9.5964016e-8, 9.3079140e-8, 9.0482540e-8, 8.7917761e-8, 8.5317122e-8, 8.3776037e-8, 8.1946300e-8, 8.0412453e-8, 7.9048531e-8, 7.8082511e-8, 7.7200627e-8, 7.6780539e-8, 7.6555496e-8, 7.6492081e-8, 7.6709085e-8, 7.7159748e-8, 7.7728575e-8, 7.8810410e-8, 7.9813916e-8, 8.1414088e-8, 8.3184083e-8, 8.5064404e-8, 8.7145967e-8, 8.9541837e-8, 9.2203014e-8, 9.5197361e-8, 9.8246621e-8, 1.0163432e-7, 1.0528027e-7, 1.0921464e-7, 1.1321161e-7, 1.1764575e-7, 1.2223338e-7, 1.2726711e-7, 1.3237650e-7, 1.3748886e-7, 1.4314700e-7, 1.4910383e-7, 1.4910383e-7]),
        "y"        : np.array([1.3268262e-7, 1.2786453e-7, 1.2307005e-7, 1.1780526e-7, 1.1356356e-7, 1.0892100e-7, 1.0505353e-7, 1.0109804e-7, 9.7278368e-8, 9.3859898e-8, 9.0732778e-8, 8.7265819e-8, 8.4585006e-8, 8.2016912e-8, 7.9708423e-8, 7.7419509e-8, 7.5101943e-8, 7.3713994e-8, 7.2071670e-8, 7.0690848e-8, 6.9458697e-8, 6.8574148e-8, 6.7759607e-8, 6.7352353e-8, 6.7120091e-8, 6.7021996e-8, 6.7173994e-8, 6.7534167e-8, 6.7992384e-8, 6.8908159e-8, 6.9745535e-8, 7.1108396e-8, 7.2621295e-8, 7.4232542e-8, 7.6019061e-8, 7.8080590e-8, 8.0376525e-8, 8.2963179e-8, 8.5596531e-8, 8.8529203e-8, 9.1685084e-8, 9.5096261e-8, 9.8559654e-8, 1.0240826e-7, 1.0639069e-7, 1.1076525e-7, 1.1520479e-7, 1.1965034e-7, 1.2456911e-7, 1.2975301e-7, 1.2975301e-7]),
        "ymin"     : np.array([1.1796105e-7, 1.1367206e-7, 1.0940640e-7, 1.0473022e-7, 1.0094995e-7, 9.6817367e-8, 9.3368040e-8, 8.9841440e-8, 8.6438729e-8, 8.3389429e-8, 8.0592033e-8, 7.7501048e-8, 7.5100969e-8, 7.2799905e-8, 7.0734091e-8, 6.8678517e-8, 6.6600127e-8, 6.5343484e-8, 6.3861008e-8, 6.2611336e-8, 6.1492668e-8, 6.0680109e-8, 5.9926419e-8, 5.9534539e-8, 5.9300399e-8, 5.9178721e-8, 5.9281434e-8, 5.9570740e-8, 5.9942127e-8, 6.0724158e-8, 6.1428979e-8, 6.2599969e-8, 6.3904432e-8, 6.5297162e-8, 6.6843778e-8, 6.8632843e-8, 7.0630564e-8, 7.2883947e-8, 7.5177425e-8, 7.7737465e-8, 8.0492081e-8, 8.3474223e-8, 8.6500375e-8, 8.9868353e-8, 9.3354036e-8, 9.7186902e-8, 1.0107609e-7, 1.0497339e-7, 1.0928445e-7, 1.1383240e-7, 1.1383240e-7]),
        "ymax"     : np.array([1.5031713e-7, 1.4486644e-7, 1.3943889e-7, 1.3346815e-7, 1.2867540e-7, 1.2342279e-7, 1.1905669e-7, 1.1458942e-7, 1.1027138e-7, 1.0641247e-7, 1.0289371e-7, 9.8977949e-8, 9.5964016e-8, 9.3079140e-8, 9.0482540e-8, 8.7917761e-8, 8.5317122e-8, 8.3776037e-8, 8.1946300e-8, 8.0412453e-8, 7.9048531e-8, 7.8082511e-8, 7.7200627e-8, 7.6780539e-8, 7.6555496e-8, 7.6492081e-8, 7.6709085e-8, 7.7159748e-8, 7.7728575e-8, 7.8810410e-8, 7.9813916e-8, 8.1414088e-8, 8.3184083e-8, 8.5064404e-8, 8.7145967e-8, 8.9541837e-8, 9.2203014e-8, 9.5197361e-8, 9.8246621e-8, 1.0163432e-7, 1.0528027e-7, 1.0921464e-7, 1.1321161e-7, 1.1764575e-7, 1.2223338e-7, 1.2726711e-7, 1.3237650e-7, 1.3748886e-7, 1.4314700e-7, 1.4910383e-7, 1.4910383e-7]),
        "channels" : [
            (r"$\mathrm{u}\bar{\mathrm{u}} + \mathrm{c}\bar{\mathrm{c}}$", np.array([1.0634997e-7, 1.0232840e-7, 9.8488999e-8, 9.4336233e-8, 9.0829972e-8, 8.7163111e-8, 8.4152495e-8, 8.1031902e-8, 7.8203491e-8, 7.5443983e-8, 7.3006568e-8, 7.0236490e-8, 6.8148391e-8, 6.6168729e-8, 6.4386166e-8, 6.2647603e-8, 6.0968245e-8, 5.9862611e-8, 5.8622760e-8, 5.7758449e-8, 5.6855750e-8, 5.6212654e-8, 5.5703319e-8, 5.5430259e-8, 5.5446303e-8, 5.5542422e-8, 5.5708682e-8, 5.6224214e-8, 5.6730599e-8, 5.7643506e-8, 5.8465496e-8, 5.9752070e-8, 6.1119329e-8, 6.2593741e-8, 6.4180037e-8, 6.6108169e-8, 6.8114131e-8, 7.0406260e-8, 7.2716285e-8, 7.5307416e-8, 7.8063327e-8, 8.0943969e-8, 8.4008157e-8, 8.7317120e-8, 9.0784410e-8, 9.4583904e-8, 9.8309896e-8, 1.0224158e-7, 1.0641016e-7, 1.1078435e-7, 1.1078435e-7])),
            (r"$\mathrm{d}\bar{\mathrm{d}} + \mathrm{s}\bar{\mathrm{s}} + \mathrm{b}\bar{\mathrm{b}}$", np.array([2.6332651e-8, 2.5536122e-8, 2.4581048e-8, 2.3469032e-8, 2.2733585e-8, 2.1757885e-8, 2.0901034e-8, 2.0066140e-8, 1.9074876e-8, 1.8415915e-8, 1.7726211e-8, 1.7029329e-8, 1.6436615e-8, 1.5848183e-8, 1.5322258e-8, 1.4771906e-8, 1.4133698e-8, 1.3851383e-8, 1.3448911e-8, 1.2932399e-8, 1.2602947e-8, 1.2361494e-8, 1.2056288e-8, 1.1922094e-8, 1.1673787e-8, 1.1479574e-8, 1.1465312e-8, 1.1309954e-8, 1.1261785e-8, 1.1264654e-8, 1.1280039e-8, 1.1356326e-8, 1.1501966e-8, 1.1638801e-8, 1.1839024e-8, 1.1972421e-8, 1.2262394e-8, 1.2556919e-8, 1.2880247e-8, 1.3221787e-8, 1.3621758e-8, 1.4152292e-8, 1.4551497e-8, 1.5091138e-8, 1.5606276e-8, 1.6181342e-8, 1.6894894e-8, 1.7408752e-8, 1.8158948e-8, 1.8968662e-8, 1.8968662e-8]))
        ],
    },
    {
        "slice_label"    : r"$\SI{8000}{GeV} < M_{\ell\bar{\ell}} < \SI{10000}{GeV}$",
        "x"        : np.array([-1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1]),
        "mid"      : np.array([-0.98, -0.94, -0.9, -0.86, -0.8200000000000001, -0.78, -0.74, -0.7, -0.66, -0.62, -0.5800000000000001, -0.54, -0.5, -0.45999999999999996, -0.42000000000000004, -0.38, -0.33999999999999997, -0.30000000000000004, -0.26, -0.22, -0.18, -0.14, -0.1, -0.06, -0.02, 0.02, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30000000000000004, 0.33999999999999997, 0.38, 0.42000000000000004, 0.45999999999999996, 0.5, 0.54, 0.5800000000000001, 0.62, 0.66, 0.7, 0.74, 0.78, 0.8200000000000001, 0.86, 0.9, 0.94, 0.98]),
        "pdf_results" : [
            (
                    r"NNPDF40\_nnlo\_as\_01180",
                    np.array([7.1861026e-10, 6.9052982e-10, 6.6277085e-10, 6.3573468e-10, 6.1024606e-10, 5.8475908e-10, 5.6096402e-10, 5.3833955e-10, 5.1558216e-10, 4.9404513e-10, 4.7398312e-10, 4.5273929e-10, 4.3405736e-10, 4.1643201e-10, 3.9946960e-10, 3.8248920e-10, 3.6635833e-10, 3.5323838e-10, 3.3921398e-10, 3.2684546e-10, 3.1412226e-10, 3.0313435e-10, 2.9192112e-10, 2.8290167e-10, 2.7525858e-10, 2.6743219e-10, 2.5992966e-10, 2.5408868e-10, 2.4832543e-10, 2.4538556e-10, 2.4100708e-10, 2.3910018e-10, 2.3747322e-10, 2.3728921e-10, 2.3669027e-10, 2.3798116e-10, 2.3947283e-10, 2.4222947e-10, 2.4587694e-10, 2.5017261e-10, 2.5579842e-10, 2.6167570e-10, 2.6866230e-10, 2.7632565e-10, 2.8548191e-10, 2.9544050e-10, 3.0598398e-10, 3.1671560e-10, 3.2904259e-10, 3.4234846e-10, 3.4234846e-10]),
                    np.array([7.1861026e-10, 6.9052982e-10, 6.6277085e-10, 6.3573468e-10, 6.1024606e-10, 5.8475908e-10, 5.6096402e-10, 5.3833955e-10, 5.1558216e-10, 4.9404513e-10, 4.7398312e-10, 4.5273929e-10, 4.3405736e-10, 4.1643201e-10, 3.9946960e-10, 3.8248920e-10, 3.6635833e-10, 3.5323838e-10, 3.3921398e-10, 3.2684546e-10, 3.1412226e-10, 3.0313435e-10, 2.9192112e-10, 2.8290167e-10, 2.7525858e-10, 2.6743219e-10, 2.5992966e-10, 2.5408868e-10, 2.4832543e-10, 2.4538556e-10, 2.4100708e-10, 2.3910018e-10, 2.3747322e-10, 2.3728921e-10, 2.3669027e-10, 2.3798116e-10, 2.3947283e-10, 2.4222947e-10, 2.4587694e-10, 2.5017261e-10, 2.5579842e-10, 2.6167570e-10, 2.6866230e-10, 2.7632565e-10, 2.8548191e-10, 2.9544050e-10, 3.0598398e-10, 3.1671560e-10, 3.2904259e-10, 3.4234846e-10, 3.4234846e-10]),
                    np.array([7.1861026e-10, 6.9052982e-10, 6.6277085e-10, 6.3573468e-10, 6.1024606e-10, 5.8475908e-10, 5.6096402e-10, 5.3833955e-10, 5.1558216e-10, 4.9404513e-10, 4.7398312e-10, 4.5273929e-10, 4.3405736e-10, 4.1643201e-10, 3.9946960e-10, 3.8248920e-10, 3.6635833e-10, 3.5323838e-10, 3.3921398e-10, 3.2684546e-10, 3.1412226e-10, 3.0313435e-10, 2.9192112e-10, 2.8290167e-10, 2.7525858e-10, 2.6743219e-10, 2.5992966e-10, 2.5408868e-10, 2.4832543e-10, 2.4538556e-10, 2.4100708e-10, 2.3910018e-10, 2.3747322e-10, 2.3728921e-10, 2.3669027e-10, 2.3798116e-10, 2.3947283e-10, 2.4222947e-10, 2.4587694e-10, 2.5017261e-10, 2.5579842e-10, 2.6167570e-10, 2.6866230e-10, 2.7632565e-10, 2.8548191e-10, 2.9544050e-10, 3.0598398e-10, 3.1671560e-10, 3.2904259e-10, 3.4234846e-10, 3.4234846e-10]),
                ),
        ],
        "qcd_y"    : np.array([7.1861026e-10, 6.9052982e-10, 6.6277085e-10, 6.3573468e-10, 6.1024606e-10, 5.8475908e-10, 5.6096402e-10, 5.3833955e-10, 5.1558216e-10, 4.9404513e-10, 4.7398312e-10, 4.5273929e-10, 4.3405736e-10, 4.1643201e-10, 3.9946960e-10, 3.8248920e-10, 3.6635833e-10, 3.5323838e-10, 3.3921398e-10, 3.2684546e-10, 3.1412226e-10, 3.0313435e-10, 2.9192112e-10, 2.8290167e-10, 2.7525858e-10, 2.6743219e-10, 2.5992966e-10, 2.5408868e-10, 2.4832543e-10, 2.4538556e-10, 2.4100708e-10, 2.3910018e-10, 2.3747322e-10, 2.3728921e-10, 2.3669027e-10, 2.3798116e-10, 2.3947283e-10, 2.4222947e-10, 2.4587694e-10, 2.5017261e-10, 2.5579842e-10, 2.6167570e-10, 2.6866230e-10, 2.7632565e-10, 2.8548191e-10, 2.9544050e-10, 3.0598398e-10, 3.1671560e-10, 3.2904259e-10, 3.4234846e-10, 3.4234846e-10]),
        "qcd_min"  : np.array([6.0733393e-10, 5.8359975e-10, 5.6013663e-10, 5.3731953e-10, 5.1577550e-10, 4.9425471e-10, 4.7413668e-10, 4.5502613e-10, 4.3579683e-10, 4.1761844e-10, 4.0067486e-10, 3.8275700e-10, 3.6698642e-10, 3.5210349e-10, 3.3780462e-10, 3.2347174e-10, 3.0987916e-10, 2.9880336e-10, 2.8698310e-10, 2.7656554e-10, 2.6585543e-10, 2.5660126e-10, 2.4717194e-10, 2.3959131e-10, 2.3316588e-10, 2.2659350e-10, 2.2029631e-10, 2.1541466e-10, 2.1059049e-10, 2.0816229e-10, 2.0451475e-10, 2.0294117e-10, 2.0162673e-10, 2.0153357e-10, 2.0108215e-10, 2.0224072e-10, 2.0356183e-10, 2.0595804e-10, 2.0910726e-10, 2.1281117e-10, 2.1763190e-10, 2.2268082e-10, 2.2865687e-10, 2.3520471e-10, 2.4302664e-10, 2.5151226e-10, 2.6051000e-10, 2.6966237e-10, 2.8016281e-10, 2.9149775e-10, 2.9149775e-10]),
        "qcd_max"  : np.array([8.5956589e-10, 8.2598102e-10, 7.9278040e-10, 7.6039124e-10, 7.2990629e-10, 6.9938974e-10, 6.7093887e-10, 6.4386082e-10, 6.1663235e-10, 5.9083445e-10, 5.6681916e-10, 5.4135215e-10, 5.1897698e-10, 4.9787381e-10, 4.7752729e-10, 4.5718800e-10, 4.3783083e-10, 4.2211740e-10, 4.0529176e-10, 3.9044215e-10, 3.7515652e-10, 3.6196347e-10, 3.4847733e-10, 3.3762360e-10, 3.2842814e-10, 3.1900151e-10, 3.0995996e-10, 3.0288807e-10, 2.9592219e-10, 2.9231820e-10, 2.8699894e-10, 2.8465892e-10, 2.8261876e-10, 2.8230239e-10, 2.8150098e-10, 2.8294006e-10, 2.8463095e-10, 2.8782554e-10, 2.9208532e-10, 2.9710983e-10, 3.0373680e-10, 3.1064043e-10, 3.1888707e-10, 3.2794310e-10, 3.3876620e-10, 3.5057090e-10, 3.6304796e-10, 3.7575647e-10, 3.9037413e-10, 4.0615136e-10, 4.0615136e-10]),
        "y"        : np.array([7.1861026e-10, 6.9052982e-10, 6.6277085e-10, 6.3573468e-10, 6.1024606e-10, 5.8475908e-10, 5.6096402e-10, 5.3833955e-10, 5.1558216e-10, 4.9404513e-10, 4.7398312e-10, 4.5273929e-10, 4.3405736e-10, 4.1643201e-10, 3.9946960e-10, 3.8248920e-10, 3.6635833e-10, 3.5323838e-10, 3.3921398e-10, 3.2684546e-10, 3.1412226e-10, 3.0313435e-10, 2.9192112e-10, 2.8290167e-10, 2.7525858e-10, 2.6743219e-10, 2.5992966e-10, 2.5408868e-10, 2.4832543e-10, 2.4538556e-10, 2.4100708e-10, 2.3910018e-10, 2.3747322e-10, 2.3728921e-10, 2.3669027e-10, 2.3798116e-10, 2.3947283e-10, 2.4222947e-10, 2.4587694e-10, 2.5017261e-10, 2.5579842e-10, 2.6167570e-10, 2.6866230e-10, 2.7632565e-10, 2.8548191e-10, 2.9544050e-10, 3.0598398e-10, 3.1671560e-10, 3.2904259e-10, 3.4234846e-10, 3.4234846e-10]),
        "ymin"     : np.array([6.0733393e-10, 5.8359975e-10, 5.6013663e-10, 5.3731953e-10, 5.1577550e-10, 4.9425471e-10, 4.7413668e-10, 4.5502613e-10, 4.3579683e-10, 4.1761844e-10, 4.0067486e-10, 3.8275700e-10, 3.6698642e-10, 3.5210349e-10, 3.3780462e-10, 3.2347174e-10, 3.0987916e-10, 2.9880336e-10, 2.8698310e-10, 2.7656554e-10, 2.6585543e-10, 2.5660126e-10, 2.4717194e-10, 2.3959131e-10, 2.3316588e-10, 2.2659350e-10, 2.2029631e-10, 2.1541466e-10, 2.1059049e-10, 2.0816229e-10, 2.0451475e-10, 2.0294117e-10, 2.0162673e-10, 2.0153357e-10, 2.0108215e-10, 2.0224072e-10, 2.0356183e-10, 2.0595804e-10, 2.0910726e-10, 2.1281117e-10, 2.1763190e-10, 2.2268082e-10, 2.2865687e-10, 2.3520471e-10, 2.4302664e-10, 2.5151226e-10, 2.6051000e-10, 2.6966237e-10, 2.8016281e-10, 2.9149775e-10, 2.9149775e-10]),
        "ymax"     : np.array([8.5956589e-10, 8.2598102e-10, 7.9278040e-10, 7.6039124e-10, 7.2990629e-10, 6.9938974e-10, 6.7093887e-10, 6.4386082e-10, 6.1663235e-10, 5.9083445e-10, 5.6681916e-10, 5.4135215e-10, 5.1897698e-10, 4.9787381e-10, 4.7752729e-10, 4.5718800e-10, 4.3783083e-10, 4.2211740e-10, 4.0529176e-10, 3.9044215e-10, 3.7515652e-10, 3.6196347e-10, 3.4847733e-10, 3.3762360e-10, 3.2842814e-10, 3.1900151e-10, 3.0995996e-10, 3.0288807e-10, 2.9592219e-10, 2.9231820e-10, 2.8699894e-10, 2.8465892e-10, 2.8261876e-10, 2.8230239e-10, 2.8150098e-10, 2.8294006e-10, 2.8463095e-10, 2.8782554e-10, 2.9208532e-10, 2.9710983e-10, 3.0373680e-10, 3.1064043e-10, 3.1888707e-10, 3.2794310e-10, 3.3876620e-10, 3.5057090e-10, 3.6304796e-10, 3.7575647e-10, 3.9037413e-10, 4.0615136e-10, 4.0615136e-10]),
        "channels" : [
            (r"$\mathrm{u}\bar{\mathrm{u}} + \mathrm{c}\bar{\mathrm{c}}$", np.array([6.7156667e-10, 6.4480795e-10, 6.1906876e-10, 5.9405103e-10, 5.7003835e-10, 5.4634372e-10, 5.2393879e-10, 5.0270170e-10, 4.8196562e-10, 4.6168043e-10, 4.4301047e-10, 4.2289518e-10, 4.0551136e-10, 3.8900214e-10, 3.7319487e-10, 3.5702887e-10, 3.4216975e-10, 3.2985593e-10, 3.1665205e-10, 3.0540505e-10, 2.9336077e-10, 2.8289755e-10, 2.7262160e-10, 2.6375604e-10, 2.5684735e-10, 2.4957103e-10, 2.4239158e-10, 2.3696269e-10, 2.3139541e-10, 2.2872965e-10, 2.2472346e-10, 2.2285342e-10, 2.2131509e-10, 2.2098926e-10, 2.2028936e-10, 2.2170613e-10, 2.2288690e-10, 2.2546879e-10, 2.2889598e-10, 2.3280137e-10, 2.3804426e-10, 2.4346151e-10, 2.4972534e-10, 2.5696862e-10, 2.6540118e-10, 2.7484717e-10, 2.8450746e-10, 2.9460883e-10, 3.0601707e-10, 3.1827665e-10, 3.1827665e-10])),
            (r"$\mathrm{d}\bar{\mathrm{d}} + \mathrm{s}\bar{\mathrm{s}} + \mathrm{b}\bar{\mathrm{b}}$", np.array([4.7043587e-11, 4.5721878e-11, 4.3702086e-11, 4.1683641e-11, 4.0207715e-11, 3.8415354e-11, 3.7025234e-11, 3.5637851e-11, 3.3616539e-11, 3.2364703e-11, 3.0972646e-11, 2.9844116e-11, 2.8546008e-11, 2.7429871e-11, 2.6274730e-11, 2.5460333e-11, 2.4188586e-11, 2.3382450e-11, 2.2561926e-11, 2.1440412e-11, 2.0761496e-11, 2.0236798e-11, 1.9299518e-11, 1.9145639e-11, 1.8411239e-11, 1.7861161e-11, 1.7538087e-11, 1.7125990e-11, 1.6930016e-11, 1.6655912e-11, 1.6283622e-11, 1.6246755e-11, 1.6158123e-11, 1.6299957e-11, 1.6400918e-11, 1.6275026e-11, 1.6585928e-11, 1.6760680e-11, 1.6980956e-11, 1.7371237e-11, 1.7754155e-11, 1.8214193e-11, 1.8936959e-11, 1.9357026e-11, 2.0080730e-11, 2.0593338e-11, 2.1476522e-11, 2.2106766e-11, 2.3025520e-11, 2.4071814e-11, 2.4071814e-11]))
        ],
    },
    {
        "slice_label"    : r"$\SI{10000}{GeV} < M_{\ell\bar{\ell}} < \SI{14000}{GeV}$",
        "x"        : np.array([-1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1]),
        "mid"      : np.array([-0.98, -0.94, -0.9, -0.86, -0.8200000000000001, -0.78, -0.74, -0.7, -0.66, -0.62, -0.5800000000000001, -0.54, -0.5, -0.45999999999999996, -0.42000000000000004, -0.38, -0.33999999999999997, -0.30000000000000004, -0.26, -0.22, -0.18, -0.14, -0.1, -0.06, -0.02, 0.02, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30000000000000004, 0.33999999999999997, 0.38, 0.42000000000000004, 0.45999999999999996, 0.5, 0.54, 0.5800000000000001, 0.62, 0.66, 0.7, 0.74, 0.78, 0.8200000000000001, 0.86, 0.9, 0.94, 0.98]),
        "pdf_results" : [
            (
                    r"NNPDF40\_nnlo\_as\_01180",
                    np.array([5.9784203e-12, 5.7334000e-12, 5.5158405e-12, 5.2821781e-12, 5.0757322e-12, 4.8721048e-12, 4.6751045e-12, 4.4854392e-12, 4.3023311e-12, 4.1272521e-12, 3.9585369e-12, 3.8041358e-12, 3.6380746e-12, 3.5024271e-12, 3.3710447e-12, 3.2510444e-12, 3.1260172e-12, 3.0105698e-12, 2.9080579e-12, 2.8034115e-12, 2.7170918e-12, 2.6321534e-12, 2.5586050e-12, 2.4978602e-12, 2.4378311e-12, 2.3877748e-12, 2.3427105e-12, 2.3080580e-12, 2.2782058e-12, 2.2586106e-12, 2.2397268e-12, 2.2391422e-12, 2.2396377e-12, 2.2520626e-12, 2.2723190e-12, 2.2981583e-12, 2.3297045e-12, 2.3715628e-12, 2.4196963e-12, 2.4771198e-12, 2.5408710e-12, 2.6118028e-12, 2.6948071e-12, 2.7819795e-12, 2.8750032e-12, 2.9793716e-12, 3.0914562e-12, 3.2108107e-12, 3.3326464e-12, 3.4695631e-12, 3.4695631e-12]),
                    np.array([5.9784203e-12, 5.7334000e-12, 5.5158405e-12, 5.2821781e-12, 5.0757322e-12, 4.8721048e-12, 4.6751045e-12, 4.4854392e-12, 4.3023311e-12, 4.1272521e-12, 3.9585369e-12, 3.8041358e-12, 3.6380746e-12, 3.5024271e-12, 3.3710447e-12, 3.2510444e-12, 3.1260172e-12, 3.0105698e-12, 2.9080579e-12, 2.8034115e-12, 2.7170918e-12, 2.6321534e-12, 2.5586050e-12, 2.4978602e-12, 2.4378311e-12, 2.3877748e-12, 2.3427105e-12, 2.3080580e-12, 2.2782058e-12, 2.2586106e-12, 2.2397268e-12, 2.2391422e-12, 2.2396377e-12, 2.2520626e-12, 2.2723190e-12, 2.2981583e-12, 2.3297045e-12, 2.3715628e-12, 2.4196963e-12, 2.4771198e-12, 2.5408710e-12, 2.6118028e-12, 2.6948071e-12, 2.7819795e-12, 2.8750032e-12, 2.9793716e-12, 3.0914562e-12, 3.2108107e-12, 3.3326464e-12, 3.4695631e-12, 3.4695631e-12]),
                    np.array([5.9784203e-12, 5.7334000e-12, 5.5158405e-12, 5.2821781e-12, 5.0757322e-12, 4.8721048e-12, 4.6751045e-12, 4.4854392e-12, 4.3023311e-12, 4.1272521e-12, 3.9585369e-12, 3.8041358e-12, 3.6380746e-12, 3.5024271e-12, 3.3710447e-12, 3.2510444e-12, 3.1260172e-12, 3.0105698e-12, 2.9080579e-12, 2.8034115e-12, 2.7170918e-12, 2.6321534e-12, 2.5586050e-12, 2.4978602e-12, 2.4378311e-12, 2.3877748e-12, 2.3427105e-12, 2.3080580e-12, 2.2782058e-12, 2.2586106e-12, 2.2397268e-12, 2.2391422e-12, 2.2396377e-12, 2.2520626e-12, 2.2723190e-12, 2.2981583e-12, 2.3297045e-12, 2.3715628e-12, 2.4196963e-12, 2.4771198e-12, 2.5408710e-12, 2.6118028e-12, 2.6948071e-12, 2.7819795e-12, 2.8750032e-12, 2.9793716e-12, 3.0914562e-12, 3.2108107e-12, 3.3326464e-12, 3.4695631e-12, 3.4695631e-12]),
                ),
        ],
        "qcd_y"    : np.array([5.9784203e-12, 5.7334000e-12, 5.5158405e-12, 5.2821781e-12, 5.0757322e-12, 4.8721048e-12, 4.6751045e-12, 4.4854392e-12, 4.3023311e-12, 4.1272521e-12, 3.9585369e-12, 3.8041358e-12, 3.6380746e-12, 3.5024271e-12, 3.3710447e-12, 3.2510444e-12, 3.1260172e-12, 3.0105698e-12, 2.9080579e-12, 2.8034115e-12, 2.7170918e-12, 2.6321534e-12, 2.5586050e-12, 2.4978602e-12, 2.4378311e-12, 2.3877748e-12, 2.3427105e-12, 2.3080580e-12, 2.2782058e-12, 2.2586106e-12, 2.2397268e-12, 2.2391422e-12, 2.2396377e-12, 2.2520626e-12, 2.2723190e-12, 2.2981583e-12, 2.3297045e-12, 2.3715628e-12, 2.4196963e-12, 2.4771198e-12, 2.5408710e-12, 2.6118028e-12, 2.6948071e-12, 2.7819795e-12, 2.8750032e-12, 2.9793716e-12, 3.0914562e-12, 3.2108107e-12, 3.3326464e-12, 3.4695631e-12, 3.4695631e-12]),
        "qcd_min"  : np.array([4.8320667e-12, 4.6340168e-12, 4.4581332e-12, 4.2694034e-12, 4.1024778e-12, 3.9380551e-12, 3.7788147e-12, 3.6256865e-12, 3.4778027e-12, 3.3363436e-12, 3.2001398e-12, 3.0753736e-12, 2.9412981e-12, 2.8318282e-12, 2.7257422e-12, 2.6288731e-12, 2.5280116e-12, 2.4349649e-12, 2.3522879e-12, 2.2679037e-12, 2.1983550e-12, 2.1299028e-12, 2.0707980e-12, 2.0218582e-12, 1.9736388e-12, 1.9333927e-12, 1.8971892e-12, 1.8695481e-12, 1.8457255e-12, 1.8302086e-12, 1.8152970e-12, 1.8151223e-12, 1.8158688e-12, 1.8262736e-12, 1.8429841e-12, 1.8642947e-12, 1.8901548e-12, 1.9243658e-12, 1.9637075e-12, 2.0105556e-12, 2.0624878e-12, 2.1202545e-12, 2.1877949e-12, 2.2587192e-12, 2.3343787e-12, 2.4192297e-12, 2.5103115e-12, 2.6072867e-12, 2.7062549e-12, 2.8174517e-12, 2.8174517e-12]),
        "qcd_max"  : np.array([7.4982907e-12, 7.1910016e-12, 6.9182003e-12, 6.6249200e-12, 6.3661005e-12, 6.1104359e-12, 5.8633804e-12, 5.6252128e-12, 5.3953612e-12, 5.1756951e-12, 4.9638200e-12, 4.7701171e-12, 4.5616013e-12, 4.3911866e-12, 4.2262299e-12, 4.0755193e-12, 3.9183857e-12, 3.7731478e-12, 3.6442761e-12, 3.5126965e-12, 3.4040659e-12, 3.2972015e-12, 3.2043850e-12, 3.1279352e-12, 3.0521448e-12, 2.9890066e-12, 2.9321201e-12, 2.8880447e-12, 2.8500918e-12, 2.8249788e-12, 2.8007072e-12, 2.7994768e-12, 2.7995189e-12, 2.8144977e-12, 2.8393388e-12, 2.8710337e-12, 2.9099934e-12, 2.9618605e-12, 3.0214992e-12, 3.0927929e-12, 3.1720739e-12, 3.2603093e-12, 3.3636599e-12, 3.4722126e-12, 3.5880952e-12, 3.7181683e-12, 3.8579312e-12, 4.0067802e-12, 4.1587643e-12, 4.3295963e-12, 4.3295963e-12]),
        "y"        : np.array([5.9784203e-12, 5.7334000e-12, 5.5158405e-12, 5.2821781e-12, 5.0757322e-12, 4.8721048e-12, 4.6751045e-12, 4.4854392e-12, 4.3023311e-12, 4.1272521e-12, 3.9585369e-12, 3.8041358e-12, 3.6380746e-12, 3.5024271e-12, 3.3710447e-12, 3.2510444e-12, 3.1260172e-12, 3.0105698e-12, 2.9080579e-12, 2.8034115e-12, 2.7170918e-12, 2.6321534e-12, 2.5586050e-12, 2.4978602e-12, 2.4378311e-12, 2.3877748e-12, 2.3427105e-12, 2.3080580e-12, 2.2782058e-12, 2.2586106e-12, 2.2397268e-12, 2.2391422e-12, 2.2396377e-12, 2.2520626e-12, 2.2723190e-12, 2.2981583e-12, 2.3297045e-12, 2.3715628e-12, 2.4196963e-12, 2.4771198e-12, 2.5408710e-12, 2.6118028e-12, 2.6948071e-12, 2.7819795e-12, 2.8750032e-12, 2.9793716e-12, 3.0914562e-12, 3.2108107e-12, 3.3326464e-12, 3.4695631e-12, 3.4695631e-12]),
        "ymin"     : np.array([4.8320667e-12, 4.6340168e-12, 4.4581332e-12, 4.2694034e-12, 4.1024778e-12, 3.9380551e-12, 3.7788147e-12, 3.6256865e-12, 3.4778027e-12, 3.3363436e-12, 3.2001398e-12, 3.0753736e-12, 2.9412981e-12, 2.8318282e-12, 2.7257422e-12, 2.6288731e-12, 2.5280116e-12, 2.4349649e-12, 2.3522879e-12, 2.2679037e-12, 2.1983550e-12, 2.1299028e-12, 2.0707980e-12, 2.0218582e-12, 1.9736388e-12, 1.9333927e-12, 1.8971892e-12, 1.8695481e-12, 1.8457255e-12, 1.8302086e-12, 1.8152970e-12, 1.8151223e-12, 1.8158688e-12, 1.8262736e-12, 1.8429841e-12, 1.8642947e-12, 1.8901548e-12, 1.9243658e-12, 1.9637075e-12, 2.0105556e-12, 2.0624878e-12, 2.1202545e-12, 2.1877949e-12, 2.2587192e-12, 2.3343787e-12, 2.4192297e-12, 2.5103115e-12, 2.6072867e-12, 2.7062549e-12, 2.8174517e-12, 2.8174517e-12]),
        "ymax"     : np.array([7.4982907e-12, 7.1910016e-12, 6.9182003e-12, 6.6249200e-12, 6.3661005e-12, 6.1104359e-12, 5.8633804e-12, 5.6252128e-12, 5.3953612e-12, 5.1756951e-12, 4.9638200e-12, 4.7701171e-12, 4.5616013e-12, 4.3911866e-12, 4.2262299e-12, 4.0755193e-12, 3.9183857e-12, 3.7731478e-12, 3.6442761e-12, 3.5126965e-12, 3.4040659e-12, 3.2972015e-12, 3.2043850e-12, 3.1279352e-12, 3.0521448e-12, 2.9890066e-12, 2.9321201e-12, 2.8880447e-12, 2.8500918e-12, 2.8249788e-12, 2.8007072e-12, 2.7994768e-12, 2.7995189e-12, 2.8144977e-12, 2.8393388e-12, 2.8710337e-12, 2.9099934e-12, 2.9618605e-12, 3.0214992e-12, 3.0927929e-12, 3.1720739e-12, 3.2603093e-12, 3.3636599e-12, 3.4722126e-12, 3.5880952e-12, 3.7181683e-12, 3.8579312e-12, 4.0067802e-12, 4.1587643e-12, 4.3295963e-12, 4.3295963e-12]),
        "channels" : [
            (r"$\mathrm{u}\bar{\mathrm{u}} + \mathrm{c}\bar{\mathrm{c}}$", np.array([5.7928985e-12, 5.5545492e-12, 5.3436676e-12, 5.1170724e-12, 4.9162902e-12, 4.7204604e-12, 4.5290894e-12, 4.3461366e-12, 4.1703059e-12, 3.9987049e-12, 3.8352353e-12, 3.6854207e-12, 3.5246703e-12, 3.3934168e-12, 3.2662258e-12, 3.1498699e-12, 3.0290848e-12, 2.9170428e-12, 2.8178193e-12, 2.7171531e-12, 2.6338379e-12, 2.5506545e-12, 2.4799421e-12, 2.4205546e-12, 2.3628084e-12, 2.3149553e-12, 2.2710462e-12, 2.2372451e-12, 2.2081619e-12, 2.1896656e-12, 2.1712047e-12, 2.1709859e-12, 2.1712641e-12, 2.1831538e-12, 2.2025217e-12, 2.2284757e-12, 2.2590442e-12, 2.2993731e-12, 2.3460771e-12, 2.4020975e-12, 2.4638450e-12, 2.5323895e-12, 2.6127955e-12, 2.6978311e-12, 2.7877094e-12, 2.8894223e-12, 2.9977314e-12, 3.1142142e-12, 3.2318652e-12, 3.3649215e-12, 3.3649215e-12])),
            (r"$\mathrm{d}\bar{\mathrm{d}} + \mathrm{s}\bar{\mathrm{s}} + \mathrm{b}\bar{\mathrm{b}}$", np.array([1.8552176e-13, 1.7885086e-13, 1.7217283e-13, 1.6510577e-13, 1.5944207e-13, 1.5164441e-13, 1.4601514e-13, 1.3930258e-13, 1.3202519e-13, 1.2854728e-13, 1.2330167e-13, 1.1871509e-13, 1.1340433e-13, 1.0901027e-13, 1.0481896e-13, 1.0117454e-13, 9.6932317e-14, 9.3526940e-14, 9.0238533e-14, 8.6258414e-14, 8.3253847e-14, 8.1498874e-14, 7.8662857e-14, 7.7305575e-14, 7.5022648e-14, 7.2819518e-14, 7.1664335e-14, 7.0812969e-14, 7.0043876e-14, 6.8944971e-14, 6.8522087e-14, 6.8156305e-14, 6.8373568e-14, 6.8908752e-14, 6.9797327e-14, 6.9682621e-14, 7.0660333e-14, 7.2189766e-14, 7.3619284e-14, 7.5022326e-14, 7.7026006e-14, 7.9413264e-14, 8.2011611e-14, 8.4148408e-14, 8.7293834e-14, 8.9949315e-14, 9.3724817e-14, 9.6596463e-14, 1.0078114e-13, 1.0464156e-13, 1.0464156e-13]))
        ],
    },
]
metadata = {
    "arxiv": r"",
    "description": r"",
    "hepdata": r"",
    "initial_state_1": r"2212",
    "initial_state_2": r"2212",
    "lumi_id_types": r"pdg_mc_ids",
    "mg5amc_repo": r"http://bazaar.launchpad.net/~maddevelopers/mg5amcnlo/3.3.1/",
    "mg5amc_revno": r"981",
    "patches": r"",
    "pineappl_gitversion": r"v0.5.2-64-g402587e-dirty",
    "results_pdf": r"MSHT20nnlo_as118",
    "runcard_gitversion": r"9b81106-dirty",
    "tau_min": r"",
    "user_cuts": r"mmllmax=500.0",
    "x1_label": r"Mll",
    "x1_label_tex": r"$M_{\ell\bar{\ell}}$",
    "x1_unit": r"\giga\electronvolt",
    "x2_label": r"costh",
    "x2_label_tex": r"$\cos \theta^*$",
    "x2_unit": r"",
    "y_label": r"dsig/dcosth",
    "y_label_tex": r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\cos \theta^*}$",
    "y_unit": r"\pico\barn",
}
# end CLI data variables


if __name__ == "__main__":
    main(plot_panels)
"#;

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["plot", "--help"])
        .assert()
        .success()
        .stdout(
            HELP_STR.replace(
                "{}",
                &thread::available_parallelism()
                    .map_or(1, NonZeroUsize::get)
                    .to_string(),
            ),
        );
}

#[test]
fn default() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "plot",
            "--threads=1",
            "../test-data/LHCB_WP_7TEV_opt.pineappl.lz4",
            "NNPDF31_nlo_as_0118_luxqed",
            "NNPDF40_nnlo_as_01180=NNPDF4.0",
        ])
        .assert()
        .success()
        .stdout(DEFAULT_STR);
}

#[test]
fn drell_yan_afb() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "plot",
            "--asymmetry",
            "--threads=1",
            "../test-data/CMS_DY_14TEV_MLL_6000_COSTH.pineappl.lz4",
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(DRELL_YAN_AFB_STR);
}

#[test]
fn drell_yan_mass_slices() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args([
            "plot",
            "--no-conv-fun-unc",
            "--threads=1",
            "../test-data/NNPDF_DY_14TEV_BSM_AFB.pineappl.lz4",
            "NNPDF40_nnlo_as_01180",
        ])
        .assert()
        .success()
        .stdout(DRELL_YAN_MASS_SLICES_STR);
}
