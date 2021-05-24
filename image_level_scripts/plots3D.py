#!/usr/bin/python3
"""
3D attrs plots
"""

# Standard library imports
import os
import pandas as pd

# Third party imports
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns

sns.set(context='notebook', style='white',
        font='sans-serif', font_scale=2,
        rc={'figure.figsize': (15, 15)})


def plot_3D(c0, c1, xlabel, ylabel, zlabel, leglabel, outname, dt, loc, num=5):
    dt_back = dt.copy()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    dt = dt_back[dt_back['class'] == 1]
    xs = dt['0D']
    ys = dt['0E']
    zs = dt['1D']
    ft = 11.
    pnt3d = ax.scatter(xs, ys, zs, s=dt['1E']*ft,
                       c='darkorange',
                       zdir='z',
                       marker='o',
                       edgecolors='w',
                       alpha=0.8)

    dt = dt_back[dt_back['class'] == 0]
    xs = dt['0D']
    ys = dt['0E']
    zs = dt['1D']
    pnt3d1 = ax.scatter(xs, ys, zs, s=dt['1E']*ft,
                        c='lightseagreen',
                        zdir='z',
                        marker='D',
                        edgecolors='w',
                        alpha=0.8)
    ax.zaxis.set_major_locator(MaxNLocator(integer=True))

    ax.set_xlabel(xlabel, labelpad=15)
    ax.set_ylabel(ylabel, labelpad=15)
    ax.set_zlabel(zlabel, labelpad=15)
    legend1 = ax.legend([pnt3d, pnt3d1],
                        [c1, c0],
                        numpoints=1,
                        title="Diagnosis",
                        loc=[0.72, 0.25])
    ax.add_artist(legend1)

    # produce a legend with a cross section of sizes from the scatter
    handles, labels = pnt3d.legend_elements(prop="sizes", num=num, alpha=.6)
    legend2 = ax.legend(handles, labels, loc=loc, title=leglabel)

    for i in range(len(legend2.get_texts())):
        text = str(legend2.get_texts()[i])
        if "$" in text:
            size_number = text.split(
                ",")[-1].split("}$")[0].split("t{")[-1]
        else:
            size_number = text.split(
                ",")[-1].split("')")[0].split("'")[-1]

        text_show = float(size_number)/ft
        legend2.get_texts()[i].set_text(str(int(text_show)))

    ax.view_init(15, 125)

    plt.savefig(os.path.join("plots", outname), dpi=400, bbox_inches='tight')


def plot_3D_mciad(c0, c1, xlabel, ylabel, zlabel, leglabel, outname,
                  dt, loc, num=5):
    dt_back = dt.copy()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    dt = dt_back[dt_back['class'] == 1]
    xs = dt['0D']
    ys = dt['0E']
    zs = dt['1D']
    ft = 11.
    pnt3d = ax.scatter(xs, ys, zs, s=dt['1E']*ft,
                       c='lightseagreen',
                       zdir='z',
                       marker='D',
                       edgecolors='w',
                       alpha=0.8)

    dt = dt_back[dt_back['class'] == 0]
    xs = dt['0D']
    ys = dt['0E']
    zs = dt['1D']

    pnt3d1 = ax.scatter(xs, ys, zs, s=dt['1E']*ft,
                        c='darkorange',
                        zdir='z',
                        marker='o',
                        edgecolors='w',
                        alpha=0.8)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.zaxis.set_major_locator(MaxNLocator(integer=True))

    ax.set_xlabel(xlabel, labelpad=15)
    ax.set_ylabel(ylabel, labelpad=15)
    ax.set_zlabel(zlabel, labelpad=15)
    legend1 = ax.legend([pnt3d1, pnt3d],
                        [c0, c1],
                        numpoints=1,
                        title="Diagnosis",
                        loc=[0.72, 0.25])
    ax.add_artist(legend1)

    # produce a legend with a cross section of sizes from the scatter
    handles, labels = pnt3d.legend_elements(prop="sizes", num=num, alpha=.6)
    legend2 = ax.legend(handles, labels, loc=loc, title=leglabel)

    for i in range(len(legend2.get_texts())):
        text = str(legend2.get_texts()[i])
        if "$" in text:
            size_number = text.split(
                ",")[-1].split("}$")[0].split("t{")[-1]
        else:
            size_number = text.split(
                ",")[-1].split("')")[0].split("'")[-1]

        text_show = float(size_number)/ft
        legend2.get_texts()[i].set_text(str(int(text_show)))

    ax.view_init(15, 125)

    plt.savefig(os.path.join("plots", outname), dpi=400, bbox_inches='tight')


def get_data(fold_range, exp):
    input_folder = os.path.join('..', 'experiment_cmpb', "hippocampus",
                                "fold", "image_train_test_data")
    for fold in range(fold_range[0], fold_range[1]):
        name_test = os.path.join(input_folder.replace("fold", str(fold)),
                                 exp, 'fts_test.csv')
        df_test = pd.read_csv(name_test)
        cols = df_test.columns
        y = df_test.loc[:, 'class']
        ft = df_test.loc[:, cols.str.contains('tissues')]
        ft.columns = [col.split("_")[0] for col in ft.columns]
        ft = pd.concat([ft, y], axis=1)
        return ft


if __name__ == "__main__":

    ft = get_data([6, 7], 'cn_ad')
    plot_3D('CN', 'AD', '# CN left', '# CN right', '# AD right',
            "# AD left", 'cn_ad_points.png', ft, [0.12, 0.57], 5)

    ft = get_data([2, 3], 'cn_mci')
    plot_3D('CN', 'MCI', '# CN left', '# CN right', '# MCI right',
            "# MCI left", 'cn_mci_points.png', ft, [0.12, 0.57], 5)

    ft = get_data([5, 6], 'mci_ad')
    plot_3D_mciad('AD', 'MCI', '# MCI left', '# MCI right', '# AD right',
                  "# AD left", 'mci_ad_points.png', ft, [0.12, 0.57], 5)
