import os
import numpy as np
import pandas as pd
import argparse
import warnings
import matplotlib
from kipoi_veff.external.concise.seqplotting_deps import add_letter_to_axis, VOCABS, letter_polygons
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

warnings.simplefilter("ignore")
warnings.filterwarnings("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

parser = argparse.ArgumentParser(description='Mutation Map Generator')
arguments_group = parser.add_argument_group("Parameters")
arguments_group.add_argument('-d', '--data_path', type=str, required=True,
                                help='data path')
arguments_group.add_argument('-o', '--output_path', type=str, required=True,
                                help='')
arguments_group.add_argument('-f', '--output_file', type=str, required=True,
                                help='')

args = vars(parser.parse_args())

data_path = args['data_path']
output_path = args['output_path']
output_filename = args['output_file']

mutation_maps = np.load(data_path)

def sequence_extractor(mutation_map):
    char_decoder = {
        0: "A",
        1: "C",
        2: "G",
        3: "U"
    }
    
    rna_sequence = ""
    rna_encoded = np.zeros_like(mutation_map)
    for position_idx in range(mutation_map.shape[0]):
        for i in range(4):
            if mutation_map[position_idx, i] == 0 and np.var(mutation_map[position_idx, :]) != 0:
                rna_encoded[position_idx, i] = 1
                rna_sequence += char_decoder[i]
                break
    return rna_encoded, rna_sequence

def center_cmap(cmap, vmax, vmin, center):
    # Centering of the colormap, taken from seaborn._HeatMapper implementation
    import matplotlib as mpl
    vrange = max(vmax - center, center - vmin)
    normlize = mpl.colors.Normalize(center - vrange, center + vrange)
    cmin, cmax = normlize([vmin, vmax])
    cc = np.linspace(cmin, cmax, 256)
    return mpl.colors.ListedColormap(cmap(cc))


def seqlogo_heatmap(letter_heights, heatmap_data, ovlp_var=None, vocab="DNA", ax=None, show_letter_scale=False,
                    cmap=None, cbar=True, cbar_kws=None, cbar_ax=None, limit_region=None, var_box_color="black",
                    show_var_id=True, box_alt=True, ref_seq=None):
    """
    Plot heatmap and seqlogo plot together in one axis.
    # Arguments
        letter_heights: "sequence length" x "vocabulary size" numpy array of seqlogo letter heights
        heatmap_data: "vocabulary size" x "sequence length" numpy array of heatmap values
    Can also contain negative values.
        vocab: str, Vocabulary name. Can be: DNA, RNA, AA, RNAStruct.
        ax: matplotlib axis
        box_alt: if true (default), variant box will be drawn on alternative sequence, otherwise on the reference sequence.
        ref_seq: str, reference sequence. If provided, will be provided as xticklabels.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.collections import PatchCollection
    if cmap is None:
        cmap = plt.cm.bwr

    seq_len = heatmap_data.shape[1]
    vocab_len = 4
    letter_rescaling = 4  # This way the heatmap and the letters are the same size on the plot

    # heatmap grid
    grid = np.mgrid[0.5:(seq_len + 0.5):1, -vocab_len:0:1].reshape(2, -1).T
    y_hm_tickpos = (np.arange(-vocab_len, 0, 1) + 0.5)[::-1]  # alphabet position with 0 on top
    y_seqlogo_tickpos = np.array([0, letter_rescaling])  # tuple of where the ticks for the seqlogo should be placed

    if ax is None:
        plt.figure(figsize=(20, 4))
        ax = plt.subplot(1, 1, 1)
    patches = []
    # add a circle
    for pos_tuple in grid:
        rect = mpatches.Rectangle(pos_tuple, 1.0, 1.0, ec="none")
        patches.append(rect)

    # Add colours to the heatmap - flip the alphabet order so that "A" is on top.
    colors = heatmap_data[::-1, :].T.reshape((seq_len * 4))
    # Centre the colours around 0
    cmap_centered = center_cmap(cmap, colors.max(), colors.min(), 0.0)
    collection = PatchCollection(patches, cmap=cmap_centered, alpha=1.0)
    collection.set_array(np.array(colors))

    # add the heatmap to the axis
    hm_ax_collection = ax.add_collection(collection)

    # rescale letters so that they look nice above the heatmap
    letter_heights_rescaled = np.copy(letter_heights)

    letter_height_scaling = (vocab_len / 2)

    if letter_heights_rescaled.min() < 0:
        lh_range = (letter_heights_rescaled.max() - letter_heights_rescaled.min()) / letter_height_scaling
        letter_y_offset = - letter_heights_rescaled.min() / lh_range
        letter_heights_rescaled = letter_heights_rescaled / lh_range
    else:
        letter_y_offset = 0.0
        letter_heights_rescaled /= letter_heights_rescaled.max() / letter_height_scaling

    assert letter_heights.shape[1] == len(VOCABS[vocab])
    x_range = [1, letter_heights.shape[0]]

    for x_pos, heights in enumerate(letter_heights_rescaled):
        letters_and_heights = sorted(zip(heights, list(VOCABS[vocab].keys())))
        y_pos_pos = letter_y_offset
        y_neg_pos = letter_y_offset
        for height, letter in letters_and_heights:
            color = VOCABS[vocab][letter]
            polygons = letter_polygons[letter]
            if height > 0:
                add_letter_to_axis(ax, polygons, color, 0.5 + x_pos, y_pos_pos, height)
                y_pos_pos += height
            else:
                add_letter_to_axis(ax, polygons, color, 0.5 + x_pos, y_neg_pos, height)
                y_neg_pos += height

    ax.set_xlim(x_range[0] - 1, x_range[1] + 1)
    ax.grid(False)
    ax.set_xticks(list(range(*x_range)) + [x_range[-1]])
    ax.set_aspect(aspect='auto', adjustable='box')

    # set the tick labels and make sure only the left axis is displayed
    if show_letter_scale:
        y_ticks = np.concatenate([y_hm_tickpos, y_seqlogo_tickpos])
        yticklabels = list(VOCABS[vocab].keys()) + ["%.2f" % letter_heights.min(), "%.2f" % letter_heights.max()]
        ax.spines['left'].set_bounds(y_seqlogo_tickpos[0], y_seqlogo_tickpos[1])
    else:
        y_ticks = y_hm_tickpos
        yticklabels = list(VOCABS[vocab].keys())
        ax.spines['left'].set_visible(False)

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(yticklabels)
    ax.axes.get_xaxis().set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.autoscale_view()

    if ref_seq is not None:
        xticklabels = list(ref_seq)
        ax.set_xticklabels(xticklabels)
        ax.axes.get_xaxis().set_visible(True)

    if ovlp_var is not None:
        # for every variant draw a rectangle
        for rel_pos, var_id, ref, alt in zip(ovlp_var["varpos_rel"], ovlp_var["id"], ovlp_var["ref"], ovlp_var["alt"]):
            # positions of ref and alt on the heatmap
            # This is for non-flipped alphabet
            # y_ref_lowlim = -vocab_len + list(VOCABS[vocab].keys()).index(ref[0])
            # y_alt_lowlim = -vocab_len + list(VOCABS[vocab].keys()).index(alt[0])
            # This is for the flipped alphabet
            y_ref_lowlim = list(VOCABS[vocab].keys()).index(ref[0]) * (-1) - 1
            y_alt_lowlim = list(VOCABS[vocab].keys()).index(alt[0][0]) * (-1) - 1
            # box drawing
            box_width = len(ref)
            # Deprecated: draw bax around ref and alt.
            # y_lowlim = min(y_ref_lowlim, y_alt_lowlim)
            # box_height = np.abs(y_ref_lowlim - y_alt_lowlim) + 1
            if box_alt:
                y_lowlim = y_alt_lowlim
            else:
                y_lowlim = y_ref_lowlim
            box_height = 1
            ax.add_patch(
                mpatches.Rectangle((rel_pos + 0.5, y_lowlim), box_width, box_height, fill=False, lw=2,
                                   ec=var_box_color))
            if show_var_id:
                # annotate the box
                ax.annotate(var_id, xy=(rel_pos + box_width + 0.5, y_lowlim + box_height / 2),
                            xytext=(rel_pos + box_width + 0.5 + 2, y_lowlim + box_height / 2),
                            arrowprops=dict(arrowstyle="->", connectionstyle="arc"),
                            bbox=dict(boxstyle="round,pad=.5", fc="0.9", alpha=0.7))

    if limit_region is not None:
        if not isinstance(limit_region, tuple):
            raise Exception("limit_region has to be tuple of (x_min, x_max)")
        ax.set_xlim(limit_region)

    if cbar:
        # Colorbar settings adapted from seaborn._HeatMapper implementation
        import matplotlib as mpl
        cbar_kws = {} if cbar_kws is None else cbar_kws
        cbar_kws.setdefault('ticks', mpl.ticker.MaxNLocator(6))
        cb = ax.figure.colorbar(hm_ax_collection, cbar_ax, ax, **cbar_kws)
        cb.outline.set_linewidth(0)
        # If rasterized is passed to pcolormesh, also rasterize the
        # colorbar to avoid white lines on the PDF rendering
        # if kws.get('rasterized', False):
        cb.solids.set_rasterized(True)
    return ax

pdf = PdfPages(os.path.join(output_path, f"{output_filename}.mutation.map.pdf"))
for i in range(mutation_maps.shape[0]):
    sequence_heatmap = mutation_maps[i]
    sequence_heatmap = np.concatenate([sequence_heatmap, np.reshape([0, -1, 0, 0], (4, 1))], axis=1)
    sequence_heights = np.absolute(sequence_heatmap).max(axis=0)
    sequence_encoded, sequence = sequence_extractor(sequence_heatmap.T)
    sequence_heights = np.multiply(sequence_encoded, (np.reshape(sequence_heights, (-1, 1))))
    plt.figure(figsize=(15, 10))
    seqlogo_heatmap(np.absolute(sequence_heights), sequence_heatmap, vocab='RNA')
    plt.savefig(pdf, format='pdf')
pdf.close()