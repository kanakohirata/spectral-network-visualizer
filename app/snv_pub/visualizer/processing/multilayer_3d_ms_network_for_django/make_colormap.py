from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import plotly.colors as pcolors


def from_named_colorscale(name, N=256, gamma=1.0):
    colorscale = pcolors.get_colorscale(name)
    if not colorscale:
        return None

    colors = pcolors.colorscale_to_colors(colorscale)
    list_rgb_255 = []
    for color in colors:
        if color.startswith('rgb'):
            unlabeled_rgb = pcolors.unlabel_rgb(color)
        else:
            unlabeled_rgb = pcolors.hex_to_rgb(color)

        list_rgb_255.append(unlabeled_rgb)

    list_rgb_0to1 = [pcolors.unconvert_from_RGB_255(rgb) for rgb in list_rgb_255]

    linear_cmap = LinearSegmentedColormap.from_list('temp_cmap', colors=list_rgb_0to1, N=N, gamma=gamma)

    return linear_cmap


def from_colorscale(colorscale, N=256, gamma=1.0):
    colors = pcolors.colorscale_to_colors(colorscale)

    list_rgb_255 = []
    for color in colors:
        if color.startswith('rgb'):
            unlabeled_rgb = pcolors.unlabel_rgb(color)
        else:
            unlabeled_rgb = pcolors.hex_to_rgb(color)

        list_rgb_255.append(unlabeled_rgb)

    list_rgb_0to1 = [pcolors.unconvert_from_RGB_255(rgb) for rgb in list_rgb_255]

    linear_cmap = LinearSegmentedColormap.from_list('temp_cmap', colors=list_rgb_0to1, N=N, gamma=gamma)

    return linear_cmap


def scale_0to1(value, max_, min_):
    if isinstance(value, str):
        return value
    scaled_value = (value - min_) / (max_ - min_)
    return scaled_value


def split_colorscale(name, _n=255):
    linear_cmap = from_named_colorscale(name, _n)
    colorscale = []

    for i in range(linear_cmap.N):
        if i == 0:
            value = i
        else:
            value = i / (linear_cmap.N - 1)

        rgb = list(map(np.uint8, np.array(linear_cmap(value)[:3]) * 255))
        labeled_rgb = pcolors.label_rgb(rgb)
        colorscale.append([value, labeled_rgb])

    return colorscale
