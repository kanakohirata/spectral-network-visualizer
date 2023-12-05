def get_text_color(rgb):
    brightness = max(rgb) / 255
    if brightness > 0.65:
        color = '#000000'
    else:
        color = '#ffffff'

    return color
