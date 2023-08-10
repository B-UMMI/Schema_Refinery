# this file stores random functions that don't 
# necessarily belong anywhere else

def hex_to_rgb(hex_color: str) -> tuple:
    # translates an hexadecimal color code
    # into RGB color code
    hex_color = hex_color.lstrip("#")
    if len(hex_color) == 3:
        hex_color = hex_color * 2
    return int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)