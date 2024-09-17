def modes_within_value(mode1, mode2, value):
    return abs(mode1 - mode2) <= 0.2 * max(mode1, mode2)
