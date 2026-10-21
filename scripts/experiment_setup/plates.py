"""384-well plate geometry."""
ROWS = "ABCDEFGHIJKLMNOP"
N_ROWS, N_COLS = 16, 24


def all_wells():
    return [f"{r}{c}" for r in ROWS for c in range(1, N_COLS + 1)]


def inner_wells():
    """The 308 non-edge wells (rows B-O, columns 2-23): edges evaporate, and over a
    10-day passage series that difference compounds."""
    return [f"{r}{c}" for r in ROWS[1:-1] for c in range(2, N_COLS)]


def split(well):
    return ROWS.index(well[0]), int(well[1:])


def neighbours(well):
    """The 8 surrounding wells."""
    r, c = split(well)
    out = set()
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if (dr or dc) and 0 <= r + dr < N_ROWS and 1 <= c + dc <= N_COLS:
                out.add(f"{ROWS[r + dr]}{c + dc}")
    return out
