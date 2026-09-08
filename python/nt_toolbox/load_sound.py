import numpy as np
import wave as wv


def load_sound(file, n0):

    x_raw = wv.open(file)
    n = x_raw.getnframes()
    x = np.frombuffer(x_raw.readframes(n), dtype="<i2").astype(float)
    x_raw.close()

    if file[::-1][:8][::-1] == "bird.wav":
        x = np.delete(
            x,
            list(range(6001))
            + list(range(12500, 15001))
            + list(range(22500, 24001))
            + list(range(32500, 34001)),
        )

    if n0 != 0 and n0 < n:
        x = x[:n0]

    return x / max(np.max(np.abs(x)), np.finfo(float).eps)
