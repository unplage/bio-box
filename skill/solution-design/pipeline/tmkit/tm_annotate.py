KD = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'D': -3.5, 'E': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9,
    'R': -4.5,
}


def kd_scan(seq, window=15, threshold=1.6):
    w = window
    n = len(seq)
    out = []
    for i in range(n):
        lo, hi = max(0, i - w // 2), min(n, i + w // 2 + 1)
        seg = seq[lo:hi]
        if not seg:
            out.append(float('-inf'))
            continue
        score = sum(KD.get(a, 0.0) for a in seg) / len(seg)
        out.append(score)
    segments = []
    i = 0
    while i < n:
        if out[i] > threshold:
            j = i
            while j < n and out[j] > threshold:
                j += 1
            segments.append((i + 1, j, max(out[i:j])))
            i = j
        else:
            i += 1
    return segments


def merge(segments, min_length=20):
    if not segments:
        return []
    merged = [list(segments[0])]
    for seg in segments[1:]:
        if seg[0] - merged[-1][1] <= 3:
            merged[-1][1] = seg[1]
            merged[-1][2] = max(merged[-1][2], seg[2])
        else:
            merged.append(list(seg))
    return [tuple(s) for s in merged if s[1] - s[0] + 1 >= min_length]


def detect_tm(seq, min_length=20, window=15, threshold=1.6):
    return merge(kd_scan(seq, window=window, threshold=threshold), min_length)


def try_deeptmhmm(seq):
    try:
        import deep_tmhmm
    except Exception:
        return None
    try:
        pred = deep_tmhmm.deeptmhmm.predict(seq)
        return pred
    except Exception as e:
        print(f'[warn] DeepTMHMM failed: {e}')
        return None


def tm_segments(seq, config, warn=True):
    if config.get('tm_segments'):
        return [(int(a), int(b)) for a, b in config['tm_segments']]
    deep = None
    if config.get('use_deeptmhmm', False):
        deep = try_deeptmhmm(seq)
        if deep is not None:
            segs = []
            for i, s in enumerate(deep['TM']):
                if s:
                    segs.append((i + 1, i + 1))
            merged = merge([[a, b, 9.9] for a, b, _ in segs], min_length=12)
            if merged:
                if warn:
                    print('[info] TM segments from DeepTMHMM:', merged)
                return [(a, b) for a, b, _ in merged]
    segs = detect_tm(seq)
    if warn:
        print('[warn] TM segments estimated by Kyte-Doolittle window (low accuracy). '
              'Provide `tm_segments` in config for reliable results.')
    return segs