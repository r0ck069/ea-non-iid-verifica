#!/usr/bin/env python3
"""
EA non-IID — verifica sorgenti (versione Python) v1.1.0-beta1

Porta del motore usato dalla pagina HTML ea-non-iid-verifica-sorgenti-v1_1_0-beta1.html:
dieci stimatori SP 800-90B (Pista Non-IID), riprodotti con gli stessi calcoli di ea_non_iid
(repository usnistgov/SP800-90B_EntropyAssessment, commit 87c104d0ed4cbc96103e7b8b38d6f2c7e0a6b289).

Uso:
  python3 ea_verifica_sorgenti.py FILE [--fmt raw|bits|hex|dec] [--bps 1..8] [--offset N]
                                  [--strict] [--min-h 1.0] [--out-dir DIR]
  python3 ea_verifica_sorgenti.py --selftest

Esito: verde (UTILIZZABILE) con exit code 0, rosso (NON UTILIZZABILE) con exit code 1,
errore di lettura o di calcolo con exit code 2.

Lo stimatore misura la min-entropia dei dati caricati. Non certifica la sorgente.
Nessuna dipendenza esterna: solo la libreria standard di Python 3.
"""
import argparse
import base64
import math
import os
import re
import struct
import sys

VERSION = "EA non-IID verifica sorgenti (Python) v1.1.0-beta1"
REFERENCE = ("repository usnistgov/SP800-90B_EntropyAssessment, "
             "commit 87c104d0ed4cbc96103e7b8b38d6f2c7e0a6b289 (26/05/2026)")
TOL = 1e-6
WARN_N = 1000000   # campioni minimi richiesti per la stima (soglia NIST)
MIN_N = 1000       # sotto questa soglia il calcolo non è definito

Z = 2.5758293035489008
EPS = sys.float_info.epsilon
DENORM_MIN = 5e-324   # Number.MIN_VALUE in JavaScript
LN2 = math.log(2.0)
LN099 = math.log(0.99)
INF = float("inf")
NAN = float("nan")

# ---------------------------------------------------------------- utilità matematiche

def ln(x):
    if x > 0:
        return math.log(x)
    if x == 0:
        return -INF
    return NAN

def log2(x):
    return ln(x) / LN2

def finite_h(p):
    if not (p > 0):
        return INF
    if p >= 1:
        return 0.0
    return -log2(p)

def upper(p, n):
    return min(1.0, p + Z * math.sqrt(max(0.0, p * (1 - p) / (n - 1))))

def pow_(base, exp):
    try:
        return math.pow(base, exp)
    except (ValueError, OverflowError):
        return NAN

# ---------------------------------------------------------------- stimatori (SP 800-90B, 6.3)

def mcv(a, k):
    n = len(a)
    if n < 2:
        return None
    c = [0] * k
    for x in a:
        c[x] += 1
    mode = 0
    for i in range(1, k):
        if c[i] > c[mode]:
            mode = i
    p = c[mode] / n
    u = upper(p, n)
    return {"h": finite_h(u), "params": {"mode": mode, "count": c[mode], "p_hat": p, "p_u": u}}

def collision(a):
    n = len(a)
    if n < 2:
        return None
    i = 0
    v = 0
    s = 0
    while i < n - 1:
        if a[i] == a[i + 1]:
            t = 2
        elif i < n - 2:
            t = 3
        else:
            break
        v += 1
        s += t * t
        i += t
    if v < 2:
        return None
    X = i / v
    sig = math.sqrt(max(0.0, (s - i * X) / (v - 1)))
    xb = max(2.0, X - Z * sig / math.sqrt(v))
    p = (0.5 + math.sqrt(max(0.0, 1.25 - 0.5 * xb))) if xb < 2.5 else 0.5
    return {"h": finite_h(p), "params": {"v": v, "sum_t": i, "X_bar": X, "sigma_hat": sig,
                                         "X_bar_prime": xb, "p": p}}

def markov(a):
    n = len(a)
    if n < 2:
        return None
    C0 = 0
    C00 = 0
    C10 = 0
    for i in range(n - 1):
        if a[i] == 0:
            C0 += 1
            if a[i + 1] == 0:
                C00 += 1
        elif a[i + 1] == 0:
            C10 += 1
    C1 = n - 1 - C0
    P00 = C00 / C0 if C0 else 0.0
    P01 = (1 - P00) if C0 else 0.0
    P10 = C10 / C1 if C1 else 0.0
    P11 = (1 - P10) if C1 else 0.0
    if a[n - 1] == 0:
        C0 += 1
    P0 = C0 / n
    P1 = 1 - P0
    H = 128.0
    if P00 > 0:
        H = min(H, -log2(P0) - 127 * log2(P00))
    if P01 > 0 and P10 > 0:
        H = min(H, -log2(P0) - 64 * log2(P01) - 63 * log2(P10))
    if P01 > 0 and P11 > 0:
        H = min(H, -log2(P0) - log2(P01) - 126 * log2(P11))
    if P10 > 0 and P00 > 0:
        H = min(H, -log2(P1) - log2(P10) - 126 * log2(P00))
    if P10 > 0 and P01 > 0:
        H = min(H, -log2(P1) - 64 * log2(P10) - 63 * log2(P01))
    if P11 > 0:
        H = min(H, -log2(P1) - 127 * log2(P11))
    return {"h": min(H / 128, 1.0), "params": {"P_0": P0, "P_1": P1, "P_00": P00, "P_01": P01,
                                               "P_10": P10, "P_11": P11, "p_max": pow_(2, -H)}}

def G(z, d, n):
    Ai = 0.0
    Aic = 0.0
    first = 0.0
    firstc = 0.0
    B = 1 - z
    Bi = B
    for i in range(2, d + 1):
        term = log2(i) * Bi
        y = term - Aic
        t = Ai + y
        Aic = (t - Ai) - y
        Ai = t
        Bi *= B
    Ad1 = Ai
    truncated = False
    for i in range(d + 1, n):
        ai = log2(i) * Bi
        y = ai - Aic
        t = Ai + y
        Aic = (t - Ai) - y
        Ai = t + y * 0
        scaled = (n - i) * ai
        if scaled > 0:
            yy = scaled - firstc
            tt = first + yy
            firstc = (tt - first) - yy
            first = tt
        else:
            truncated = True
            break
        Bi *= B
    y = (n - d) * Ad1 - firstc
    t = first + y
    firstc = (t - first) - y
    first = t
    if not truncated:
        ai = log2(n) * Bi
        yy = ai - Aic
        tt = Ai + yy
        Aic = (tt - Ai) - yy
        Ai = tt
    return z * (z * first + (Ai - Ad1)) / (n - d)

def comExp(p, k, d, n):
    q = (1 - p) / (k - 1)
    return G(p, d, n) + (k - 1) * G(q, d, n)

def compression(a):
    b = 6
    d = 1000
    n = len(a) // b
    if n <= d:
        return None
    dic = [0] * 64
    X = 0.0
    Xc = 0.0
    S = 0.0
    Sc = 0.0
    for i in range(d):
        block = 0
        for j in range(b):
            block |= (a[i * b + j] & 1) << (b - j - 1)
        dic[block] = i + 1
    v = n - d
    for i in range(d, n):
        block = 0
        for j in range(b):
            block |= (a[i * b + j] & 1) << (b - j - 1)
        z = log2(i + 1 - dic[block])
        y = z - Xc
        t = X + y
        Xc = (t - X) - y
        X = t
        zz = z * z
        yy = zz - Sc
        tt = S + yy
        Sc = (tt - S) - yy
        S = tt
        dic[block] = i + 1
    X /= v
    sigma = 0.5907 * math.sqrt(max(0.0, S / (v - 1) - X * X))
    Xp = X - Z * sigma / math.sqrt(v)
    p = 1 / 64
    if comExp(p, 64, d, n) > Xp:
        lo = p
        hi = 1.0
        lv = INF
        hv = -INF
        mid = (lo + hi) / 2
        for _ in range(100):
            mid = (lo + hi) / 2
            val = comExp(mid, 64, d, n)
            if abs(val - Xp) <= abs(Xp) * EPS * 4 + DENORM_MIN:
                break
            if Xp < val:
                lo = mid
                lv = val
            else:
                hi = mid
                hv = val
            if lo >= hi or not (min(lv, hv) <= Xp <= max(lv, hv)):
                break
        p = mid
    else:
        p = 1 / 64
    h = (-log2(p) / 6) if p > 1 / 64 else 1.0
    return {"h": h, "params": {"X_bar": X, "sigma_hat": sigma, "X_bar_prime": Xp, "p": p}}

def build_sa(a):
    """Array dei suffissi (unico per un dato testo): raddoppio dei prefissi."""
    n = len(a)
    if n == 0:
        return []
    rank = list(a)
    sa = sorted(range(n), key=lambda i: rank[i])
    k = 1
    while True:
        keys = [(rank[i], rank[i + k] if i + k < n else -1) for i in range(n)]
        sa.sort(key=lambda i: keys[i])
        new = [0] * n
        for idx in range(1, n):
            new[sa[idx]] = new[sa[idx - 1]] + (1 if keys[sa[idx - 1]] != keys[sa[idx]] else 0)
        rank = new
        if rank[sa[-1]] == n - 1:
            break
        k <<= 1
        if k >= n:
            break
    return sa

def kasai(a, sa):
    n = len(a)
    rank = [0] * n
    lcp = [0] * n
    for i in range(n):
        rank[sa[i]] = i
    h = 0
    for i in range(n):
        ri = rank[i]
        if ri == 0:
            continue
        j = sa[ri - 1]
        while i + h < n and j + h < n and a[i + h] == a[j + h]:
            h += 1
        lcp[ri] = h
        if h > 0:
            h -= 1
    return lcp

def q_and_s(lcp, n):
    Q = [0.0] * (n + 1)
    S = [0.0] * (n + 1)
    parent = list(range(n))
    sz = [1] * n
    head = [-1] * (n + 1)
    nxt = [-1] * n
    for i in range(1, n):
        x = lcp[i]
        nxt[i] = head[x]
        head[x] = i
    cur_max = 1
    for t in range(n - 1, 0, -1):
        S[t] = S[t + 1] if t < n - 1 else 0.0
        e = head[t]
        while e != -1:
            a_ = e - 1
            b_ = e
            ra = a_
            rb = b_
            while parent[ra] != ra:
                ra = parent[ra]
            while parent[rb] != rb:
                rb = parent[rb]
            if ra != rb:
                if sz[ra] < sz[rb]:
                    ra, rb = rb, ra
                S[t] += sz[ra] * sz[rb]
                parent[rb] = ra
                sz[ra] += sz[rb]
                if sz[ra] > cur_max:
                    cur_max = sz[ra]
            e = nxt[e]
        Q[t] = cur_max
    return Q, S

def tuple_lrs(a):
    n = len(a)
    if n < 2:
        return {"tuple": None, "lrs": None}
    sa = build_sa(a)
    lcp = kasai(a, sa)
    Q, S = q_and_s(lcp, n)
    u = 1
    while u <= max(1, n - 1) and Q[u] >= 35:
        u += 1
    v = 0
    for i in range(1, len(Q)):
        if Q[i] >= 2:
            v = i
    if u > v:
        return {"tuple": None, "lrs": None}
    Pmax = 0.0
    for i in range(1, u):
        cp = Q[i] / (n - i + 1)
        pm = pow_(cp, 1 / i)
        if pm > Pmax:
            Pmax = pm
    pu = upper(Pmax, n)
    tup = {"h": finite_h(pu), "params": {"t": u - 1, "p_hat_max": Pmax, "p_u": pu}}
    Pmax = 0.0
    for i in range(u, v + 1):
        den = (n - i + 1) * (n - i) / 2
        cp = S[i] / den
        if cp > 0:
            pm = pow_(cp, 1 / i)
            if pm > Pmax:
                Pmax = pm
    pu = upper(Pmax, n)
    lrs = {"h": finite_h(pu), "params": {"u": u, "v": v, "p_hat": Pmax, "p_u": pu}}
    return {"tuple": tup, "lrs": lrs}

def prediction_estimate(C, N, maxrun, k):
    if not (N > 1) or not (k >= 2):
        return None
    pg = C / N
    pgp = min(1.0, pg + Z * math.sqrt(pg * (1 - pg) / (N - 1))) if pg > 0 else 1 - pow_(0.01, 1 / N)
    cur = max(1 / k, pgp)
    pl = None

    def pef(p):
        q = 1 - p
        x = 1.0
        last = 0.0
        i = 0
        while i <= 65 and (x - last) > EPS * x:
            last = x
            x = 1 + q * pow_(p, maxrun + 1) * pow_(x, maxrun + 2)
            if not (x >= last and p * x <= 1):
                break
            i += 1
        return ln(1 - p * x) - ln((maxrun + 2 - (maxrun + 1) * x) * q) - (N + 1) * ln(x)

    if cur < 1 and pef(cur) > LN099:
        lo = cur
        hi = 1.0
        p = (lo + hi) / 2
        lv = INF
        hv = -INF
        for _ in range(1076):
            if pef(p) > LN099:
                lo = p
                lv = pef(p)
            else:
                hi = p
                hv = pef(p)
            np_ = (lo + hi) / 2
            if np_ == p:
                break
            p = np_
            if not (min(lv, hv) <= LN099 <= max(lv, hv)):
                break
        pl = p
        cur = max(cur, pl)
    return {"h": finite_h(cur), "params": {"N": N, "P_global": pg, "P_global_prime": pgp,
                                           "P_local": pl, "C": C, "r": maxrun + 1}}

def multi_mcw(a, k):
    W = [63, 255, 1023, 4095]
    L = len(a)
    if L < W[3] + 1:
        return None
    score = [0] * 4
    cnt = [[0] * k for _ in range(4)]
    pos = [[0] * k for _ in range(4)]
    freq = [0] * 4
    maxc = [0] * 4
    for j in range(4):
        for i in range(W[j]):
            x = a[i]
            cnt[j][x] += 1
            pos[j][x] = i
            if maxc[j] <= cnt[j][x]:
                maxc[j] = cnt[j][x]
                freq[j] = x
    winner = 0
    C = 0
    run = 0
    maxrun = 0
    N = L - W[0]
    for i in range(W[0], L):
        if freq[winner] == a[i]:
            C += 1
            run += 1
            if run > maxrun:
                maxrun = run
        else:
            run = 0
        for j in range(4):
            if freq[j] == a[i] and i >= W[j]:
                score[j] += 1
                if score[j] >= score[winner]:
                    winner = j
        for j in range(4):
            if i >= W[j]:
                old = a[i - W[j]]
                x = a[i]
                cnt[j][old] -= 1
                cnt[j][x] += 1
                pos[j][x] = i
                if old != freq[j] and maxc[j] <= cnt[j][x]:
                    maxc[j] = cnt[j][x]
                    freq[j] = x
                elif old == freq[j]:
                    maxc[j] -= 1
                    mp = i - W[j]
                    for z in range(k):
                        if maxc[j] < cnt[j][z] or (maxc[j] == cnt[j][z] and mp <= pos[j][z]):
                            maxc[j] = cnt[j][z]
                            freq[j] = z
                            mp = pos[j][z]
    return prediction_estimate(C, N, maxrun, k)

def lag(a, k):
    L = len(a)
    if L <= 2:
        return None
    D = 128
    score = [0] * D
    rings = [[] for _ in range(k)]
    winner = 0
    high = 0
    C = 0
    run = 0
    maxrun = 0
    rings[a[0]].append(0)
    for i in range(1, L):
        if i - winner - 1 >= 0 and a[i] == a[i - winner - 1]:
            C += 1
            run += 1
            if run > maxrun:
                maxrun = run
        else:
            run = 0
        rb = rings[a[i]]
        cut = i - D if i >= D else 0
        for z in range(len(rb) - 1, -1, -1):
            if rb[z] < cut:
                del rb[0:z + 1]
                break
            off = i - rb[z] - 1
            score[off] += 1
            sc = score[off]
            if sc >= high:
                winner = off
                high = sc
        if len(rb) >= D:
            rb.pop(0)
        rb.append(i)
    return prediction_estimate(C, L - 1, maxrun, k)

def best_post(v):
    b = 0
    for i in range(1, len(v)):
        if v[i] > v[b] or (v[i] == v[b] and i > b):
            b = i
    return b, v[b]

def multi_mmc(a, k):
    D = 16
    MAX = 100000
    L = len(a)
    if L < 4:
        return None
    maps = [dict() for _ in range(D)]
    entries = [0] * D
    score = [0] * D
    winner = 0
    C = 0
    run = 0
    maxrun = 0
    for d in range(D):
        if not (d < L - 1):
            break
        key = bytes(a[0:d + 1])
        v = maps[d].get(key)
        if v is None:
            v = [0] * k
            maps[d][key] = v
        v[a[d + 1]] = 1
        entries[d] = 1
    for i in range(2, L):
        found = False
        cur_winner = winner
        for d in range(D):
            if not (d <= i - 2):
                break
            key = bytes(a[i - d - 1:i])
            m = maps[d]
            v = None
            if d == 0 or found:
                v = m.get(key)
                found = v is not None and best_post(v)[1] > 0
            else:
                found = False
            if found:
                sym, _count = best_post(v)
                if sym == a[i]:
                    score[d] += 1
                    if score[d] >= score[winner]:
                        winner = d
                    if d == cur_winner:
                        C += 1
                        run += 1
                        if run > maxrun:
                            maxrun = run
                elif d == cur_winner:
                    run = 0
                if v[a[i]] != 0:
                    v[a[i]] += 1
                elif entries[d] < MAX:
                    v[a[i]] = 1
                    entries[d] += 1
            elif entries[d] < MAX:
                v = [0] * k
                v[a[i]] = 1
                m[key] = v
                entries[d] += 1
    return prediction_estimate(C, L - 2, maxrun, k)

def lz78y(a, k):
    B = 16
    MAX = 65536
    L = len(a)
    if L <= B + 2:
        return None
    D = [dict() for _ in range(B)]
    size = 0
    for j in range(1, B + 1):
        key = bytes(a[B - j:B])
        v = D[j - 1].get(key)
        if v is None:
            v = [0] * k
            D[j - 1][key] = v
        v[a[B]] = 1
        size += 1
    C = 0
    run = 0
    maxrun = 0
    N = L - B - 1
    for i in range(B + 1, L):
        have = False
        pred = 0
        max_count = 0
        for j in range(B, 0, -1):
            key = bytes(a[i - j:i])
            m = D[j - 1]
            v = m.get(key)
            if v is not None:
                sym, count = best_post(v)
                if count > max_count:
                    max_count = count
                    pred = sym
                    have = True
                v[a[i]] += 1
            elif size < MAX:
                v = [0] * k
                v[a[i]] = 1
                m[key] = v
                size += 1
        if have and pred == a[i]:
            C += 1
            run += 1
            if run > maxrun:
                maxrun = run
        else:
            run = 0
    return prediction_estimate(C, N, maxrun, k)

# ---------------------------------------------------------------- analisi completa

NAMES_SYMBOL = ["MCV", "Collisione", "Markov", "Compressione", "t-Tuple", "LRS",
                "MultiMCW", "Lag", "MultiMMC", "LZ78Y"]

def analyze(symbols, bps):
    k = 1 << bps
    vals = list(symbols)
    obs = set(vals)
    if len(obs) == 1:
        return {"rows": [], "H_original": 0, "H_bitstring": None if bps == 1 else 0, "final": 0,
                "bitRows": [], "error": "Alfabeto di un solo simbolo: nessuna entropia."}
    mp = {}
    nxt = 0
    for v in range(k):
        if v in obs:
            mp[v] = nxt
            nxt += 1
    mapped = [mp[v] for v in vals]
    bits = [(v >> (bps - 1 - j)) & 1 for v in vals for j in range(bps)]

    def calc(a, kk, bit):
        rows = []

        def add(name, x):
            rows.append({"name": name,
                         "status": "calcolato" if x else "non eseguito: campioni insufficienti",
                         "h": x["h"] if x else None, "params": (x or {}).get("params", {})})

        def na(name):
            rows.append({"name": name,
                         "status": "n/a: stimatore definito dal NIST solo per binario",
                         "h": None, "params": {}})

        add("MCV", mcv(a, kk))
        if bit or kk == 2:
            add("Collisione", collision(a))
        else:
            na("Collisione")
        if bit or kk == 2:
            add("Markov", markov(a))
        else:
            na("Markov")
        if bit or kk == 2:
            add("Compressione", compression(a))
        else:
            na("Compressione")
        tl = tuple_lrs(a)
        add("t-Tuple", tl["tuple"])
        add("LRS", tl["lrs"])
        add("MultiMCW", multi_mcw(a, kk))
        add("Lag", lag(a, kk))
        add("MultiMMC", multi_mmc(a, kk))
        add("LZ78Y", lz78y(a, kk))
        hs = [r["h"] for r in rows if r["h"] is not None]
        return {"rows": rows, "H": min(hs) if hs else None}

    o = calc(mapped, k, False)
    b = calc(bits, 2, True)
    h_orig = min(bps, bps if o["H"] is None else o["H"])
    h_bit = min(1, 1 if b["H"] is None else b["H"])
    final = min(h_orig, bps * h_bit)
    return {"rows": o["rows"], "H_original": h_orig,
            "H_bitstring": None if bps == 1 else h_bit, "final": final,
            "bitRows": b["rows"], "samples": len(symbols), "alphabet": len(obs),
            "bitsLength": len(bits), "bps": bps}

# ---------------------------------------------------------------- lettura del file grezzo

def read_wav(buf):
    if len(buf) < 12 or buf[0:4] != b"RIFF" or buf[8:12] != b"WAVE":
        return None
    p = 12
    fmt = None
    data = None
    while p + 8 <= len(buf):
        cid = buf[p:p + 4]
        sz = struct.unpack_from("<I", buf, p + 4)[0]
        body = p + 8
        if cid == b"fmt " and body + 16 <= len(buf):
            code, ch, _rate, _br, _align, bits = struct.unpack_from("<HHIIHH", buf, body)
            sub = struct.unpack_from("<H", buf, body + 24)[0] if (code == 0xFFFE and body + 26 <= len(buf)) else 0
            fmt = {"format": code, "channels": ch, "bits": bits, "subformat": sub}
        elif cid == b"data":
            data = (body, min(len(buf), body + sz))
        p = body + sz + (sz & 1)
    if not fmt or not data:
        return None
    fmt["dataStart"], fmt["dataEnd"] = data
    return fmt

def extract(buf, o):
    mask = (1 << o["bps"]) - 1
    notes = []
    vals = []
    used = 0
    skipped = 0
    source = ""
    wav = read_wav(buf)
    if o["fmt"] == "raw":
        if wav:
            pcm = wav["format"] == 1 or (wav["format"] == 0xFFFE and wav["subformat"] == 1)
            flt = wav["format"] == 3 or (wav["format"] == 0xFFFE and wav["subformat"] == 3)
            if not pcm and not flt:
                raise ValueError("WAV formato %d non supportato: ammessi PCM e virgola mobile" % wav["format"])
            if flt and wav["bits"] not in (32, 64):
                raise ValueError("WAV float a %d bit non supportato (ammessi 32, 64)" % wav["bits"])
            if not flt and (wav["bits"] < 8 or wav["bits"] > 32 or wav["bits"] % 8):
                raise ValueError("WAV PCM a %d bit non supportato (ammessi 8, 16, 24, 32)" % wav["bits"])
            source = ("WAV float " if flt else "WAV PCM ") + "%d bit, %d canali, blocco data" % (wav["bits"], wav["channels"])
            if wav["channels"] > 1:
                notes.append("Canali interleaved: i campioni dei canali sono concatenati, non separati.")
            bp = wav["bits"] // 8
            for i in range(wav["dataStart"], wav["dataEnd"] - bp + 1, bp):
                vals.append(buf[i] & mask)
            used = len(vals)
            if o["offset"]:
                notes.append('Il campo "Salta" non si applica al WAV: si usa il blocco data.')
        else:
            source = "byte grezzi"
            if o["offset"] >= len(buf):
                raise ValueError("Offset oltre la fine del file")
            vals = [b & mask for b in buf[o["offset"]:]]
            used = len(vals)
    else:
        text = buf.decode("utf-8", errors="replace")
        lines = re.split(r"\r?\n", text)[o["offset"]:]
        source = {"bits": "testo: bit 0/1", "hex": "testo: byte esadecimali",
                  "dec": "testo: byte decimali"}[o["fmt"]]
        bit_parts = []
        for idx, ln in enumerate(lines):
            t = ln.strip()
            if t == "":
                continue
            ok = False
            if o["fmt"] == "bits":
                s = re.sub(r"\s+", "", t)
                ok = re.fullmatch(r"[01]+", s) is not None
                if ok:
                    bit_parts.append(s)
            else:
                toks = [x for x in re.split(r"[\s,;:_-]+", t) if x]
                if o["fmt"] == "hex":
                    ok = len(toks) > 0 and all(re.fullmatch(r"[0-9a-fA-F]{2}", x) for x in toks)
                else:
                    ok = len(toks) > 0 and all(re.fullmatch(r"[0-9]{1,3}", x) and int(x) <= 255 for x in toks)
                if ok:
                    for x in toks:
                        vals.append((int(x, 16) if o["fmt"] == "hex" else int(x)) & mask)
            if ok:
                used += 1
            elif o["skipBad"]:
                skipped += 1
            else:
                raise ValueError('Riga %d non conforme: "%s"' % (idx + 1 + o["offset"], t[:60]))
        if o["fmt"] == "bits":
            bs = "".join(bit_parts)
            n = len(bs) // o["bps"]
            rem = len(bs) - n * o["bps"]
            if rem:
                notes.append("%d bit finali scartati per completare l'ultimo simbolo." % rem)
            for i in range(n):
                vals.append(int(bs[i * o["bps"]:(i + 1) * o["bps"]], 2))
    return {"sym": vals, "notes": notes, "used": used, "skipped": skipped, "source": source}

# ---------------------------------------------------------------- autotest contro i riferimenti NIST

EMBED_DATA = {
"rand1_short": "AQEBAAABAAEAAAABAQAAAAEAAQABAAEBAAEAAAABAAEAAAEBAAEBAAAAAQAAAAEBAAEBAQABAQEAAAEAAQAAAAABAQAAAQEBAAEBAQEAAAABAAEBAQEBAAAAAAABAAEAAAAAAAEAAAABAQABAAEAAQAAAAEAAQEBAQEBAAABAQABAQAAAQAAAAEAAQABAQABAQEAAQEAAQAAAQAAAQEBAQABAQAAAAEBAQEBAAABAAEBAQEAAAABAQAAAQABAAEBAQEBAQEBAAABAQEAAQEAAQEAAAAAAAEBAAABAAEBAQAAAQEBAQEAAAAAAQEAAAEBAQAAAAEBAQEAAAEBAAEAAAEAAQEBAAABAAEBAAEAAAABAAEBAAABAAEAAAAAAQEAAAEBAQEAAQABAAEAAAEAAAEAAQEBAAEBAQABAAEAAAEAAQABAQAAAAEAAAAAAAEAAQAAAQEBAAEAAAABAAABAAEBAAEAAAEBAQEBAAEAAQABAAAAAAEAAAABAQAAAQEAAAEBAAEAAQEAAQEAAQEBAQAAAQAAAQAAAAABAAABAQEBAQABAQEBAAEBAAABAAEAAAAAAQEAAAABAAEBAAABAQEAAQEBAAEBAQEAAAAAAAEAAAEBAQEAAQAAAAABAQABAAEBAAABAAEAAQEBAQAAAQAAAQEAAQABAQEBAQABAAEAAAAAAAEAAQEBAAEAAQABAQABAAEBAQEAAQAAAQAAAAEAAQABAAABAAEAAQEBAQAAAQEAAQAAAAEBAQABAAEAAAEBAAAAAAAAAQEBAQABAQABAQAAAAEAAAEAAQAAAQAAAAAAAAABAAEBAQABAAEAAAABAQABAAEAAQABAQEAAAAAAAAAAAABAAEBAQAAAQABAQEBAQEAAAEBAQAAAAABAAEBAAEAAQABAQABAQEBAQEBAQEAAQEBAAEAAAEBAAEAAAAAAAAAAAABAAAAAQABAQAAAQEBAAEAAQABAAABAQAAAAABAQAAAAEAAQEBAAABAQEBAQEAAQAAAAAAAAEBAQAAAAEBAAABAAEBAQABAAAAAQAAAAEBAQAAAAEBAAABAAAAAQEBAAEBAQABAAEBAAAAAAEAAQEAAQAAAQABAAEBAAEAAAAAAAEAAAAAAQABAAAAAAEBAQEBAAABAAAAAQEAAQAAAAEBAQEAAQABAQABAAEAAQEAAAEAAAEAAQEBAQEBAAAAAQAAAAABAAEAAAABAAABAAEBAQAAAQEAAQABAAABAQEAAQAAAAEBAQABAAEAAAABAAEAAQABAQEBAQEBAQEBAAEBAAAAAQEBAQABAAAAAAAAAQABAQABAQEAAAABAAABAAEBAQEBAAABAAAAAQEBAQEBAAABAAABAQEBAQEAAQAAAQAAAAABAQEAAQEAAAABAAAAAQEAAQABAQEAAQABAQAAAAAAAQEAAQEAAAABAQEBAQEAAQEAAAAAAAABAAABAQAAAAEBAAAAAAEAAAABAQEBAQABAQEBAQEAAQEBAQABAAABAAAAAAEAAAABAQABAAEAAQABAQABAQEBAAEAAAAAAAEBAAAAAQABAAAAAQEAAAEAAQAAAAABAQEAAAEAAQEBAQABAQAAAAAAAQEAAAABAAEBAQABAAEAAQABAAABAQEBAAEAAAABAQABAAABAQABAQEBAQEBAQAAAQABAQAAAAAAAQEBAAEBAAAAAAAAAAEBAQEAAQABAQEBAQEBAAEBAQEBAAAAAQEBAAABAQEAAQEAAAAAAQABAAAAAAEAAQABAQABAAABAAABAQABAAEAAAEAAQAAAAABAAEBAQAAAQAAAAAAAAAAAQABAAEAAQEAAQABAQEBAAABAAEAAAEBAQABAAABAAEBAAAAAQABAAEAAQEAAQEAAQEAAQABAQEBAQABAQABAQABAAEAAQABAAEAAQABAAEAAQEAAAEAAAEBAAAAAAEBAQEBAAAAAQABAAAAAAEBAAEAAQAAAAABAQEAAAEBAAABAQAAAAEBAQEAAQAAAQEAAQABAAEBAQABAQEAAQAAAQEBAAABAAAAAAAAAQABAQEAAAAAAQEAAQEBAQAAAQAAAAEBAQEBAQABAAEAAAABAQEAAQEBAQEBAQAAAAEAAQEAAAABAAAAAQEBAAAAAAAAAAAAAQEBAQAAAAEBAAABAAEBAQEBAQEAAAEAAQEAAQEAAAEBAQEBAAEBAQAAAAEBAQABAQEBAAEBAQAAAQEAAQABAAAAAAEAAAAAAQEBAQAAAQAAAQAAAAEBAAEAAAEBAAEBAQAAAAABAQEAAAABAAAAAQAAAAAAAAEAAQABAAEBAAEAAAABAQEAAAEBAQEBAQEAAQEBAAEBAAABAQEBAAEBAAEAAAAAAAABAAABAQAAAQEAAQEBAAEAAQABAQEBAAABAAAAAQAAAQABAAEAAQABAQEBAQAAAQEAAAEBAAABAAAAAAABAQEAAQAAAQEAAQEBAQEBAAEAAAEAAAABAQEBAQEAAAEAAAAAAAAAAQABAQEAAQEAAAABAAEAAAEAAAEBAQABAQABAQABAAEAAAAAAAEAAQABAAABAQABAAEAAAABAQABAAAAAQABAQEBAQAAAAAAAAAAAAEBAQABAQABAAEBAQEBAAEBAAEBAAABAQEAAQEAAAEBAAAAAAAAAQEAAQEBAAEBAQAAAQEAAQEAAAAAAQABAQEAAAABAQABAAEBAAAAAAABAAAAAQEAAAEBAQEBAAEAAQAAAAEAAQEAAQAAAQEAAAAAAQAAAQAAAQAAAQEBAAABAAEAAAEAAQABAQEBAAEAAAEBAAABAAEAAAEBAQEAAQAAAAAAAQAAAQEAAQABAQEBAAAAAQAAAAEBAQEAAAEAAAAAAQEAAQEAAAEAAAABAAABAQAAAQEBAQEAAAEBAQEAAAAAAQAAAAAAAAEAAQAAAQEAAAEBAAAAAQAAAAEBAQAAAAABAAEBAAEAAAEBAAEBAAAAAAAAAQABAAEBAAEAAAEBAAABAQABAQEBAAEBAQAAAAEAAQEAAAABAQAAAQABAQEBAAAAAQEBAQEBAAEAAAAAAAEBAQABAQABAQAAAQABAAAAAAABAQEBAAEBAQEAAAABAQEAAQABAQEBAQEAAQAAAQABAAEAAQEBAAEBAQAAAQABAAEAAQEAAQABAAABAAAAAQEAAQEAAQEAAAEBAAEBAAEBAAABAAEAAAABAAEBAQAAAAAAAAEAAQABAQEAAQAAAAEAAQEBAAABAQEAAQAAAQEBAAABAAABAQAAAAEAAQABAQABAAEBAAAAAQAAAQEBAQAAAQEBAQEAAAEAAAABAQEBAQEAAAEBAQABAAEBAAAAAQABAAEBAQABAQEBAQEBAQEAAQABAAEAAAEBAAEBAAAAAAAAAQEBAAEAAQAAAQEBAQEBAQEBAQAAAAEAAQABAQEAAAEAAAEAAQEAAAEAAAEBAAABAQAAAAEAAAEAAAEAAQAAAAAAAQABAQEAAQEBAAAAAQAAAAABAAABAQAAAAEAAAABAQAAAAEBAAABAQEAAQABAQABAAEAAAEAAQABAQAAAAABAAABAQABAAAAAQAAAQAAAQABAAEBAAABAQEBAQAAAQAAAAEBAAABAAEBAQEAAQEBAAAAAQEBAQEBAQEAAQEBAAEBAQEBAAAAAQEBAQAAAAEAAAEAAAABAQABAAEAAQABAQABAQABAAAAAQABAQAAAAEBAQEAAAEBAQEBAAABAAEAAAABAQAAAQAAAAABAAEBAQABAAEBAQAAAAABAQEBAAEBAQAAAAABAAEAAQABAAEAAAABAQAAAQEAAAABAQABAAABAQAAAAAAAAABAAAAAAEAAAEBAQEBAQABAQEBAQAAAAABAQAAAAABAAABAAAAAQEBAAEAAQAAAAAAAQAAAAEAAQEBAAABAQABAAEBAAAAAAEAAQABAAABAQABAAEAAQABAQEBAAEBAAABAAEBAAABAAABAAEAAAAAAQABAQEBAAEBAQAAAQABAQAAAQAAAQEBAAAAAAEAAAEAAQEBAAABAQABAQEAAQAAAQEBAQAAAQAAAAEBAAABAAABAQEAAQABAQEAAAEBAQABAQEBAQABAQEBAAEBAAAAAAABAAEBAQEAAQABAQEAAQAAAQABAAEBAQAAAAEBAQEAAQAAAQEAAAEBAAABAAEBAQEBAQAAAQAAAQEAAAAAAQEAAQAAAAEAAQEAAAAAAQEBAQEAAQEBAQABAAAAAAAAAAABAAABAQAAAQEBAAEBAAEBAAABAQEBAQAAAAEBAAEAAQAAAAABAQEBAAAAAQABAAEBAQEBAQEAAQAAAQAAAQABAQEBAAEAAQAAAAABAQEBAAEBAQAAAAAAAAAAAQAAAQABAAEAAAABAQEAAAEBAQEAAAEBAAEAAQEBAAAAAAABAQAAAAEBAQEBAQEBAQAAAQEBAQABAAABAQEBAQEAAQEBAAEAAQEBAAEAAQEBAQEAAQAAAAEAAAABAQEBAQEAAQEAAQEAAQAAAQABAQABAQEAAQABAQABAAEAAQABAAABAAAAAQEAAAABAAEAAQABAQABAQABAAABAQEBAAAAAAEBAQEBAAABAAAAAAAAAQAAAQAAAQEAAQEAAQABAQAAAAABAAABAQEAAQABAQEBAAABAAEBAQABAQEAAQABAQEAAAEBAAEBAAAAAAAAAQEBAQABAAAAAAAAAAEAAQABAQEAAAABAQAAAAAAAAAAAAEBAAABAAEAAAEAAAAAAQEAAQEBAQABAAEAAQABAQEBAAEBAAAAAAABAQABAQEAAAEAAAAAAAEAAQEAAQABAAAAAQEAAQABAQABAQEAAAABAQEAAQABAAABAQABAAEBAQEBAAABAAEAAQAAAAAAAQABAQABAAABAAABAAEBAAAAAAEBAAEBAQAAAQAAAAAAAQEAAQEAAAEBAAAAAAEBAQABAQAAAAEAAQEBAAEBAQABAAEBAAABAAEBAQEAAQEAAQEAAQEAAQEBAQEAAQEAAAABAQAAAAABAAAAAQABAAABAQEBAQEBAQAAAQAAAAEAAQAAAAEBAQEAAQABAAEAAQEBAAABAQEBAAAAAAAAAAABAQAAAAEAAQAAAAEAAAABAQAAAQAAAQEBAQABAAEBAAEAAAABAAAAAQEAAQEBAQEBAQABAQEAAAEBAAEAAAEBAAEBAQEBAAAAAAABAAEAAQEBAAEBAQEAAQEAAAEBAAABAAEBAQEBAAAAAAEBAAEAAAEBAQABAAABAQEAAQEBAQABAAEBAAABAQEAAQABAAEBAQEBAQEBAQEAAQEBAAEAAQAAAQEAAQEAAAABAAABAAABAAABAQABAAEAAQEAAQAAAQAAAQABAQEAAAEBAQAAAQEBAQEBAAEBAAABAAAAAQEBAAEAAAEAAQABAAEBAQAAAAABAAAAAAEBAQEAAQABAAEBAAAAAAEAAAEAAAABAQEAAAABAAEAAAEBAAAAAQABAQEBAQABAAABAAEAAAAAAQAAAAABAQEAAQABAQEBAAAAAQEBAQAAAQABAAAAAQEBAQABAAAAAQAAAQEAAAABAQEAAAAAAQAAAAABAQEAAAAAAQAAAQABAQAAAAABAQEBAAABAQABAQEBAQEBAAEBAQEBAQAAAQEAAQEBAQEBAAAAAQAAAAEBAAABAQAAAAEBAQEAAAEBAQEAAQEBAQEBAQAAAAABAAEAAQAAAAEAAQABAAAAAAABAQAAAQAAAQAAAAEAAQAAAQAAAQEAAAAAAQABAAEBAAEBAQABAAEBAQEAAQEAAAEAAAABAQEAAAAAAAEBAQEBAQABAAABAAABAAABAAABAQABAQEBAAEBAAEBAAEBAAABAAEAAQEBAAABAQAAAQABAAEAAAABAQABAAABAQABAQEBAQEBAQAAAAEBAAABAQABAQEBAQEBAQEBAAAAAQABAAEBAAEBAQABAQEBAQAAAQEBAAABAAABAQABAQEBAQEAAQABAQAAAAEAAQAAAQABAQEAAAEBAQEBAAEAAAEBAQABAQAAAQEBAQEAAAEBAAABAQEAAQEAAAABAQAAAQEBAQABAQEBAQABAAABAAABAAEAAQAAAAAAAAABAAEAAAEBAQEBAQABAAAAAQEAAQABAQABAQEBAAEBAQABAAAAAQEBAAEBAAABAAEAAAEBAQEBAAEBAQABAAEAAQABAAEAAAAAAQEAAQAAAAAAAAEBAQABAAEBAQEBAQAAAQAAAAAAAAAAAAAAAQEBAAAAAAABAAEBAAEAAAAAAAAAAAEBAQABAQAAAQAAAAEAAAAAAQEAAQAAAAABAQEBAQAAAQABAAEAAQAAAAAAAQEAAAEAAQABAQAAAAABAAABAAEAAAAAAAEAAQEAAQEBAQABAAAAAQEAAAABAAEAAAAAAQABAQEAAAAAAQAAAQEAAQEAAAABAQABAAEAAQEBAQABAQAAAQAAAQEBAQEBAAEBAQAAAQABAQABAAEAAAEBAAEBAAEAAQAAAAABAAABAAAAAQEBAQAAAQEBAQABAQAAAQABAAEAAQEBAQAAAAABAAEBAQABAQAAAAEBAAEBAQEBAAEBAQEBAAEAAQEAAAEBAAEAAAEBAQABAAABAAEBAAEAAQAAAQEAAAEAAQEAAAEBAAABAAAAAAEAAAABAAEAAQAAAAEBAQEAAAEBAAEAAQAAAAEAAQEAAAEAAAEBAQAAAQEBAAEBAQABAQAAAQABAAEBAQEBAAABAQAAAQAAAAEAAQEAAAEBAAEAAQEAAQEAAAAAAQABAAABAQAAAQAAAQAAAAEBAAAAAAAAAAABAQEAAQEAAQEBAAAAAQEAAQAAAQABAQABAQEBAAAAAAABAQABAAEBAQAAAQEAAAEAAAEAAQABAQEBAQABAQEBAQAAAAEAAAAAAAEBAAEBAAEBAQEBAAEAAAABAAABAAABAAAAAQAAAQABAAEAAAABAQEBAQEAAAEBAAAAAAEAAAEBAAABAQEBAAEAAAEBAQEAAQEBAQABAAEAAQEBAAEAAQEBAAAAAQAAAQEAAAEAAAABAQABAQABAQAAAQABAQEAAQABAQEAAQAAAAEBAAEBAQAAAAEAAQABAQABAAEBAAEAAQAAAAAAAQEBAAEBAQEBAQEAAAAAAQAAAQABAAEBAAEBAQEBAQEAAAEBAAEBAAAAAQABAAABAAABAQEAAAABAQEAAQEBAQEBAQAAAQEAAAABAAEAAAAAAQAAAAABAAAAAAEAAAEAAQEBAQEBAQEAAQABAQABAAAAAAAAAAABAQEBAAABAAEAAAAAAAABAQEAAQABAAEBAAABAAEAAAEAAQABAAABAAEAAAAAAQEBAAEBAAABAAAAAAEAAAEBAAABAAABAAABAQEAAQEAAQABAAABAQEBAAABAQEAAAAAAAEBAQAAAAABAQABAAAAAAEAAAAAAAAAAQAAAAABAQEAAAAAAQEAAQAAAQABAQABAQEAAAAAAQEBAAABAAAAAAEBAAEBAQAAAAEBAAAAAQABAAABAQEAAAABAAEAAQEAAQEBAAABAQEBAAABAQABAAEBAAEBAAAAAAABAQAAAQEBAQAAAAABAQAAAQAAAAEAAAEBAQAAAAEAAQAAAQABAAAAAAEAAQEBAQAAAQEAAQEBAAEBAQEAAQAAAQABAQEBAQAAAAABAAAAAAAAAQAAAAAAAAEAAQABAAAAAAABAAEAAQEAAQEBAQABAQEBAQAAAQEAAQAAAAAAAQABAAEAAQABAQEBAAAAAAABAAEAAQEBAQABAQEBAAAAAAAAAAEBAQEAAAABAAEBAAEBAAEBAAABAQAAAQAAAQAAAAEBAQABAQABAAABAAABAQEBAQEBAQEAAQABAAEBAQAAAQEAAQEAAQABAAABAAAAAQEBAQAAAQEBAAEBAAEAAQEAAQAAAQEAAQABAAEAAAEAAQEAAAEBAAABAAEBAQEAAAEBAAEAAQAAAQABAQEBAQABAQAAAQABAQEAAQEAAAAAAQEAAAEAAAAAAQABAAEBAAAAAAEAAQEBAQEAAAEBAAEBAAEBAQAAAQEAAAAAAQABAQAAAAABAQAAAQAAAAEAAQEBAAAAAAEBAQEAAAABAQABAAAAAQABAQEBAQABAQAAAAAAAAEBAAAAAAAAAAABAAEBAAAAAAEAAQAAAQEAAQAAAAAAAQAAAAEAAQAAAAEAAQAAAAAAAAEBAQABAAEAAQABAAEBAAABAQABAAEBAAEBAQEAAQAAAQEAAAABAQEAAQAAAQEAAQAAAAAAAQEBAAEAAQABAAABAAAAAQAAAQABAQAAAQABAAABAQEBAQAAAAAAAAEAAAEBAAEAAQABAAEAAQAAAAEBAAEAAAABAQABAAAAAQEAAAAAAQAAAAEBAQAAAAABAAABAAAAAAEBAQABAQAAAQAAAAEBAAAAAQEBAQAAAAABAQABAAABAQEBAQEBAQEBAQEAAAEBAQEBAQEAAAEBAQEAAQEBAAABAAAAAAAAAQABAAEBAQEBAQAAAAEAAQEAAAEAAAEAAAEBAAEAAAABAQEAAAEAAQEBAQAAAQAAAAAAAAEBAQABAAEAAAAAAQEBAQEBAAAAAAAAAAAAAAABAAABAQABAAEAAAAAAAEBAAABAAEBAQEBAAAAAQEAAAEBAAEAAAAAAQAAAAABAAEBAAAAAQAAAAAAAQAAAQEAAQEAAAAAAQEAAAABAQEAAAEAAAEAAQAAAAEAAAAAAAABAAEAAAABAQAAAQEAAAEAAQAAAAAAAQEAAQEBAQEBAQEAAAAAAQABAQEBAAEAAQEAAQAAAQAAAAABAQABAQEBAAEAAQEBAQEAAQEAAAABAQABAQABAAABAAAAAQEAAAEAAQEBAQEBAQEAAAEAAQABAQEAAAABAQEBAQEBAQEAAQABAAEAAQAAAQAAAAEBAAAAAQEAAQAAAAABAAABAQAAAQEAAQAAAQEBAQAAAQEBAQABAAEBAQABAAAAAAAAAAAAAQAAAAABAQAAAQEBAQEBAQAAAQABAQAAAQEAAAEBAAAAAQAAAAEBAQABAQEBAAEBAAAAAAAAAQAAAQAAAAEBAAAAAQEBAAEAAAAAAAABAQAAAQEBAAEAAAEBAAAAAAEBAQAAAQABAQABAQAAAQEBAQAAAAEAAQAAAAABAQABAAEBAAEAAAEBAAEBAAEAAQABAQAAAAABAAAAAAEAAAAAAAABAAEBAQAAAAAAAAABAQAAAAABAAAAAAEBAAEBAQABAAABAAEAAQEAAAEBAAAAAAAAAQABAQEAAQAAAQEBAAEBAQEBAQAAAAABAAAAAAAAAQEAAQEAAQAAAAEBAQAAAQEAAQAAAAEAAQEAAAEBAQEAAQEBAAAAAAAAAAEAAAEBAAEAAQEAAQABAQEBAAABAAAAAQEBAQEAAQEBAQAAAQAAAAABAAAAAQABAAEAAAEBAAAAAAEBAQABAQEAAQEAAQAAAAEAAAEAAAAAAQABAQAAAQABAQABAAEBAQEBAAAAAQABAAEBAAEAAQAAAQAAAQAAAAABAQEAAQABAAABAAAAAQEAAAABAQEBAAEBAAAAAAEBAAEBAAABAQEBAQEAAQABAAEBAQAAAAAAAAAAAAABAAAAAQEBAQABAQABAAABAQABAAEAAQEBAQEAAQAAAAEBAAEAAQAAAAEBAAEBAAEAAAEAAQABAQEAAAEBAQEBAAAAAAAAAAEBAAAAAQEAAAEAAAEAAQAAAAAAAAEAAQABAAEBAQAAAQEAAQAAAAEBAAAAAQABAQEBAAEAAQEAAAEAAAEAAAABAAABAQABAAEAAQABAAABAAEAAQEAAAEBAAABAAAAAAEBAQAAAQEBAQABAAEAAQEAAAEBAQEBAQAAAQAAAAEBAAAAAQABAQABAQEBAQABAQABAQEBAQABAQAAAAEAAQEBAAEBAAAAAAEBAQAAAQEBAQAAAAAAAAEBAQABAQABAQEBAQEAAQABAAABAAEBAAEAAQABAQABAAABAQABAQEAAQEAAQAAAQEAAAABAAAAAAEAAAAAAAEAAAEBAAABAAEBAQEAAAEBAQABAQEAAQABAQEAAAAAAAEAAAABAQABAAABAAEBAQAAAQEAAAAAAQEBAQAAAAEAAAABAAEAAQABAAEAAQABAQAAAQAAAAAAAQEBAAABAQAAAAEAAAABAAEBAQAAAAAAAAEAAQAAAQAAAAEBAQABAQABAAEBAAABAQAAAQABAAAAAAABAAAAAQAAAAEAAQEAAAAAAAAAAAEBAAEAAAEAAAAAAAAAAAAAAQABAAABAAAAAAAAAAEBAQEBAAABAQEAAAAAAAABAAAAAAAAAAEAAQEAAAABAQAAAAEAAAEAAAABAAABAAEAAQAAAAEBAAABAQEBAQEAAQABAAEAAAEAAQAAAAEBAQAAAQAAAAEAAQAAAAEAAAABAQEBAQEAAAEBAQAAAQEBAQEBAAEBAQAAAAEBAAEBAAABAAABAQEAAAABAAEBAQAAAQEBAAEAAAABAAEBAQAAAQEBAQABAAEBAQAAAAEBAAEAAAEAAQEAAAABAQEAAQABAQEBAAABAAEBAQEBAQEBAAEBAAEBAQABAAAAAAABAAEAAAAAAQEBAAEAAQEBAAABAQEAAAAAAAABAAEBAQEBAQAAAAEBAAEBAAEBAAEBAAABAAABAQAAAAAAAQABAAEBAQAAAQABAQABAAABAQAAAAABAQABAAEBAQEAAQAAAQABAQAAAAAAAAABAAAAAQAAAQABAAEAAQABAQEBAAEBAAAAAAEAAAEAAAABAAEAAAEBAAEBAAEAAQEBAQEBAQEAAAEAAQAAAQEAAAEBAQEBAQAAAQEBAQEBAQEAAAABAQEAAQEBAAEAAAEBAAABAQABAQEBAAAAAQAAAAAAAQEAAAAAAAAAAAEAAAAAAQEAAAEAAQEBAAAAAQAAAQABAQEAAAEBAAEBAQEBAQEBAAEBAAEBAQEBAAABAAEBAQEBAQABAAEBAAAAAAEAAQABAQEAAAAAAAEAAQABAQABAQEBAAABAQAAAQEAAAAAAAEAAAAAAAEBAAEAAAEAAAABAAEBAAEAAQAAAQABAQEAAAAAAQABAAAAAAABAQEBAAEBAAABAQEAAQABAAAAAAEAAAEBAAEBAQEAAQEAAQABAQEBAAEBAQEBAQAAAAABAQEBAAEAAQEAAQEBAQABAQEAAQAAAQEBAQABAAEAAAAAAQEAAAABAQABAAAAAQABAAABAQAAAQEAAQEBAAABAAAAAAEAAQABAAAAAAEAAAEAAQABAQEAAAABAAAAAAEAAAEAAQAAAQEBAAABAQABAAAAAAAAAQABAAABAAEAAQABAAAAAQABAAEAAQAAAQEAAQAAAQEBAAABAQAAAQAAAAEAAQABAAEBAQAAAAEBAAEAAQABAQAAAAAAAQEBAQEAAAEAAAEBAAEBAQAAAAAAAAEBAQAAAAAAAAAAAAAAAAEAAAEAAAAAAQEAAAABAQAAAQEAAAEBAQEBAQEBAQEBAQAAAQEBAAABAAEBAAEAAAEAAQEBAAAAAAEBAQAAAAEAAAABAQEBAQEAAAABAQAAAAEBAAEAAQABAQEBAQEAAQEAAAABAQAAAQABAAEBAQABAAABAAEAAQABAAEAAAEBAQAAAQEAAAEBAAEAAQABAQEAAAEBAQEAAQABAAAAAAEBAQEAAQAAAQEAAQABAQEAAQAAAQAAAQABAQAAAAEBAQAAAAAAAAABAAAAAAEBAQEBAAAAAQAAAAEAAAEBAAABAQEAAAEAAAABAAEBAQABAQEBAAEAAQAAAAEAAAABAAEBAAABAQEBAQABAQAAAQEBAAEAAQAAAAEAAQEAAQEBAQABAAAAAQEBAQAAAQEAAQEBAAABAQABAAEAAAEAAQABAQAAAQABAAEBAAAAAQEAAAAAAQEBAAABAQEBAQEBAQABAQAAAQEBAAAAAQAAAQEAAAEAAQABAAABAQAAAQAAAAAAAQABAAEBAQAAAAAAAAEBAAAAAAEBAQAAAQEAAQEAAQEBAQABAAEBAQAAAAAAAAEBAQAAAAEAAAEBAAABAAEAAQABAAEAAAEBAAAAAQAAAQEBAAEBAAEBAQABAQEBAQAAAQAAAQEAAAEBAQABAQEBAAEBAAABAQAAAAEAAAEAAQAAAQABAAAAAAAAAQEBAAEAAAABAQABAAEBAAABAQAAAAAAAAABAQAAAQABAAEBAAEBAQEBAAAAAAAAAQABAQABAAABAQABAQAAAAAAAQAAAAAAAAEBAAABAQEAAAEAAAEAAQAAAQEAAQABAAEBAQABAQABAQEAAAAAAAABAAAAAAABAQEBAQEBAAAAAAABAQAAAQABAQAAAAAAAQEBAQABAAEAAQAAAAAAAAEAAAAAAQAAAAEAAAEAAQAAAAEAAQEBAAABAQEAAQEBAQABAQEAAAABAQAAAQAAAQABAQAAAAABAAEAAAEAAAEBAQEAAAEAAQAAAQAAAAEBAQAAAAEAAAEAAQEBAAABAQEBAQEBAAEAAAEBAAABAAEBAQAAAQAAAAAAAAABAAEBAQAAAQEAAAEBAAAAAAAAAQAAAQABAAEAAAAAAAEAAAEAAQAAAAAAAAAAAQEAAAEAAQEBAAEAAQEAAAEAAAEBAQEAAQEBAQEAAAEBAAEBAQEBAAEBAQABAAEAAAAAAQABAAEAAQABAQABAAABAAABAQEAAQAAAQAAAAABAQEAAAEBAAABAQEAAQAAAQAAAQABAQABAAEAAQEAAQEAAQEBAQAAAQEAAQEBAQAAAAEBAAAAAAEAAQEBAQEBAAEAAQEAAQABAAAAAQAAAAAAAAEBAQEAAQABAQEAAQABAQAAAQABAAABAAEAAAEBAQEBAQEBAAAAAQEAAQEAAQAAAAEBAQAAAQABAAEBAQEBAAABAQAAAQABAAEBAAEAAQEAAQABAAEAAQEBAAABAQEBAQEAAQAAAQEBAAAAAAEAAAAAAAEBAAABAAEBAAEAAAABAAABAQAAAQABAQEAAAEBAQEAAQABAAEAAQABAAEBAQAAAQABAQEAAQEBAAEBAAEAAAABAQEAAQAAAAAAAQABAQABAAABAQABAQEBAQAAAAEAAAAAAAEBAAEBAQABAQEAAQEAAAAAAQAAAAAAAAABAAEAAQEBAAABAAEBAQEBAAAAAAABAAEBAQEBAQEBAAAAAAAAAQEAAQEBAAAAAAABAQEBAQEBAQAAAAAAAAEBAAEBAQEAAAAAAQEBAQABAAEBAAEBAAABAAAAAAEAAAEBAQAAAQEBAAEBAQABAQABAAAAAQABAQEAAAABAAAAAQABAAEBAQAAAQEBAAABAQEAAQABAQABAQEBAQABAAEBAAAAAAABAAEAAA==",
"rand4_short": "DQYBAA8NBAEIAAcCCggIAw4DCQcCCwkKAw0LAwQFDwIACwkOAAYCDA8DDQQECw4OAAICCggGDQYIBwEPDgQLAAAGBQsDBgIFCw0GAQIECw8FBwkBBQIOCwsEDAwCAw4HBwcDDAkLCwAGAQcKCw4ACgsEBAYICg4BBQYOBAsCBQwCCQcMBgsHBwIDAg8FCgsKAQEMDgIPCQ8AAQALDgoOBQwCCQQLCA0CCgoGDAQBDwECDQEICAcODgMLAgQCDwYCBgUBDQkGDwgKCQoICQ4MBQ0ADAwDDAoNAwMCDwoCBgcMCQEJDg0LCAQPCQIADAUKBg8HCAUODQUDCAIOCAILAAYEAQkFDAAFDw4CBwIHDwcNCwIKCAsNDwoPDwsPAQ0MCwoBCQkIAQQCCwkEDw4AAwEDCgUJBAYKDgsNBwEABA0HCQMCDgMDBwoADAMAAwoIAQAPAQIOCAYDDAgNAAgNDwINBAkLBA4NDQEBBw4ACg8BDAcLCAkDCwcCCQ8ICQkIDw4BBAULDQoLAQIMAQUGAgsLBgwIBgQPDAQCDw8BCgIFBgMCDg0ODA4KAAgNDw4IBwQPBQEPAwwFBA4GCwANBg0IAg0KBA4JAQIICwsHBwYLAwUGCgEPAgQADAAGDQANAw4ODwMKDAYMBQoNDgMLDQ8EAgEMAgcJDQoCAQgBCw0DBAwDDA4IAQICCA4HAQQCDgcODw4KDgUEBQwEDQgPBwoABQ8KCAgGBwEMBQIKCg8ICgoEDQwMDA0NAgsCCwIICAIGCwAADAQHCQ0ODQgJCgcNBQwMAA4DCwkLAgQIAggPDAwIDwAFCAoJAgAFCAIEAw8JDQcCAwUMAQQBCQ0ICw8DBAwLBAgGCw4BCg8HBQcGCgoBCg0MCw4OAgAPDAcCDAYFAw4HAggJBQwFAgENBQoBAgYLBgIOAA0ICA4KDAEACwcMBAAEDwwOCAgPBgANDgMIBgYKAwkDDgMGAwQOBgcMBwMGCgMEBwYDCgYKAw8GCgYODg0PBgsHCwAIDAUBBQgABA0KDAgKAwUNAgAGDwwIAQgDAAUKAgMMBwoOCAwOAAwLAQ4CDA4OBQQGCwEKDAMLAAEEAQoIAQ8LAQgFCw0OBAoPDwcPBwYEBwgDBgwOCwgDDA4NCQINCAsJAAMHCgUEAg8ACwIJCgUPDg8NBgsBDw0HAQEDBAEEDw4LAQ4ODQUJCAgPCQgABA8OAAEECAQCCwcFDw0GCQwIAQUMAQQODwQJCgoIAgsLCAEJBQANDAoEBwMKCA4IBgsBAgIKDQ0CCg4CDwUJBAAEDgoJBgsKAAgABwgPBQwCCQwKAQoIDwQBDAwOBwUCBgEFDAsGBAAICQkFCQkKBAwNAwMOCQYBDgkIDgMABQAHCAgGCgEACgYPAwgNDwsACg4ICgAKCQkHAQQKDQcPAwwEDgYPDQ4ACwQJAQYBBg4PCgYPCAsPAA8BCAwJBQcHCAICDwAIDwILCwgODQENAgsKBQUFBQ4PCQYKDQMCAAsJBAwPAw4MCAQKDQgOBA0LAgQKBQsFAQkODAIHCwwAAw4BCAQMBQABCAgHDgIEAwYABgMHCAcBDg0EAAULCQIKCg8OAQ8LBAgHBQ0HBA8BCwEODgcCBAYICgsJBAsBCgYCAA0IAAMJAgAJCwwFDQsNBQIOCgEFCQsPAggGDQ8DCwsPCQkBDA4KDAQMBwsHCw0GDg0HBgcOCgENDgUODA0BBgYJCAAHDg4BBwwCDQ4AAQIHAwQACgQJAwwMBAYICg8KCg8BBAwJDAsCDQQGDAkIBAUMAg8HCwwKCAIOCgwMAQALCwwIDAsDDQwBDAYEDQEIDQQOAgYBAAMFDAIPCw4KDgQKCwcGBwoKAwQLDgcBCwoEAQUDCQEHBgkPCQ8KCQYJDg8KDQsJBQcCCAIACggGCwIHAgMBAwEICQUHCAUIDwkFDAgHBwkNAQoCBA8IAgwKAwAIBAEODwIIBggPCAsNAQ4PBA4CDwcOBQwPCwoKCgUMAgsDAQUMBQAFDwoOCQMFAQEDDQ8IBQ0GDwQKDAgECwMICQ8OAQYDBAwMAQ0HCQ4FCQQIAQMPBAENCgYKAQQFCwEACAsNDg4EDQ4BDQYBBgQADAYABQYDDA8JAwAIDgEFAwMDBAQACAADBgkJCQMHAwYCAAgNAAEOBgECDg4HBwoLAgEMDQ4PCw4NBgkFCAICBQgACAUFAgAGAgIGBAYCDQQMDgYDCQcGAAcLBwkAAAQMBg4JBQMDDQgIAQMICwUGBwkPBAkGAwYBBQwLAwgNDgYIDAgGAwQPDQsDDAAPCAEDDgwJBwIJAwUDDAgBBgIHBA0KCQUOAgoECQIBCwYKBwgNDgEODA0GDw4HAwsMDw8NAwgJDg8LBgkMCgkBBAYIAgQFDQoGCwgGAw8ICgQICgMGDw4IAQQBDgEKCQMLCQkACwkGDwAFAw8MAwsCDQIECQgHDQoIBQINAgEBCAoHBAUABAsKDgEFDQkPDA4PBAAEAwwHDQQNAgoMCg4PDwwGBwkADA0HCwYJAQQBDgsBCgwBBwIPCwsLDwkABQgKCwgBAQMNDAEPBAcGBQQOAwwEBQQNCgoCCQgMBA8DAwAMDwYOCgsGBQcGAAEOBQoODA4NCwUKDwwCDQ0LAgoPDAYNDw8DBAQMCg0LCQkKDwkFCQsABgUICgEPBwgADQAOCgIECAMCAQgBAAoOBAcACwwNBAMADgYLAAoCDgEBAwoDCQsDBAoPCgkDCgsPCgYPCAAEAgsLBQMLCQ8DAg4NAAcOAwYLBAoABQ8FCQEKDAUIDA0JDwwODwEMBA4NBgkPAggEDAQIAQcHCQoHDgsKCQ4OAwQHCggEAAUKCQ4DAAkLDwgIDQ0ICQEAAAEDAwQJDwcKCQcADgkNBQEPDwkJAgUPDAcDDQAHDA8EBwUEDgMPDQ4PCA0ODAgPBQQEBwkJBgQHCwMNCwUOBA0IAwEGAAkNCwgEBAsIBQcNCA4JDw4IAAYGBwMFBQcCBQEODAUCAQgABAkIAQILBwMFCgAAAg8KAQYBDwQMAgUNAQIGBAoEBwcPAA0OAgAOBgAOBQMMBgwKCAMAAAAOCQoJDA4ODgoDCAQOCg4LBQAJDgkFBQ4ACgwADQ8NAAoABQUJDwIOAgALBgcNCAUIBg4HCQABBwAJBwsIAQ0LBgEMCgYPAAoJDAkJDAEJBwwKAwEMCwgJBAYFAQ8DDAsMAwoACQACCQIIDgwHAAkHAQ4GAAMJAgEPCQgHDgIBDAQGCwwODggIDAoDAgMNCQgFCwoFCwgEDQMLCQoGDgQMAQkGAgYPDA0ABwgNAQ0NBg0DAwoOBQAFBwkNDwwADQMKDwUGCggBBgIFAg0BBQYHBgMKBA0ICwoJDwoJBAwLAAQEBQICCwEMDgwPCQUFDQIJDgYFBQQPCwkICgMJBgoMBA4ECwYDDgUHCwYIBw0JDA8BDggABA4IBAUMBwEMAgwHBwsPDgQCDgkIAg4JAAcDCQwHBwAKBAMKDgALBwEDCAgNAgsECAkODw8JCQUIDgcNCgAJAA4ACwELAQQGDgcEBAwOAAcFAQgBDAsHAgQMAQwAAgwNBAECAwMFBgYFCgwECgMMDwENDwMGBA4ODg0HDgENAAsICQIODgEGBQ0EBQIBDwYBCgEGBwMAAQ0CBQQMCAoABwUABwwJBQcJAwkHDgAJBggOAgcFDwgLDgMCBgMIBAUNBgkMCwkFCQ4NCAsGDw4DDwsKCgsLBQADDAQGDAEODQMJCw0BAQgMDwcNCgkJBAkFAQ8HAwgOAwYMCQwBBAUHCwAFBQQPAQkEDAsKBgQNAAgJDQcHBA8PDAMCCAsMAA0MBw0CAQUGBggOCQ0BBggLDwACDgMHAAcGDAANCAQDDQ4GAA0HCQADAwAPCQ8HAAYHBQIODA8CCgwODA4GDwAPAAIMDgcODAgOAgMGBw4EAQQGBQYKDw0BCAgKDAUNCAcGCAMBBAoIDAgICggBCg4JBAoGDQgHAgsCBgcMCgwADQUGAwMEAwYPAwkPBQ0BBwAKCgQJCg4HCQALAg0OBQQABg8EDAUCDwsBDA8FBQkECQABAQUPDAIJDgMJDQ0CDggNAg0BDw4BBwcDAwQJDwYCCQoDCQ8KAQUEBQAIBQsNDwoJCQsDBw4BCQoDAQ8ECQEEAwoBDAANCwIFCgUGAQADBQQEAgIDAAgDCQEIDgQIBwEPDQIMAA8LBQUJBQILBQUCAgAOBAAPCwAMAAsCAwcCBQQMDQ4GDw0CAwkFCwsOBgUNDQgMCwEFAAILBQsMDAwKDQYBBwcHAggKCwEKDQEGDQkDBQQMAgsICQQCCgIBDwkEAgcEAAwNCgYCCAINAgALCwsMAAMIAg8LCQsCCwkMDA8KBQEHBwECBQkLDAoOAAYICQ4BBwgEBwAFBgwIDAkBDgYGCQUBBQsIAQEGDQsFAAAJDgYKCQ0DAQ4GBAEDBgEJAwoIDgkGDgkPAQQMAwkDBAMOAQIKAQIDCQ0ECQ4ECQcMAQwBBQcFAg4ABQcODgcBBw4IDQwFAAQGCAAHBQoABAgBBAUMCQwJCQAJDggNCwwCDgUBDAgHBAIFCQYNBwAJCA8PCwgICQsMAgcLDgYEDgkLAQsCCgAPBAAKDwQJDQoADQABBQIBBwgNAAUOCQ8FCgYPDQMOBgUICAYGBAUODwMMDwsLCgECAQAOCAwMBwsIDA4ICgMGDQUMBwwKCg8ABAgACwwEDwkJCwsGCAULCw0HAAgMAgkEDQYKBQAJBw4HDgsCBQkGDAkFBwAIBAwODQgFCgILCw0MBwAICwoIDQIDDQABAwwACgoEAg8JCQsCDgIHCAcBDgUDCQQOAAsBBQQBDwEGDg0BDAMJDwUPCgQDCA8LDQcGCgMNDw8CDwQFDAYFDwgMCw4CBAkNBQ8HBgAACwMHBQYHAQYHDwYFBwgHBAEJAQ8OAQcEBAkMBQEBCw0OAgQNBA4JCAwFCAoNBgkDCQIEAQ0HAQoLDgIFCQoPDQkPDAgCCgcLBAsGAwMHDgsCBwsFDQ4KCQ4JDQMHCQwEBAYDDQcABA0FBgUABgsJCAYMCgEBCwQOBgECDwgEAw8HAwgBBQcPBAIHAgwPBAkIBwQBAQgEAQ8IDQoODQ8KCAUGDw4DBAMJDAwMAQIEAQkCDgoODggFDQEMBAUCBw4HDwUNBA4KAQEFAwsADAMKDwEPAgkBDwEJAgYHDg0MCwIOCAcGBAINBgIFCAgAAwoBCAUDDA0BAggCBAwBDgYMBgUGDgQFDw4FBAsFAQIHCQsGBA0GDg0CBg0CAQkLAQUJBw0IBQQIBgIHAgkDDAQKCwcOAAsNAAMFBQgFCAUHDgUPAgwLCgQDCQMNDg4MCAIFDQ8LBgkKDQkBAwkKCAQOCwAHBg0GAw4IDAoLCgQNDAoNBgMDBAQJDwwEBwIGCgUPCgcBBw4KDw8IAAEIAwYEBwsBAQ4MAwgOBwIIDAENBwMHAgQCAQ4JBQQCBA8OCwkLDAAPDQgIAAAGBgIMAQMGDgUKCwwPBgIBDgsIDQAJAwMIBQ4GAQcBAAMDAgMEBgUEDAcKDQ4EBA8PAQAHBwsKBQYPCA8ABQ8EBQQOAQAHCAMFAg4KBwsGDgwGAQUJCQUODwYODQoJAgcACgoLBQUBCQoODAcDDgcPCw0JBwgOBQ8DCAYOAgYLDwEFDgYDDgoGDgECBQcMBAgCBwIDBQgHBw0EDQsOCAgPCAQBBgcEAg8LAQQPAA0HCQcIDQ8IAg0ADAwDCQMBDAYCCA0CCgANCwoOAgYDAQoEDAoDCw8KDQkACwMLBwoJAQUKCgsODAwNDwsNBwcBAwMFAgcJBgYJAA0KAwkNBAYFAAsBAA0EDQ4JDAgNAgMIDQQFCAYHAwANDgQGAAMJCQgJAgoBBgYCAgsIDAkJCwMGAggBAQIBDQYGBgkCDAwPBQMODwoCDwIGBgcOCQINBgMACQkBCwANBwMBDgsJAgkBDgEDAgwODQQLAw4HBgcMAQwBCg4ABwEBCQkJDw0DBAkFAQwBCggIDgIBCg4PBAkCAA8GCQ0HCQELBAQKBgsLBgEKAQAKAgYFBgAOBAcOAAMOCAEMCgwLCAQECwsKCAoEBwEGDQkCBQgBCwAEAg0PDgMBBgcFBgoJBwUCCgsGDwcADgoFAwYKBA8FDA0LAAICDQ4GDwYOAQsKAwoEBgACBAIHCAAMAQ0PDgwPAA8IAwIKBgQKDAYBBgUBDQwHAQcEDwsPCAwLDwsPBQUGCAMECwQODwoCDAkPDAMDAAwACw0FDg8NAQMPAgQHAw4LAAMEDgUIBQALCAcFAgIBCQAIDQMGBA0ICQ8MAwALDgkMBAkDAQsKAQ8DDAsNDAAJAwIPAwQIAAYECQkPCg0FDgsFDw8DAgYDDwQKBw4ACwAIBgcNDgUPAAwPCQMMDg0DBgUMBAANBg8CBAYLBQoJCgcOAQwDBgQIAQ4NCgEKCwEKCggPCQANAQgFDwkLAwUJDwoNAwUABQQACA8NCwMGBAMLAQ4CBQQHAQoNAggCBQ8DDQQBBwoLBQ4FBQgAAAQJBA0FCQMGDgQLAwENCAoOCQUODw0GBg0OBA0IDQIHBwIFAA0HCA4FCQQKBAsNAA0NBAIDCQ8DAQ4DCAIPCQIEAwYNAwwBAwcPBAkEAwMJBA0HDQALAwYJCQkPCAYFBwoADQgHBw4PAgoDBwoMAwoKBgcACgEGCw0CDgMBDQQDBggCCAoNBgkFBAwBAQ8IDgwHCw0IDgEDDAkCCA4HAQ0PDgkICwAPCgcGBwMCBQsPBgALBQsMCAUCAAAECwcMCAwGBQUICgIHDQEBAAsPAgUMBgADBgQCAQMKAgIJAAQPDAsNAgcFDwUKAQcECg0IBAUDBAEIDQ8OCQkJDQsOBg0JAQADDAYLCwwICwIGAwkCBwUBDQkFBAkDDwcGDwMNDAMGCAsEBQMHDgcFAA0DAA8GBQUOBQUIDwoBCwwHBAIACw4NDwsCAg4FBwYPDAYJDg0DDQUOAggJBwUPAQ0JDw0HBgoBDwcPDgwIDQYLBAkDBAANDQgODgwKCgcPBQMBAgoDCgwCDg4LDQgEAwkBAw0CDgsHBA0IAggFAQQNDw0JBAUDAQcHDwQDAgYCAwINCg4JBQEADw8PCQAEAgwFDw4GCw4NDQoCBg4EAAUNDAIKDw8IBQQCBwwJCA8KCgEEDAgCDQoPDQUCBQUNDQMMDgoLAwQCCwYICA8ICQsDCwAODwwLAAEEAwIEDgoBAQIFDQYPDQAJBQQNDwsABwMEAwMFCgUGDAsGCAUADAwIAAsNAQYDCgQGBAINBAQJCwIPAwkIBgsKCw8FAw4MCAMPCQUGCA4FAAMJAAgHAAUIAwICBwwGAwQCDAwCAQ4PDQ0ICggGCwgBDQYBDA8JAQkIDAwEAA4BBgcACwMBCQgODwsLBwkMCA8DBwALCgsICwsKCg0LDAwABQoMDQMIDggKDQMCBgYMCw8DCQEKCwQJDA8PAAIGAAoCCgwDDAkGAwwICQ8NAAMOAg0DCAsGBAkPCg4MBgoGBwsNCQUCBAMFBQMFDAoFBw0GDAwPDAQNAQwHAgEJBgwMCwwCDQQEAgsLDAIHAwYKCQAEBwYPBQcHCgQKBwgDCAUIAQYNAggAAgUDBQ0OAA0FBQEPCgwLCgQCDgIBDQ8EAAQJAAAPBAcMCggEAwcEAw4LDQQFDwwEBAACCw4IAg0EAQ0CCgQKDgYGDQAMBgkDDQwPDAMEDQMFBAoMBQcNAQoJAwsDDAwCCAYLDgENBg8JCg4GCA8IBAsACwECBQ0DCwwOCA0GDAYJAAIOAwQBBg8DDQMABQMECA4CCggGDQAGCwAKCwMGDgoCCgUFBgYMCQcGCwEGAgwECQYEBQ4ADgAGBwcBBA8LBgoFDgAPDQoNDwoCCAYMCgkMCQsIAQ4ADgMGDw0DCwkDDgEEAA0GCwQFDAYBAQ4IBA0GCw4HAgcADQcFBQUFAggGDgMMDQQCDgYMBgIDDAcLBg0ACwoHCwcLCgwBAwYOBA8MCA4DAgYBAgEEBwkNDAkCDwYJDgcHDQIICA8ADQ8FBQgDCwMBBgsMCQAHCw4PBQ0OAAMDDwYKDQIABAANBgMNAwQDAgsGCgICBQcBBg8ECA8LBAUKCwsHAAwHDAwJDwINDAwLBAEABg8DDQIIBQ8LAwsPCQMADAEDDgMMCA0HBQAPAwAABAkCBgsFBAUHBQMNDwkFCwQNCgwCCw4JCA0ODQEACgsPBgELAA8JAgoDCQcPAQYDBQAPAwMLBgYNCAYLCggBCAQKCwwLDQsHBQ8FAQwADQsKCw8JAQUGDgMEDAMDBAkDDwUNDwoIAQsFAwoPDgcLAwIKAA0ACAAFAQABAwwPAQUEAAAMBgwEAwUBBgMHBwsCBwcBAAsLDAkOCgQJBQQLDgsPBAoJAwwPDQ4OBggDCAcBCgcOAwcBCQUJCgEKBgUKCQ0DAAEOCgkBBQoIDAUJCQEIAgECDwsNDwMIAwQDDQ8FAwIJBQENAQcICwoHBg4PBA0DBwEMCAUKAAYMAwsOAg4MCwEIAg4CAgsHAgMMBA8EBAIKAQcJBwoIAw8LDQINDQ0CAAMIDwQHBAwICAwOAwgFCgwBBA4PDAsPDw0JCAsOBAYOAQYBDA8IBQYMCwYCAwUEBAgIAQcDAwAFDAIGBAMMDQsABwsIAggEBgUIDQYPDwMLBQoICAAKBgYECw4IBwIIDgQDDwIDAAgMDAgNAA4ICAQGCQ8MBg4ABwkHDQsPBwMICgEDBg0FBg4JCg4IAA8ECQsKCgYJDAYMAgELDgcPDg8HBwwFCg8KCg0GCAIDDg0NDgoLCgEECAMOCQwGDg8JCAsFBwgMAQkPDAIDBg0KDgYFDgIHCwUIDgABAQ4DCwgKBgIJCgkCCg0AAggCAAEECQ8LDgwEBwoCDAMLDQYIBgcHDwsCBwgGBQsEBg4IBgwHDg4JAQAMBQYJBQAGBwUKCggOAAAPAwADDAMADQwDAwQOBgAECgENDgIOAgsPBggOAQQIBgYDDwINAAMLDgMACAYBCQAFBQ4NDwYJAwQLCw8LDQQPBQALBw4NBAAHDgUKBAkICwMHDwQKBgMDDgUJDAIIAwoFDwMADwIKAQMHCQcODgABCgMIDQ4CCQcPDwELBQALCAcJCgcGCA4BAAUMAAIODwUGDw0DCgACDgsAAwQJBAsEAwkMDwAEAAwJCAkBCQMHAwgMAAMGBggHCwMFDwAMAA0ODAoBAQsBCwgFBwkICg0AAAYKBw8MAAMECgQBDgkLBgMLDgYFAA0JDgUCBg8MCwcOCAwLCAcPAgIDDgEICgsIDAoBAgAMDQIICgcADAUDCg4PBA0OCgYDCQoGAAQDCwQADAIAAAEKBAUMBQkHDgsOCQIBCwsGBAcHDgIJBgIDCgUKDQADCggIDAwGCAAHAgAKCw4HAw0FDwwMCwkNCgcGDgoMBgkADggLBw0DBwIICAsDDQwKDgQLAQABDQQODA4GBgUBBwAMAgwBAwkGAAIDAQcACQcLDgcLBAMNAgkICgABCwwKBAAMBQ0IAwwNDg0FDwECDwgHCwINAAkFAAQKDQUBAg4ADgQODAENBQ8MAAgECwgACwcJCAgHDQ4BDggBCgQBBwsGCQABDAYPDQkJAQ0MBQwDDwUBAAUDDgALCAUGBAYBBwEEDwYPBQAFCg4IBwUOCQkNAAgDBwQMBwYBCwAEBQYABwIEAwoDAQALAAMHBw4MAQcBAAYBBQMPBQUABQcPBAoJCAwHDQcFAgYODwoPAQEPBwIMCA4HCwgOAwwFBgcGCAIPBQMGBAgBDw0EDg0ICw4FBAgBCAQDBAIOAQQMCQ4JDgwJAAoMBQUEDAINBgoOAgsDBQMMAw0OBAsPDggGDgsJAgUOCQcNCggNBgIGCgsGDg8EAwIPAQUJCwEKCwcIDA0KCQ8ACQoBDgAIBwcHAwICBQQJBQ8CBQoGDg0ICgwOAAsKCAkNAwYACgQFCg8GAA8CCQgNAQUMAgkDDgwEBgoEBwsJBAMLBAEMCwcDDAYLBwYOBwwJDwoOAg0PBwEKCgsCAQoEBgQJDQEDCwAIAwkIDwwJAgQJBA8BBw4BBQQACQoDBAsEDQYBCQcPDwEMBwsNBwYHAgsBDAEACgAKCQQABQgDDAMEAQ0GCAgJDQsLDAILCwgFCA0FCQcEDwAMCAkMBQQAAQMHBA4CAQcOCQANBw0KAwwKAgQIDgcBDAcLAAAODAgIBQ8AAAkEBg4MBAgEBAwDAwYGDwIECwYDDwsIDwcEBwMDAQMODwwOAAMPCA8LAwYCDgYBCA8GCQ0DDwgIDAAFBwMLAwgIBwIIDQ0ICwQGDgQMCAoADwIHDAgKDAgNAA4LAgwPDw0FDAMFDQ0NCgIIAgQCDAILAAILCgwPAwINBg8EDg8ECwYEBgQIDw0PDA8ODgYJAQkIAAkNBQQPCQ4FBgoCBAkGDggLAwEJBwYBAw8ODwQOBQkLCwsKBAYGBQAABgINCAYPBAwCCwkCAQ8LCQEGAggFAAUODQgLAAYFCwIEBAYMBQYDAAwKCQYAAgcGCgEKCwMOAwkBCgYKCQENDwICDA4OCgwNCwMDBwwFCAkEDwgHAg0NDA4LDQ4GDwQJCgoMAAMAAAgOBAUCAQoJAQcHBwoNBwUKBAECAA0NCwcGCwQBDgQHDwsMDggCBwYCAgwEDQQFCgwCDAwCCQsEDwsCDgMBBw4KCg0IBQwABgMNCQgPBwcMBwQCCAQHAgQICg8BDwYEBQ8BAgkFDwkFCAMACwoHAwkEBAsFAAILBAkOCwUNAAYHAgkBAg0KBgEIDAIPDQMAAQ0HDAUCAgkJBQMECQIKCgIJBgEIBwQFCgkPAQoCDwoIAwoBCAQGDQwBDQUFBwMECgsABAsDBQIMBAwHBwwPBA4MCAMHBQ4ACw0LBQUNDwQLBQcMAgMBCg0IAAECCA0OAQoGBQkCDwUNBA4GAgEGDw8BCAkNCQoKAQcKAwoKDwADCwoMCQwOAQYMDAsMAwELAQMCAgQGDgoGCwgNDg8PAQIDAA8LBA4PBAQDDAEJBAIPBw0HDwAFAwkBAwEHBQcECQgBDgELCAAAAQwNCwwIBAgBDgUJDwoCDAcKDQoMCQgOAwgHBwkLBw0BBAIFDQILCQUECQQCAQUKAwkOBgIFCQsIAwgNBQ0KBw0PBQIDDQkKBggFBQIIBQMECQ0FBg8GBAIAAwMLCQYDCA8FBQcFCAcDAg8FBwsCBQAPAgYPBg0GAwQJAAAHDwIICggFBAcLCw0EDAYIBg8NAgkMAwQGDAgPDwYPDA4DDgMACgIOBwcKAQMNCw8ABAQMDQQCAwQGCAYEBQsBCwcDBQAJCAkPCgsACgkAAQcDDgUGAwcNBQsHDgQICw0ECQwJCwYGBQsEBAoNAwIHCQ8BAQcHCwgADQkKAQYDCAEDDw4FDAkCAgQHDAAIAAYJBQwEDwsDBAwDDwsIDAYDCA0AAA8MCgUDBAMBBgQICg4DCwoFAQkFBgIGDgoBDQwMDQIGAwABCAoFBAIABA4MAgUKDwcEDA0CAQcEDAUPCwwCAQYKDgELCwMBDQMMCAYJBwUADgsJCQkOBgIAAQoIDAQFAgsPBQsKAw0FCggOAwsOCAoMBQsGCwUDAAEOBg8IDgkADw4OAwkEDQMNBQoEBAIBDg0EBggNAwQHCA4FAQoEBgwLDQYDBQMNAgILDQ8EAgQMDgsNAgwDAwEJBQwHAw8DCAABAwAPCwYJAQ0KBAAACggABwUABw4JCw8OBgYACAYEBwgGDw0OAgcBBQYCAwIMAQMACwUEAAwBBw4FBwAKAQUNAgAFAggMDQIGDwkPDgwMBAoPAgkPBQkBBQwBBwcMCQUNAgMLBQoMDQYGCAgOAQcOBAkACwYEAgsKBggIAgsNDwEJBwgIAgYECwQIBQEBCAMOCgYBAQoCBQgGBgsHCgYKAw8HBQcADQ4EBAUPAQsJAwQJBgMACgwLBAcPDQMGBw4JCgQDBAwEDAYGDwAHDAoBCgwACQkEAQ0MBgoFBggFDwECDQQMAw8DDA0KAgsHAAEDDwECAwgNBwUIDgMNAAIIBAYKAAcMDQsICg4FBwUJCAkNAQMNAgoKAggEBQ8HAg4LDAUNDQQCCQoKBA4KCAANAQsIAAUHDA4LBgwHAQABCwgABwsMCQoNAwkPDAABAQANAA4GBA8KDQcACAwHAAgCDQAICQwKBw4JAwUKDQkBDA8DBA8IAgsCDQcJBw8NAQAKBQgHAA8PDgcEBAMLDA0GDw8FBQEIAAELDwIDAg8EDg4MBQgJCwkBAgMCCQAODQYIDAgLCA8OCwgADwgOCQIJDgEJCAsECwMFBQgAAQMFDQADCwQLBQYFBgMGAAEACQwMBwQHDwEADAAABgQIAw0HBAEMBQQPAA0DAQUNCAgFDQ4ABgUODgACCgEHCAAGAQwIAQQMDgMNBQsMBwQFCA0LDwEIDAoKAg4LAAoMAQwIDQ8LCwEGBQsPCAUICwEPCQYKCAsBCgoGAQ0BDwIABw4CBAgKDAoGDg8PAQMADA8IAAAFBw0FAgMHAgIEAAYCAAgDAA8MCQQAAgwCCAUPBwMLBwIIDAoMAgsLAAMABwgAAQ0EBwULDQoOAwULBQUDDAoLAA4CBQ4EAg0FBQEEAAEJAA0JAQUFDwECDQwNBwYMBAQMBg8NDw4MCgwFDgMODAMFCAANBgUPBQMLCgAPDw8BCgEADwYKBAQOCQIPBgwABgwEAgsGAgACCwIHAgkADgILDQ0JAgsMCgEEBQ4NAgQADAkFDAgMCw0MCQUGAgoGCwUACg4IAQcFAAgACgQMAAwDCwkFAgIADQQKCgILBgQCDQsEBgcOBAQBDwUIBAYJCgEOBg4KDAgMBgIEDgoEDQcPBAYFBgAACgMBAA0OAQIBBAwKBgcIAQABCwIBAw4ADQkMAQsCDAoAAQ4JDwwNAQkKCgkPCQgHDQkOCwALAwUDCwsICwYNDwwJAw4JDwoICAwFCAMOAgcHBgQCDgkIBwEGCgEHDQQOCwEPDgsODwYCBQwDDQ8NDQ0DBQkJDwIBBwkLAwsJDQoLDQQHBw8ABQMCDgQOAAwCCgcGBQoBBgMMCwYPAA8PCAkCAAEFDQgGAQIIDg4EBwcDBgAIBAAEAwQNAAIPDgEKAwcKDAEBDg4GCAgPAAgCDAoFBgoNAgIJDAMPDAMDCg4DBg4FAAcDBAAFCg0OCgkOCQgMCA8MCwUCCAkKDQIHAAQLDgYNAQgOBAUNBAEGAwAJCw8ACwYDBAkGDQMECgYJBAIOCQ8CBQwNBwcMBg8OCQYDDw8KAAsHDQcGAwIBBQgIBQYFBAsKDg0MBQwFDQANDgIJAg8OCAICDAsADAsNCwYPAA==",
"rand8_short": "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"
}
EMBED_REFS = {
 "rand1_short": {
  "sym": {
   "MCV": 0.9610588257005508,
   "Collisione": 0.6914641209972466,
   "Markov": 0.987596104459409,
   "Compressione": 0.6117162047939451,
   "t-Tuple": 0.8676244308170737,
   "LRS": 0.9626258038327872,
   "MultiMCW": 0.9526180605326265,
   "Lag": 0.9433337065757712,
   "MultiMMC": 0.9616166678288808,
   "LZ78Y": 0.961446245952478
  },
  "bit": {},
  "final": 0.6117162047939451
 },
 "rand4_short": {
  "sym": {
   "MCV": 3.790037390213974,
   "t-Tuple": 3.5674726723995995,
   "LRS": 3.833525522232983,
   "MultiMCW": 3.8669546824826764,
   "Lag": 3.7836506125537386,
   "MultiMMC": 3.884655279493375,
   "LZ78Y": 3.882495664055963
  },
  "bit": {
   "MCV": 0.9791894829624022,
   "Collisione": 0.8981794483810188,
   "Markov": 0.9906168077709475,
   "Compressione": 0.8038720669691841,
   "t-Tuple": 0.8987772293903774,
   "LRS": 0.9329693144953363,
   "MultiMCW": 0.986560864592317,
   "Lag": 0.9826418217359891,
   "MultiMMC": 0.977696512087203,
   "LZ78Y": 0.980145173356401
  },
  "final": 3.2154882678767365
 },
 "rand8_short": {
  "sym": {
   "MCV": 7.010454037736041,
   "t-Tuple": 7.010454037736041,
   "LRS": 7.289198671720653,
   "MultiMCW": 7.375192249729911,
   "Lag": 6.636441287083911,
   "MultiMMC": 7.327627679188088,
   "LZ78Y": 7.353355393835384
  },
  "bit": {
   "MCV": 0.9833867846591503,
   "Collisione": 0.832052982215248,
   "Markov": 0.9977249767279653,
   "Compressione": 0.7326117180606562,
   "t-Tuple": 0.9107864457354141,
   "LRS": 0.9819303577363743,
   "MultiMCW": 0.994537115514506,
   "Lag": 0.9896934938879549,
   "MultiMMC": 0.9878145440422941,
   "LZ78Y": 0.9880818037999031
  },
  "final": 5.860893744485249
 }
}

def selftest(quiet=False):
    all_ok = True
    report = []
    for name in EMBED_DATA:
        data = base64.b64decode(EMBED_DATA[name])
        mx = max(data)
        bps = max(1, mx.bit_length())
        ref = EMBED_REFS[name]
        res = analyze(list(data), bps)
        checks = []
        for key, rv in ref["sym"].items():
            r = next((x for x in res["rows"] if x["name"] == key), None)
            h = r["h"] if r else None
            ok = h is not None and abs(h - rv) <= TOL
            checks.append((key, rv, h, ok))
        for key, rv in ref["bit"].items():
            r = next((x for x in res["bitRows"] if x["name"] == key), None)
            h = r["h"] if r else None
            ok = h is not None and abs(h - rv) <= TOL
            checks.append(("bit " + key, rv, h, ok))
        f_ok = res.get("final") is not None and abs(res["final"] - ref["final"]) <= TOL
        checks.append(("minimo finale", ref["final"], res.get("final"), f_ok))
        if not all(c[3] for c in checks):
            all_ok = False
        report.append((name, len(data), bps, checks))
    if not quiet:
        for name, n, bps, checks in report:
            ok_n = sum(1 for c in checks if c[3])
            print("  %s (%d simboli, %d bit): %d/%d ok" % (name, n, bps, ok_n, len(checks)))
            for label, rv, h, ok in checks:
                if not ok:
                    print("    NO  %s  NIST %s  calcolato %s" % (label, fmt6(rv), fmt6(h)))
    return all_ok

def fmt6(x):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "—"
    return "%.6f" % x

# ---------------------------------------------------------------- esecuzione

def run_file(path, args):
    with open(path, "rb") as f:
        buf = f.read()
    o = {"fmt": args.fmt, "bps": args.bps, "offset": args.offset, "skipBad": not args.strict}
    ex = extract(buf, o)
    sym = ex["sym"]
    if len(sym) == 0:
        raise ValueError("Nessun simbolo estratto.")
    if len(sym) < MIN_N:
        raise ValueError("Simboli estratti: %d (minimo %d). Il calcolo non è definito." % (len(sym), MIN_N))
    distinct = len(set(sym))
    if distinct < 2:
        raise ValueError("Alfabeto di un solo simbolo: nessuna entropia da stimare.")
    res = analyze(sym, args.bps)
    return buf, ex, res

def verdict_of(res, n, bps, min_h, selftest_ok):
    reasons = []
    h_bit = res["final"] / bps
    if not selftest_ok:
        reasons.append("l'autotest del calcolo contro i riferimenti NIST non è superato")
    if n < WARN_N:
        reasons.append("campioni estratti: %d (il NIST richiede almeno %d)" % (n, WARN_N))
    if not math.isfinite(min_h):
        reasons.append("soglia di entropia non valida")
    elif not (h_bit >= min_h):
        reasons.append("min-entropia per bit %.6f inferiore alla soglia %s" % (h_bit, min_h))
    return {"ok": not reasons, "reasons": reasons, "h_bit": h_bit, "min": min_h, "n": n}

def paint(text, ok):
    if not sys.stdout.isatty():
        return text
    color = "\033[1;97;42m" if ok else "\033[1;97;41m"
    return color + text + "\033[0m"

def main():
    ap = argparse.ArgumentParser(description=VERSION)
    ap.add_argument("file", nargs="?", help="file grezzo della sorgente")
    ap.add_argument("--fmt", choices=["raw", "bits", "hex", "dec"], default="raw")
    ap.add_argument("--bps", type=int, default=1, help="bit per simbolo (1..8)")
    ap.add_argument("--offset", type=int, default=0, help="byte (grezzi) o righe (testo) da saltare")
    ap.add_argument("--strict", action="store_true", help="errore invece di saltare le righe non conformi")
    ap.add_argument("--min-h", type=float, default=1.0, help="min-entropia minima per bit (predefinita 1.00)")
    ap.add_argument("--out-dir", default=None, help="cartella per rapporto e sorgente estratta")
    ap.add_argument("--selftest", action="store_true", help="esegue solo l'autotest contro i riferimenti NIST")
    args = ap.parse_args()

    if not 1 <= args.bps <= 8:
        print("Errore: bit per simbolo fuori da 1..8", file=sys.stderr)
        return 2

    print(VERSION)
    print("Autotest del calcolo contro i riferimenti NIST ufficiali:")
    st_ok = selftest()
    print("  Esito autotest: %s" % ("PASS" if st_ok else "FAIL"))
    if args.selftest:
        return 0 if st_ok else 1
    if not args.file:
        print("Errore: indicare il file da verificare (oppure --selftest).", file=sys.stderr)
        return 2

    try:
        buf, ex, res = run_file(args.file, args)
    except (ValueError, OSError) as e:
        print("Errore: %s" % e, file=sys.stderr)
        return 2

    n = len(ex["sym"])
    v = verdict_of(res, n, args.bps, args.min_h, st_ok)
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", os.path.splitext(os.path.basename(args.file))[0]) or "sorgente"
    out_dir = args.out_dir or (os.path.dirname(os.path.abspath(args.file)))
    os.makedirs(out_dir, exist_ok=True)

    lines = []
    lines.append(VERSION)
    lines.append("Rapporto di verifica — stima della min-entropia, non certificazione della sorgente.")
    lines.append("")
    lines.append("File grezzo: %s" % os.path.basename(args.file))
    lines.append("Dimensione: %d byte" % len(buf))
    lines.append("")
    lines.append("Estrazione: %s" % ex["source"])
    lines.append("Formato: %s   bit per simbolo: %d   salta: %d   ignora righe non conformi: %s"
                 % (args.fmt, args.bps, args.offset, "no" if args.strict else "sì"))
    lines.append("Righe usate: %d   scartate: %d" % (ex["used"], ex["skipped"]))
    lines.append("Simboli: %d   alfabeto osservato: %d" % (n, len(set(ex["sym"]))))
    for note in ex["notes"]:
        lines.append("Nota: %s" % note)
    lines.append("")
    lines.append("Stimatori (H per simbolo):")
    for r in res["rows"]:
        lines.append("  %s: %s" % (r["name"], r["status"] if r["h"] is None else fmt6(r["h"])))
    if args.bps > 1 and res.get("bitRows"):
        lines.append("Stimatori (H per bit):")
        for r in res["bitRows"]:
            lines.append("  %s: %s" % (r["name"], r["status"] if r["h"] is None else fmt6(r["h"])))
    lines.append("")
    lines.append("Min-entropia finale: %s bit per simbolo" % fmt6(res["final"]))
    lines.append("Autotest del calcolo vs riferimenti NIST: %s" % ("PASS" if st_ok else "FAIL"))
    lines.append("Min-entropia per bit: %s   soglia richiesta: %s" % (fmt6(v["h_bit"]), v["min"]))
    lines.append("ESITO: %s" % ("SORGENTE UTILIZZABILE" if v["ok"] else "SORGENTE NON UTILIZZABILE"))
    for r in v["reasons"]:
        lines.append("  motivo: %s" % r)
    lines.append("Riferimento ea_non_iid: %s" % REFERENCE)
    lines.append("")
    lines.append("Avvertenze: il valore misura la min-entropia dei dati caricati, non la qualità della sorgente nel suo complesso.")
    lines.append("Il file della sorgente estratta è candidato alla concatenazione solo dopo la revisione di chi lo usa.")

    report_path = os.path.join(out_dir, "rapporto_%s.txt" % stem)
    source_path = os.path.join(out_dir, "sorgente_estratta_%s.bin" % stem)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    with open(source_path, "wb") as f:
        f.write(bytes(ex["sym"]))

    print("")
    print("Estrazione: %s — simboli: %d, alfabeto: %d" % (ex["source"], n, len(set(ex["sym"]))))
    print("Min-entropia finale: %s bit per simbolo; per bit: %s (soglia %s)" % (fmt6(res["final"]), fmt6(v["h_bit"]), v["min"]))
    print("")
    banner = " SORGENTE UTILIZZABILE " if v["ok"] else " SORGENTE NON UTILIZZABILE "
    print(paint(banner.center(60), v["ok"]))
    for r in v["reasons"]:
        print("  motivo: %s" % r)
    print("")
    print("Rapporto:  %s" % report_path)
    print("Sorgente:  %s" % source_path)
    return 0 if v["ok"] else 1

if __name__ == "__main__":
    sys.exit(main())
