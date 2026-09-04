#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
農業部水產試驗所「衛星海面水溫圖」PNG -> 定量 SST 場 (numpy-only, self-calibrating)

關鍵設計
--------
1. 逐張偵測色階條 (colorbar) 位置與色彩序列, 建立 color -> 正規化位置 的查找表。
2. 逐張自動判定色階值域上限 (vmax): 圖上紅色等溫線一律繪於 5 的倍數,
   以此為內部真值檢定 30 / 35 兩種 CPT, 取殘差較小者。
   (本檔案庫 2018 年為 -2~35 之前的 -2~30 版, 2019 年起改為 -2~35,
    但圖面色標刻度直到 2020/21 才同步更新, 直接照色標判讀會低估 3~5 °C。)
3. 逐張以虛線經緯格線最小平方擬合像素<->經緯度轉換。
4. 以每格中位數聚合到 0.25° 網格, 對等溫線/地名字/海岸線等疊加物具強健性。
"""
import os, sys, json
import numpy as np
from PIL import Image

LON0, LON1, DLON = 116.0, 128.0, 0.25
LAT0, LAT1, DLAT = 18.0, 32.0, 0.25
NLON = int(round((LON1 - LON0) / DLON))   # 48
NLAT = int(round((LAT1 - LAT0) / DLAT))   # 56
GRID_LONS = [118, 120, 122, 124, 126]
GRID_LATS = [30, 28, 26, 24, 22, 20]
TOL = 26.0
VMIN = -2.0
VMAX_CANDIDATES = (30.0, 35.0)

_CACHE = {}


def _longest_run(mask):
    best = (0, 0, 0); s = None
    n = len(mask)
    for y in range(n):
        v = mask[y]
        if v and s is None: s = y
        if (not v) and s is not None:
            if y - s > best[0]: best = (y - s, s, y - 1)
            s = None
    if s is not None and n - s > best[0]: best = (n - s, s, n - 1)
    return best


def find_colorbar(a):
    h, w, _ = a.shape
    colored = (a.max(2) - a.min(2)) > 50
    xlo, xhi = int(w * 0.03), int(w * 0.22)
    info = {}
    for x in range(xlo, xhi):
        L, ya, yb = _longest_run(colored[:, x])
        if 0.20 * h < L < 0.55 * h and ya < 0.35 * h:
            info[x] = (L, ya, yb)
    if not info: return None
    xs = sorted(info)
    groups, cur = [], [xs[0]]
    for x in xs[1:]:
        if x - cur[-1] <= 2: cur.append(x)
        else: groups.append(cur); cur = [x]
    groups.append(cur)
    groups.sort(key=len, reverse=True)
    g = groups[0]
    if len(g) < 12: return None
    cx = g[len(g) // 2]
    _, y0, y1 = info[cx]
    return g[0], g[-1], y0, y1


def build_lut(a, cb):
    x0, x1, y0, y1 = cb
    strip = a[y0 + 2:y1 - 1, x0 + 3:x1 - 2].astype(np.float32)
    med = np.median(strip, axis=1)
    n = med.shape[0]
    pos = 1.0 - np.arange(n, dtype=np.float32) / (n - 1.0)
    return med, pos


def fit_axes(a, cb):
    dark = a.sum(2) < 250
    rf = dark.mean(1); cf = dark.mean(0)
    rows = np.where((rf > 0.2) & (rf < 0.8))[0]
    cols = np.where((cf > 0.2) & (cf < 0.8))[0]
    if cb:
        x0, x1, _, _ = cb
        cols = cols[(cols < x0 - 4) | (cols > x1 + 4)]
    cols = cols[cols > 150]

    def group(v):
        if len(v) == 0: return []
        g, cur = [], [v[0]]
        for x in v[1:]:
            if x - cur[-1] <= 3: cur.append(x)
            else: g.append(float(np.mean(cur))); cur = [x]
        g.append(float(np.mean(cur)))
        return g

    rows = group(rows); cols = group(cols)
    sx_nom, x116_nom = 100.167, 73.3
    sy_nom, y32_nom = 100.10, 47.8

    def refine(det, labels, nom_s, nom_o, is_lon):
        pred = {L: (nom_o + (L - LON0) * nom_s if is_lon else nom_o + (LAT1 - L) * nom_s)
                for L in labels}
        pairs = []
        for d in det:
            L = min(labels, key=lambda L: abs(pred[L] - d))
            if abs(pred[L] - d) < 8: pairs.append((L, d))
        if len(pairs) < 2: return nom_s, nom_o, len(pairs)
        X = np.array([p[0] for p in pairs], float)
        Y = np.array([p[1] for p in pairs], float)
        m, b = np.linalg.lstsq(np.vstack([X, np.ones_like(X)]).T, Y, rcond=None)[0]
        if is_lon: return abs(m), m * LON0 + b, len(pairs)
        return abs(m), m * LAT1 + b, len(pairs)

    sx, x116, nc = refine(cols, GRID_LONS, sx_nom, x116_nom, True)
    sy, y32, nr = refine(rows, GRID_LATS, sy_nom, y32_nom, False)
    return dict(sx=sx, x116=x116, sy=sy, y32=y32, nfit_lon=nc, nfit_lat=nr)


def _lookup(key, codes_u, lut, pos):
    c = _CACHE.setdefault(key, {'codes': np.zeros(0, np.int64), 'vals': np.zeros(0, np.float32)})
    n = len(c['codes'])
    out = np.full(len(codes_u), np.nan, np.float32)
    if n:
        idx = np.clip(np.searchsorted(c['codes'], codes_u), 0, n - 1)
        hit = c['codes'][idx] == codes_u
        out[hit] = c['vals'][idx[hit]]
    else:
        hit = np.zeros(len(codes_u), bool)
    miss = ~hit
    if miss.any():
        mc = codes_u[miss]
        rgb = np.stack([(mc >> 16) & 255, (mc >> 8) & 255, mc & 255], 1).astype(np.float32)
        vals = np.full(len(mc), np.nan, np.float32)
        for s in range(0, len(mc), 20000):
            blk = rgb[s:s + 20000]
            d = ((blk[:, None, :] - lut[None, :, :]) ** 2).sum(2)
            k = d.argmin(1)
            dist = np.sqrt(d[np.arange(len(k)), k])
            v = pos[k].copy()
            v[dist >= TOL] = np.nan
            vals[s:s + 20000] = v
        out[miss] = vals
        allc = np.concatenate([c['codes'], mc]); allv = np.concatenate([c['vals'], vals])
        o = np.argsort(allc); c['codes'] = allc[o]; c['vals'] = allv[o]
    return out


def _dilate(m, it):
    """numpy 版二值膨脹 (4-鄰域重複 it 次)。"""
    out = m
    for _ in range(it):
        o = out.copy()
        o[1:, :] |= out[:-1, :]; o[:-1, :] |= out[1:, :]
        o[:, 1:] |= out[:, :-1]; o[:, :-1] |= out[:, 1:]
        out = o
    return out


def pick_vmax(a, P, cb):
    """以紅色等溫線 (必為 5 的倍數) 判定色階上限。回傳 (vmax, score30, score35, n)。"""
    R, G, B = a[:, :, 0], a[:, :, 1], a[:, :, 2]
    red = (R > 225) & (G < 70) & (B < 45)
    red[:, :max(cb[1] + 20, 200)] = False
    ring = _dilate(red, 5) & (~_dilate(red, 2)) & np.isfinite(P)
    p = P[ring]
    if len(p) < 500:
        return None, None, None, len(p)
    sc = {}
    for vm in VMAX_CANDIDATES:
        v = VMIN + p * (vm - VMIN)
        sc[vm] = float(np.median(np.abs(v - np.round(v / 5.0) * 5.0)))
    best = min(sc, key=sc.get)
    return best, sc[30.0], sc[35.0], int(len(p))


def decode(path, vmax_hint=None):
    a = np.asarray(Image.open(path).convert('RGB')).astype(np.int16)
    h, w, _ = a.shape
    cb = find_colorbar(a)
    if cb is None: raise RuntimeError('colorbar not found')
    x0, x1, y0, y1 = cb
    cbh = y1 - y0
    lut, pos = build_lut(a, cb)
    ax = fit_axes(a, cb)

    sat = (a.max(2) - a.min(2)) > 45
    codes = (a[:, :, 0].astype(np.int64) << 16) | (a[:, :, 1].astype(np.int64) << 8) | a[:, :, 2].astype(np.int64)
    cu, inv = np.unique(codes[sat], return_inverse=True)
    p = _lookup(cbh, cu, lut, pos)
    P = np.full((h, w), np.nan, np.float32)
    P[sat] = p[inv]
    P[max(y0 - 4, 0):y1 + 5, max(x0 - 4, 0):x1 + 5] = np.nan

    vmax, s30, s35, nring = pick_vmax(a, P, cb)
    auto = vmax is not None
    if vmax is None:
        vmax = vmax_hint if vmax_hint else 35.0

    lons = LON0 + (np.arange(NLON) + 0.5) * DLON
    lats = LAT1 - (np.arange(NLAT) + 0.5) * DLAT
    px = ax['x116'] + (lons - LON0) * ax['sx']
    py = ax['y32'] + (LAT1 - lats) * ax['sy']
    hx = ax['sx'] * DLON / 2.0; hy = ax['sy'] * DLAT / 2.0

    out = np.full((NLAT, NLON), np.nan, np.float32)
    for j in range(NLAT):
        ya = max(int(round(py[j] - hy)), 0); yb = min(int(round(py[j] + hy)), h)
        if yb - ya < 3: continue
        rowP = P[ya:yb]
        for i in range(NLON):
            xa = max(int(round(px[i] - hx)), 0); xb = min(int(round(px[i] + hx)), w)
            if xb - xa < 3: continue
            blk = rowP[:, xa:xb]
            vv = blk[np.isfinite(blk)]
            if vv.size < 8: continue
            out[j, i] = np.median(vv)
    T = np.where(np.isfinite(out), VMIN + out * (vmax - VMIN), np.nan).astype(np.float32)
    nsat = int(np.nansum(T >= vmax - 0.2))
    if vmax >= 34.0:
        # -2~35 版色階不可能出現真實 35 度海水, 此類格點為紅色等溫線污染, 剔除
        T[T >= vmax - 0.2] = np.nan
    return T, dict(cbh=int(cbh), vmax=float(vmax), auto=bool(auto),
                   s30=s30, s35=s35, nring=int(nring), nsat=nsat,
                   nvalid=int(np.isfinite(T).sum()),
                   sx=round(float(ax['sx']), 3), sy=round(float(ax['sy']), 3),
                   x116=round(float(ax['x116']), 2), y32=round(float(ax['y32']), 2),
                   nfit=[int(ax['nfit_lon']), int(ax['nfit_lat'])])


if __name__ == '__main__':
    import time
    for p in sys.argv[1:]:
        t = time.time(); T, info = decode(p); v = T[np.isfinite(T)]
        print('%-44s n=%d min=%.2f max=%.2f mean=%.2f vmax=%.0f(auto=%s s30=%.3f s35=%.3f) sat=%d %.2fs'
              % (os.path.basename(p)[:44], len(v), v.min(), v.max(), v.mean(), info['vmax'],
                 info['auto'], info['s30'] or -1, info['s35'] or -1, info['nsat'], time.time() - t))
