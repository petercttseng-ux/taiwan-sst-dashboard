#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""衛星反演 SST × CTD/SADCP 15 弧分氣候圖集：綜合分析資料集產製。"""
import os, json, datetime as dt
import numpy as np, pandas as pd
import insitu_core as C

BASE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.environ.get('OUTDIR', os.path.join(BASE, 'site', 'data'))
CUBE = os.environ.get('CUBE', os.path.join(BASE, 'out', 'cube.npz'))
os.makedirs(OUT, exist_ok=True)

LON0, LAT1, D = 116.0, 32.0, 0.25
NLON, NLAT = 48, 56
LEGACY_END = dt.date(2019, 6, 2)
SEA = {13: [12,1,2], 14: [3,4,5], 15: [6,7,8], 16: [9,10,11],
       17: [10,11,12,1,2,3], 18: [4,5,6,7,8,9], 0: list(range(1,13))}
TPS = [0, 13, 14, 15, 16, 17, 18]
# 圖集網格（與衛星同解析度，但落在格線交點，故需內插）
GLON = np.arange(116.0, 128.001, 0.25)
GLAT = np.arange(18.0, 32.001, 0.25)

slon = LON0 + (np.arange(NLON) + 0.5) * D
slat = LAT1 - (np.arange(NLAT) + 0.5) * D


def bilin(F, lon, lat):
    x = (lon - slon[0]) / D; y = (slat[0] - lat) / D
    x0 = np.floor(x).astype(int); y0 = np.floor(y).astype(int)
    out = np.full(len(lon), np.nan)
    ok = (x0 >= 0) & (x0 < NLON - 1) & (y0 >= 0) & (y0 < NLAT - 1)
    fx = (x - x0)[ok]; fy = (y - y0)[ok]; a = x0[ok]; b = y0[ok]
    out[ok] = (F[b, a] * (1-fx) * (1-fy) + F[b, a+1] * fx * (1-fy)
               + F[b+1, a] * (1-fx) * fy + F[b+1, a+1] * fx * fy)
    return out


def grid_of(df, col):
    """把逐點診斷量攤成 GLAT×GLON 陣列。"""
    A = np.full((len(GLAT), len(GLON)), np.nan)
    ix = np.round((df.lon.values - GLON[0]) / 0.25).astype(int)
    iy = np.round((df.lat.values - GLAT[0]) / 0.25).astype(int)
    ok = (ix >= 0) & (ix < len(GLON)) & (iy >= 0) & (iy < len(GLAT))
    A[iy[ok], ix[ok]] = df[col].values[ok]
    return A


def enc(A, k=2):
    return [[None if not np.isfinite(v) else round(float(v), k) for v in row] for row in A]


def main():
    z = np.load(CUBE); dates = sorted(z.files)
    dob = [dt.date(*map(int, d.split('-'))) for d in dates]
    keep = [i for i, d in enumerate(dob) if d > LEGACY_END]
    cube = np.stack([z[dates[i]] for i in keep])
    mon = np.array([dob[i].month for i in keep])
    print('satellite days used: %d (%s -> %s)' % (len(keep), dates[keep[0]], dates[keep[-1]]))

    satF = {}
    for tp in TPS:
        with np.errstate(invalid='ignore'):
            satF[tp] = np.nanmean(cube[np.isin(mon, SEA[tp])], axis=0)

    ctd = C.load_ctd(); adcp = C.load_adcp()
    seasons, fields, scatter = [], {}, {}

    for tp in TPS:
        a = C.ctd_diagnostics(ctd, tp)
        a = a[(a.lon >= GLON[0]) & (a.lon <= GLON[-1])
              & (a.lat >= GLAT[0]) & (a.lat <= GLAT[-1])].copy()
        a['sat'] = bilin(satF[tp], a.lon.values, a.lat.values)
        a['bias'] = a.sat - a.t10
        a['swi'] = a.t10 - a.tm100          # 表層增溫指標：10 m 溫度 − 0-100 m 平均
        m = a.dropna(subset=['bias'])

        # 依混合層深度分組之偏差
        bins = [(0, 20, '<20'), (20, 30, '20–30'), (30, 45, '30–45'), (45, 1e4, '>45')]
        binned = []
        mm = m.dropna(subset=['mld'])
        for lo, hi, lab in bins:
            g = mm[(mm.mld > lo) & (mm.mld <= hi)]
            binned.append(dict(lab=lab, n=int(len(g)),
                               bias=None if len(g) < 3 else round(float(g.bias.mean()), 2),
                               dT=None if len(g) < 3 else round(float(g.dT.median()), 2)))

        d = m.bias.values
        seasons.append(dict(
            tp=tp, name=C.SEASON[tp], n=int(len(m)),
            bias=round(float(d.mean()), 3), rmse=round(float(np.sqrt((d**2).mean())), 3),
            median=round(float(np.median(d)), 3), sd=round(float(d.std(ddof=1)), 3),
            r=round(float(np.corrcoef(m.sat, m.t10)[0, 1]), 4),
            mld=None if m.mld.isna().all() else round(float(m.mld.median()), 1),
            dT=None if m.dT.isna().all() else round(float(m.dT.median()), 2),
            satdays=int(np.isin(mon, SEA[tp]).sum()), binned=binned))

        GX, GY = np.meshgrid(GLON, GLAT)
        satg = bilin(satF[tp], GX.ravel(), GY.ravel()).reshape(GY.shape)
        fields[tp] = dict(sat=enc(satg),
                          mld=enc(grid_of(a.dropna(subset=['mld']), 'mld'), 0),
                          dT=enc(grid_of(a.dropna(subset=['dT']), 'dT')),
                          t10=enc(grid_of(a.dropna(subset=['t10']), 't10')),
                          swi=enc(grid_of(a.dropna(subset=['swi']), 'swi')),
                          bias=enc(grid_of(m, 'bias')))
        if tp in (13, 15, 0):
            s = m.sample(min(len(m), 900), random_state=1)
            scatter[tp] = [[round(float(x), 2), round(float(y), 2)]
                           for x, y in zip(s.t10, s.sat)]
        print('tp=%2d n=%4d bias=%+.2f rmse=%.2f r=%.3f' % (tp, len(m), d.mean(),
              np.sqrt((d**2).mean()), np.corrcoef(m.sat, m.t10)[0, 1]))

    # ---------- 海流 ----------
    cur, axis = {}, {}
    for tp in [0, 13, 15, 14, 16]:
        g = C.adcp_diagnostics(adcp, tp, 100)
        g = g[(g.lon >= GLON[0]) & (g.lon <= GLON[-1])
              & (g.lat >= GLAT[0]) & (g.lat <= GLAT[-1])]
        cur[tp] = dict(u=enc(grid_of(g, 'u'), 3), v=enc(grid_of(g, 'v'), 3))
        k = g[(g.lon >= 121.0) & (g.lon <= 125.0) & (g.lat >= 21.5) & (g.lat <= 26.5)]
        ax = []
        for la, gg in k.groupby('lat'):
            gg = gg[gg.v > 0.15]
            if not len(gg): continue
            i = gg.v.idxmax()
            lo = float(gg.loc[i, 'lon'])
            sst = float(bilin(satF[tp], np.array([lo]), np.array([la]))[0])
            wst = float(bilin(satF[tp], np.array([lo - 1.0]), np.array([la]))[0])
            ax.append(dict(lat=float(la), lon=lo, v=round(float(gg.loc[i, 'v']), 3),
                           spd=round(float(gg.loc[i, 'spd']), 3),
                           sst=None if not np.isfinite(sst) else round(sst, 2),
                           sstW=None if not np.isfinite(wst) else round(wst, 2)))
        axis[tp] = ax

    meta = dict(generated=dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC'),
                glon=[float(x) for x in GLON], glat=[float(x) for x in GLAT],
                seasons=seasons, baseline=[dates[keep[0]], dates[keep[-1]]],
                nctd=int(len(ctd)), nadcp=int(len(adcp)))
    json.dump(dict(meta=meta, fields=fields, cur=cur, axis=axis, scatter=scatter),
              open(os.path.join(OUT, 'insitu.json'), 'w', encoding='utf-8'),
              ensure_ascii=False, separators=(',', ':'))
    print('insitu.json %.2f MB' % (os.path.getsize(os.path.join(OUT, 'insitu.json')) / 1e6))


if __name__ == '__main__':
    main()
