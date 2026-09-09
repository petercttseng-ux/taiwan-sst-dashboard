# -*- coding: utf-8 -*-
"""CTD/SADCP 15 弧分網格氣候圖集之剖面診斷量。"""
import os
import pandas as pd, numpy as np

U = os.environ.get('INSITU_DIR', '.') + '/'
SEASON = {0: '年平均', 13: '冬 DJF', 14: '春 MAM', 15: '夏 JJA',
          16: '秋 SON', 17: '東北季風期', 18: '西南季風期'}

def load_ctd():
    d = pd.read_csv(os.path.join(U, 'ctd_grid15moa', 'ctd_grid15moa.csv'), encoding='utf-8-sig')
    d.columns = ['lon','lat','p','t','tsd','s','ssd','sig','sigsd','tp','item']
    return d

def load_adcp():
    d = pd.read_csv(os.path.join(U, 'sadcp_grid15moa', 'sadcp_grid15moa.csv'), encoding='utf-8-sig')
    d.columns = ['lon','lat','z','u','usd','v','vsd','dir','spd','tp','item']
    return d

def ctd_diagnostics(ctd, tp):
    """逐格點：10 m 溫鹽、混合層深度、100 m 溫度、層化、0-100 m 熱含量代理。"""
    c = ctd[ctd.tp == tp].sort_values(['lon','lat','p'])
    recs = []
    for (lo, la), g in c.groupby(['lon','lat'], sort=False):
        p = g.p.values.astype(float); t = g.t.values.astype(float); s = g.s.values.astype(float)
        if len(p) < 4 or p.min() > 10 or p.max() < 60:
            continue
        t10 = float(np.interp(10, p, t)); s10 = float(np.interp(10, p, s))
        # 混合層深度：自 10 dbar 起向下，溫度低於 10 m 值 0.5 °C 之深度。
        # 僅在 10 dbar 以深搜尋，避免表層溫度逆轉（t(0) < t(10)）被誤判為 MLD = 0。
        below = np.where((t <= t10 - 0.5) & (p > 10))[0]
        if len(below):
            i = below[0]
            if i > 0:
                mld = p[i-1] + (t[i-1] - (t10 - 0.5)) / (t[i-1] - t[i]) * (p[i] - p[i-1])
            else:
                mld = p[i]
        else:
            mld = p.max() if p.max() >= 100 else np.nan   # 觀測深度內未達判準
        t100 = float(np.interp(100, p, t)) if p.max() >= 100 else np.nan
        sel = p <= 100
        tm100 = (float(np.trapezoid(t[sel], p[sel]) / (p[sel].max() - p[sel].min()))
                 if sel.sum() > 2 and p[sel].max() >= 80 else np.nan)
        recs.append(dict(lon=lo, lat=la, t10=t10, s10=s10, mld=float(mld),
                         t100=t100, dT=t10 - t100 if np.isfinite(t100) else np.nan,
                         tm100=tm100, nlev=len(p), pmax=float(p.max())))
    return pd.DataFrame(recs)

def adcp_diagnostics(adcp, tp, zmax=100):
    """逐格點：0-zmax m 深度平均流速（u,v,速度,流向）。"""
    a = adcp[(adcp.tp == tp) & (adcp.z <= zmax)]
    g = a.groupby(['lon','lat']).agg(u=('u','mean'), v=('v','mean'),
                                     n=('z','size'), zmax=('z','max')).reset_index()
    g = g[g.n >= 3]
    g['spd'] = np.hypot(g.u, g.v)
    g['dir'] = (np.degrees(np.arctan2(g.u, g.v)) + 360) % 360   # 流去向，自北順時針
    return g
