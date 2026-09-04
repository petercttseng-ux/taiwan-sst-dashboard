#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""由反演出的每日 SST 網格立方，產生儀表板所需的資料檔。"""
import os, json, datetime as dt
import numpy as np
from PIL import Image

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ.get('CUBE', os.path.join(BASE, '..', 'out', 'cube.npz'))
LOG = os.environ.get('QCLOG', os.path.join(BASE, '..', 'out', 'log.jsonl'))
OUTDIR = os.environ.get('OUTDIR', os.path.join(BASE, '..', 'data'))
os.makedirs(OUTDIR, exist_ok=True)

LON0, LAT1, D = 116.0, 32.0, 0.25
NLON, NLAT = 48, 56
NC = NLON * NLAT
SCALE, OFFSET, NODATA = 0.15, 0.0, 255
MIN_VALID = 700          # 有效格點過少者剔除
TMIN, TMAX = 5.0, 34.0   # 物理合理範圍
# 2019-06-03 圖檔填色調色盤由 -2~30 改為 -2~35，前後為兩個不同世代之產製流程，
# 域平均在該日出現約 +0.3~0.5 °C 之不連續。氣候基期僅採用其後之同質期間。
LEGACY_END = dt.date(2019, 6, 2)

REGIONS = [
    ("all",    "臺灣周邊全域",     "Full domain",          116.0, 128.0, 18.0, 32.0),
    ("nts",    "臺灣海峽北部",     "N. Taiwan Strait",     119.0, 121.0, 24.5, 26.0),
    ("sts",    "臺灣海峽南部",     "S. Taiwan Strait",     118.5, 120.5, 22.5, 24.5),
    ("penghu", "澎湖海域",         "Penghu",               119.0, 120.0, 23.0, 24.0),
    ("ne",     "東北部海域",       "NE off Taiwan",        121.5, 123.5, 25.0, 26.5),
    ("east",   "東部海域(黑潮)",   "E. Taiwan / Kuroshio",  121.5, 123.0, 22.5, 24.5),
    ("sw",     "西南部海域(高屏)", "SW off Taiwan",        119.5, 120.5, 21.5, 22.5),
    ("bashi",  "巴士海峽",         "Bashi Channel",        120.5, 122.0, 20.5, 22.0),
    ("nscs",   "南海北部(東沙)",   "N. South China Sea",   116.5, 118.5, 19.5, 21.5),
    ("secs",   "東海南部",         "S. East China Sea",    122.0, 126.0, 27.5, 30.5),
    ("nwp",    "西北太平洋",       "NW Pacific",           124.0, 127.0, 22.5, 25.5),
]

lons = LON0 + (np.arange(NLON) + 0.5) * D
lats = LAT1 - (np.arange(NLAT) + 0.5) * D


def region_mask(w, e, s, n):
    return np.outer((lats >= s) & (lats <= n), (lons >= w) & (lons <= e))


def doy_noleap(d):
    k = d.timetuple().tm_yday
    if d.year % 4 == 0 and (d.year % 100 != 0 or d.year % 400 == 0):
        if k == 60: return 59
        if k > 60: return k - 1
    return k


def circ_smooth2d(a, half):
    """沿 axis0 (365) 之圓形移動平均，忽略 NaN。a: (365, m)"""
    n = a.shape[0]
    ok = np.isfinite(a)
    v = np.where(ok, a, 0.0)
    s = np.zeros_like(v); c = np.zeros(a.shape, np.float32)
    for k in range(-half, half + 1):
        s += np.roll(v, k, axis=0); c += np.roll(ok.astype(np.float32), k, axis=0)
    return np.where(c > 0, s / np.maximum(c, 1), np.nan)


def circ_smooth1d(x, half):
    return circ_smooth2d(x.reshape(-1, 1).astype(np.float32), half).ravel()


def main():
    z = np.load(SRC)
    dates_all = sorted(z.files)
    keep, dropped = [], []
    nlowsd = 0
    for d in dates_all:
        T = z[d]
        T = np.where(np.isfinite(T) & (T >= TMIN) & (T <= TMAX), T, np.nan)
        if np.isfinite(T).sum() < MIN_VALID:
            dropped.append((d, 'too few valid cells')); continue
        sd = float(np.nanstd(T))
        if sd < 0.05:                     # 整片單色 = 原始圖檔產製失敗（無資料）
            dropped.append((d, 'uniform field, sd=%.3f' % sd)); continue
        if sd < 0.6:
            nlowsd += 1                   # 結構偏少（粗糙重繪版），保留但計數
        keep.append((d, T.astype(np.float32)))
    for d, why in dropped:
        print('  dropped %s: %s' % (d, why))
    dates = [k[0] for k in keep]
    cube = np.stack([k[1] for k in keep])
    nd = len(dates)
    print('days kept %d / %d (low-structure images kept: %d)   %s -> %s'
          % (nd, len(dates_all), nlowsd, dates[0], dates[-1]))

    dobj = [dt.date(*map(int, d.split('-'))) for d in dates]
    doys = np.array([doy_noleap(d) for d in dobj])
    legacy = np.array([d <= LEGACY_END for d in dobj])
    nleg = int(legacy.sum())
    print('legacy segment (-2~30 palette): %d days, baseline uses the other %d' % (nleg, nd - nleg))

    # ---------- 逐格 日序氣候值 (±7 日窗，僅用同質期間) ----------
    flatc = cube.reshape(nd, -1)
    clim = np.full((365, NC), np.nan, np.float32)
    for k in range(365):
        dd = np.abs(doys - (k + 1)); dd = np.minimum(dd, 365 - dd)
        sel = (dd <= 7) & (~legacy)
        if sel.sum() >= 3:
            with np.errstate(invalid='ignore'):
                clim[k] = np.nanmean(flatc[sel], axis=0)
    clim = circ_smooth2d(clim, 5)
    print('climatology done (homogeneous baseline)')

    # ---------- 輸出 PNG 立方 ----------
    def to_png(arr2, path):
        q = np.full(arr2.shape, NODATA, np.uint8)
        ok = np.isfinite(arr2)
        q[ok] = np.clip(np.round((arr2[ok] - OFFSET) / SCALE), 0, 254).astype(np.uint8)
        Image.fromarray(q, mode='L').save(path, optimize=True)
        return os.path.getsize(path)

    s1 = to_png(flatc, os.path.join(OUTDIR, 'sst_cube.png'))
    s2 = to_png(clim, os.path.join(OUTDIR, 'clim_cube.png'))
    print('sst_cube %.2f MB   clim_cube %.2f MB' % (s1 / 1e6, s2 / 1e6))

    # 以量化後之值重建，確保前端顯示與統計一致
    qs = np.round(flatc / SCALE); qs = np.where(np.isfinite(qs) & (qs >= 0) & (qs <= 254), qs, np.nan) * SCALE
    qc_ = np.round(clim / SCALE); qc_ = np.where(np.isfinite(qc_) & (qc_ >= 0) & (qc_ <= 254), qc_, np.nan) * SCALE

    # ---------- 分區統計 ----------
    series = {}
    for rid, zh, en, w, e, s, n in REGIONS:
        m = region_mask(w, e, s, n).ravel()
        with np.errstate(invalid='ignore'):
            mean = np.nanmean(np.where(m, qs, np.nan), axis=1)
            cmean = np.nanmean(np.where(m, qc_, np.nan), axis=1)
        an = mean - cmean[doys - 1]

        p90 = np.full(365, np.nan)
        for k in range(365):
            dd = np.abs(doys - (k + 1)); dd = np.minimum(dd, 365 - dd)
            sel = (dd <= 7) & np.isfinite(mean) & (~legacy)
            if sel.sum() >= 10:
                p90[k] = np.nanpercentile(mean[sel], 90)
        p90 = circ_smooth1d(p90, 5)
        thr = p90[doys - 1]

        # 熱浪偵測僅在同質期間進行（前段色階世代不同，不可比）
        over = np.isfinite(mean) & np.isfinite(thr) & (mean > thr) & (~legacy)
        events, i = [], 0
        while i < nd:
            if over[i]:
                jx = i
                while jx + 1 < nd and over[jx + 1] and (dobj[jx + 1] - dobj[jx]).days <= 2:
                    jx += 1
                dur = (dobj[jx] - dobj[i]).days + 1
                if dur >= 5:
                    sl = slice(i, jx + 1)
                    events.append(dict(start=dates[i], end=dates[jx], days=int(dur),
                                       imax=round(float(np.nanmax(mean[sl] - thr[sl])), 2),
                                       imean=round(float(np.nanmean(mean[sl] - thr[sl])), 2),
                                       tmax=round(float(np.nanmax(mean[sl])), 2)))
                i = jx + 1
            else:
                i += 1

        R = lambda a, k=2: [None if not np.isfinite(v) else round(float(v), k) for v in a]
        series[rid] = dict(mean=R(mean), anom=R(an), clim=R(cmean), p90=R(p90), mhw=events)
        print('%-7s n=%4d mean=%5.2f  MHW %d 事件 / %d 日'
              % (rid, int(np.isfinite(mean).sum()), np.nanmean(mean),
                 len(events), sum(e['days'] for e in events)))

    # ---------- 年×月 距平矩陣 ----------
    years = sorted({d.year for d in dobj})
    ymi = np.array([years.index(d.year) for d in dobj])
    mo = np.array([d.month - 1 for d in dobj])
    monthly = {}
    for rid, *_ in REGIONS:
        arr = np.array([np.nan if v is None else v for v in series[rid]['anom']])
        Msum = np.zeros((len(years), 12)); Cnt = np.zeros((len(years), 12), int)
        ok = np.isfinite(arr)
        np.add.at(Msum, (ymi[ok], mo[ok]), arr[ok])
        np.add.at(Cnt, (ymi[ok], mo[ok]), 1)
        Mv = np.where(Cnt > 0, Msum / np.maximum(Cnt, 1), np.nan)
        monthly[rid] = dict(
            v=[[None if not np.isfinite(x) else round(float(x), 2) for x in r] for r in Mv],
            n=[[int(x) for x in r] for r in Cnt])

    # 不連續量測：以同質氣候值為基準，前後段之平均距平差
    with np.errstate(invalid='ignore'):
        amean = np.nanmean(qs - qc_[doys - 1], axis=1)
    jump = float(np.nanmean(amean[legacy]) - np.nanmean(amean[~legacy]))

    meta = dict(
        generated=dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC'),
        grid=dict(lon0=LON0, lat1=LAT1, d=D, nlon=NLON, nlat=NLAT),
        scale=SCALE, offset=OFFSET, nodata=NODATA,
        dates=dates, years=years,
        legacy_n=nleg, legacy_end=LEGACY_END.isoformat(),
        baseline=[dates[nleg], dates[-1]],
        legacy_offset=round(jump, 2),
        regions=[dict(id=r[0], zh=r[1], en=r[2], box=[r[3], r[4], r[5], r[6]]) for r in REGIONS],
    )
    print('legacy segment offset vs homogeneous baseline: %+.2f °C' % jump)

    # ---------- QC 統計 ----------
    qc = None
    if os.path.exists(LOG):
        rows = {}
        for line in open(LOG, encoding='utf-8'):
            try: r = json.loads(line)
            except Exception: continue
            rows[r['date']] = r
        agg = {}
        for d in dates:
            r = rows.get(d)
            if not r: continue
            y = int(d[:4]); a = agg.setdefault(y, dict(y=y, n=0, ok=0, v30=0, v35=0, auto=0, nvalid=0, tmean=0.0, tn=0))
            a['n'] += 1
            if r.get('ok'):
                a['ok'] += 1
                if r.get('auto'): a['auto'] += 1
                if float(r.get('vmax', 35)) < 32: a['v30'] += 1
                else: a['v35'] += 1
                a['nvalid'] += r.get('nvalid', 0)
        for y, a in agg.items():
            k = max(a['ok'], 1); a['nvalid'] /= k
        for y, a in agg.items():
            sel = [i for i in range(nd) if dobj[i].year == y]
            a['tmean'] = float(np.nanmean(qs[sel])) if sel else float('nan')
            a.pop('tn', None)
        yl = [agg[y] for y in sorted(agg)]
        tot = dict(n=sum(a['n'] for a in yl), ok=sum(a['ok'] for a in yl),
                   v30=sum(a['v30'] for a in yl), v35=sum(a['v35'] for a in yl),
                   auto=sum(a['auto'] for a in yl),
                   nvalid=float(np.mean([a['nvalid'] for a in yl])) if yl else 0.0,
                   tmean=float(np.nanmean(qs)))
        qc = dict(years=yl, total=tot, dropped=len(dates_all) - nd,
                  dropped_list=[{'date': a, 'why': b} for a, b in dropped],
                  lowsd=nlowsd, min_valid=MIN_VALID, trange=[TMIN, TMAX],
                  legacy_n=nleg, legacy_end=LEGACY_END.isoformat(),
                  legacy_offset=round(jump, 2))

    def dump(name, obj):
        p = os.path.join(OUTDIR, name)
        json.dump(obj, open(p, 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
        print('%-14s %.2f MB' % (name, os.path.getsize(p) / 1e6))

    dump('meta.json', meta); dump('series.json', series); dump('monthly.json', monthly)
    if qc: dump('qc.json', qc)


if __name__ == '__main__':
    main()
