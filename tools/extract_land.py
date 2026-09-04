#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""由原始水溫圖萃取全解析度陸地遮罩。

作法：取跨年代的多張影像，各自以其經緯格線擬合的轉換投影到共同的 0.01 度網格，
逐點統計「非海色」（低飽和度＝灰底陸地／黑色文字與線條）出現的比例。
等溫線逐日移動，陸地不動，因此高比例者即為陸地；門檻設 0.85。
"""
import os, sys, glob, random
import numpy as np
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sst_core as C

LON0, LON1 = 116.0, 128.0
LAT0, LAT1 = 18.0, 32.0
RES = 0.01
NX = int(round((LON1 - LON0) / RES))   # 1200
NY = int(round((LAT1 - LAT0) / RES))   # 1400
THRESH = 0.85

ROOT = '/mnt/user-data/uploads/衛星水溫圖'


def sample_files(n=28):
    pats = []
    for y in ['2018', '2019', '2020', '2021', '2022', '2023']:
        pats += sorted(glob.glob(os.path.join(ROOT, y, '*.png')))
    pats += sorted(glob.glob(os.path.join(ROOT, 'SST 2025*', '*domain-taiwan*contour.png')))
    pats += sorted(glob.glob(os.path.join(ROOT, 'SST2026*.png')))
    pats = [p for p in pats if 'contour_only' not in p.lower() and 'onlly' not in p.lower()]
    random.seed(7)
    return random.sample(pats, min(n, len(pats)))


def main():
    lons = LON0 + (np.arange(NX) + 0.5) * RES
    lats = LAT1 - (np.arange(NY) + 0.5) * RES
    acc = np.zeros((NY, NX), np.float32)
    used = 0
    for p in sample_files():
        try:
            a = np.asarray(Image.open(p).convert('RGB')).astype(np.int16)
            cb = C.find_colorbar(a)
            if cb is None:
                continue
            ax = C.fit_axes(a, cb)
            h, w, _ = a.shape
            px = np.round(ax['x116'] + (lons - LON0) * ax['sx']).astype(int)
            py = np.round(ax['y32'] + (LAT1 - lats) * ax['sy']).astype(int)
            if px.min() < 0 or px.max() >= w or py.min() < 0 or py.max() >= h:
                continue
            blk = a[np.ix_(py, px)]                      # (NY, NX, 3)
            nonocean = (blk.max(2) - blk.min(2)) <= 45   # 灰／黑＝非海色
            acc += nonocean
            used += 1
        except Exception as e:
            print('skip', os.path.basename(p), e)
    print('images used:', used)
    frac = acc / max(used, 1)
    land = frac >= THRESH

    # 形態學清理：先去除海上孤立雜點，再補回陸地內部小孔
    def dil(m, it):
        o = m
        for _ in range(it):
            q = o.copy()
            q[1:, :] |= o[:-1, :]; q[:-1, :] |= o[1:, :]
            q[:, 1:] |= o[:, :-1]; q[:, :-1] |= o[:, 1:]
            o = q
        return o
    def ero(m, it):
        return ~dil(~m, it)
    land = dil(ero(land, 1), 1)     # opening：去雜點
    land = ero(dil(land, 2), 2)     # closing：補孔（等溫線壓在陸地邊緣造成的缺口）

    # 圖框黑線會被判為陸地，去掉最外圈
    b = 4
    land[:b, :] = False; land[-b:, :] = False
    land[:, :b] = False; land[:, -b:] = False

    # 色階條疊在陸地上，會在陸地中挖出一塊「海」。凡是無法從圖緣連通到的
    # 非陸地區域皆屬此類封閉孔洞，一律填為陸地。
    from collections import deque
    sea = ~land
    reach = np.zeros_like(sea)
    q = deque()
    for x in range(NX):
        for y in (0, NY - 1):
            if sea[y, x] and not reach[y, x]: reach[y, x] = True; q.append((y, x))
    for y in range(NY):
        for x in (0, NX - 1):
            if sea[y, x] and not reach[y, x]: reach[y, x] = True; q.append((y, x))
    while q:
        y, x = q.popleft()
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ny_, nx_ = y + dy, x + dx
            if 0 <= ny_ < NY and 0 <= nx_ < NX and sea[ny_, nx_] and not reach[ny_, nx_]:
                reach[ny_, nx_] = True; q.append((ny_, nx_))
    holes = sea & ~reach
    print('filled enclosed holes: %d px' % holes.sum())
    land |= holes

    print('land fraction: %.1f %%' % (100 * land.mean()))
    out = np.where(land, 255, 0).astype(np.uint8)
    im = Image.fromarray(out, mode='L').convert('1')
    for d in ['site/data', 'repo/data']:
        os.makedirs(d, exist_ok=True)
        im.save(os.path.join(d, 'land_mask.png'), optimize=True, bits=1)
    print('land_mask.png %.1f KB' % (os.path.getsize('site/data/land_mask.png') / 1e3))


if __name__ == '__main__':
    main()
