#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
增補資料：將水試所網站「衛星海溫影像 → 衛星水溫影像圖」下載之 PNG
(https://www.tfrin.gov.tw/theme_list.php?theme=satellite&option=112，
 檔案位址 /files/satellite/<id>/A01_1.png，檔名不含日期) 反演並併入 out/cube.npz。

資料日期以 OCR 讀取圖面印製之 YYYY-MM-DD（需 tesseract）；網站公布日與資料日
常差 1~10 日（連假補發），不可用公布日代替。

用法:
    python3 tools/ingest_web.py <png 資料夾> [--overwrite]

已存在於立方中的日期預設不覆蓋（原始檔案庫優先），僅回報兩者差異作為一致性檢核。
之後執行 tools/build_web.py 重建儀表板資料。
"""
import os, sys, re, json, glob, subprocess, tempfile
import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sst_core as C

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
OUT = os.environ.get('SST_OUT', os.path.join(BASE, 'out'))
DATE = re.compile(r'(20\d\d)-(\d\d)-(\d\d)')


def read_date(path):
    """OCR 圖面左上方之資料日期（1350×1500 版面）。"""
    im = Image.open(path).convert('L')
    w, h = im.size
    box = (int(150 * w / 1350), int(420 * h / 1500), int(520 * w / 1350), int(540 * h / 1500))
    crop = im.crop(box).resize((740, 240))
    with tempfile.NamedTemporaryFile(suffix='.png') as f:
        crop.save(f.name)
        t = subprocess.run(['tesseract', f.name, '-', '--psm', '7',
                            '-c', 'tessedit_char_whitelist=0123456789-'],
                           capture_output=True, text=True).stdout
    m = DATE.search(t)
    return '%s-%s-%s' % m.groups() if m else None


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    src, overwrite = sys.argv[1], '--overwrite' in sys.argv
    cp = os.path.join(OUT, 'cube.npz')
    cube = dict(np.load(cp)) if os.path.exists(cp) else {}
    logp = os.path.join(OUT, 'log.jsonl')
    log = open(logp, 'a', encoding='utf-8')
    seen, nadd, nskip = set(), 0, 0
    for fp in sorted(glob.glob(os.path.join(src, '*.png'))):
        d = read_date(fp)
        if not d:
            print('  ?? date not readable:', os.path.basename(fp)); continue
        if d in seen:
            continue
        seen.add(d)
        try:
            T, info = C.decode(fp)
        except Exception as e:
            print('  !! %s %s' % (d, e)); continue
        if d in cube and not overwrite:
            a = cube[d]; m = np.isfinite(a) & np.isfinite(T)
            dd = T[m] - a[m]
            print('  = %s exists; web-vs-archive bias %+.3f rmse %.3f (n=%d)'
                  % (d, dd.mean(), np.sqrt((dd ** 2).mean()), m.sum()))
            nskip += 1; continue
        cube[d] = T.astype(np.float32)
        log.write(json.dumps(dict(date=d, ok=True, file=os.path.basename(fp), src='web', **info),
                             ensure_ascii=False) + '\n')
        if info['vmax'] < 34:
            print('  ** %s decoded with -2~30 palette — check for a production change!' % d)
        nadd += 1
    np.savez_compressed(cp, **cube)
    print('added %d, skipped %d existing, cube now %d days' % (nadd, nskip, len(cube)))


if __name__ == '__main__':
    main()
