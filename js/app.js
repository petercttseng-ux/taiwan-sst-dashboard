/* 臺灣周邊海域衛星遙測海面水溫時空分布儀表板
   Taiwan Seas Satellite SST Spatiotemporal Dashboard */
'use strict';

var M = null, SST = null, CLIM = null, SER = null, MON = null, QC = null;
var ND = 0, NC = 0, NLON = 0, NLAT = 0, SCALE = 0.15, NODATA = 255;
var NLEG = 0;   // 前段 (-2~30 色階世代) 日數，索引 0..NLEG-1
var LANDCV = null;   // 由原始圖檔萃取之全解析度陸地圖層
var DOY = null, YEAR = null, MONTHI = null;
var DAYIDX = null, SPAN = 1;   // 真實時間軸：距首日之日數

/* ---------- palettes ---------- */
function clamp(x, a, b) { return x < a ? a : x > b ? b : x; }
function lerp(a, b, t) { return a + (b - a) * t; }

// 海溫色階：深藍→青→綠→黃→橙→紅 (thermal-like, perceptually ordered)
var SST_STOPS = [
  [0.00, 12, 24, 92], [0.12, 20, 82, 178], [0.26, 32, 150, 200],
  [0.40, 60, 190, 170], [0.53, 150, 214, 110], [0.66, 238, 220, 80],
  [0.79, 246, 160, 48], [0.90, 226, 84, 40], [1.00, 150, 24, 32]
];
// 距平色階：藍-白-紅 diverging
var ANO_STOPS = [
  [0.00, 20, 60, 150], [0.22, 60, 130, 205], [0.42, 165, 205, 235],
  [0.50, 244, 244, 240], [0.58, 246, 196, 170], [0.78, 226, 105, 66],
  [1.00, 140, 18, 26]
];
function ramp(stops, t) {
  t = clamp(t, 0, 1);
  for (var i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      var a = stops[i - 1], b = stops[i];
      var u = (t - a[0]) / (b[0] - a[0] || 1);
      return [lerp(a[1], b[1], u) | 0, lerp(a[2], b[2], u) | 0, lerp(a[3], b[3], u) | 0];
    }
  }
  var l = stops[stops.length - 1]; return [l[1], l[2], l[3]];
}
function sstColor(v, lo, hi) { return ramp(SST_STOPS, (v - lo) / (hi - lo)); }
function anoColor(v, r) { return ramp(ANO_STOPS, 0.5 + v / (2 * r)); }
function css(c) { return 'rgb(' + c[0] + ',' + c[1] + ',' + c[2] + ')'; }

/* ---------- data access ---------- */
function sstAt(d, c) { var p = SST[d * NC + c]; return p === NODATA ? NaN : p * SCALE; }
function climAt(k, c) { var p = CLIM[k * NC + c]; return p === NODATA ? NaN : p * SCALE; }
function anomAt(d, c) { var a = sstAt(d, c), b = climAt(DOY[d] - 1, c); return a - b; }

function fmt(v, n) { return (v === null || v === undefined || !isFinite(v)) ? '—' : v.toFixed(n === undefined ? 2 : n); }
function sgn(v, n) { return !isFinite(v) ? '—' : (v >= 0 ? '+' : '') + v.toFixed(n === undefined ? 2 : n); }

function meanOf(a) { var s = 0, n = 0; for (var i = 0; i < a.length; i++) if (isFinite(a[i]) && a[i] !== null) { s += a[i]; n++; } return n ? s / n : NaN; }
function movavg(a, w) {
  var out = new Array(a.length), h = w >> 1;
  for (var i = 0; i < a.length; i++) {
    var s = 0, n = 0;
    for (var k = -h; k <= h; k++) { var j = i + k; if (j >= 0 && j < a.length && a[j] !== null && isFinite(a[j])) { s += a[j]; n++; } }
    out[i] = n >= w * 0.4 ? s / n : null;
  }
  return out;
}
/* 對「實際經過日數」做最小平方迴歸；off 為 y[0] 在完整序列中的索引。
   斜率單位 = °C/日，故乘 3652.5 即為 °C/10 年。 */
function linfit(y, off) {
  off = off || 0;
  var sx = 0, sy = 0, sxx = 0, sxy = 0, n = 0;
  for (var i = 0; i < y.length; i++) {
    if (y[i] === null || !isFinite(y[i])) continue;
    var x = DAYIDX[i + off] - DAYIDX[off];
    sx += x; sy += y[i]; sxx += x * x; sxy += x * y[i]; n++;
  }
  if (n < 10) return null;
  var d = n * sxx - sx * sx;
  if (!d) return null;
  var m = (n * sxy - sx * sy) / d;
  return { m: m, b: (sy - m * sx) / n, n: n };
}
function spanDays(i0, i1) { return DAYIDX[i1] - DAYIDX[i0]; }


/* 前段（色階世代不同）遮罩與標註 */
function shadeLegacy(g, x0, y0, w, h) {
  if (!NLEG) return;
  var xe = x0 + (tOf(NLEG - 1) + tOf(NLEG)) / 2 * w;
  g.save();
  g.fillStyle = 'rgba(255,179,0,.09)';
  g.fillRect(x0, y0, xe - x0, h);
  g.strokeStyle = 'rgba(255,179,0,.55)'; g.setLineDash([4, 3]); g.lineWidth = 1;
  g.beginPath(); g.moveTo(xe, y0); g.lineTo(xe, y0 + h); g.stroke();
  g.setLineDash([]);
  g.fillStyle = 'rgba(255,179,0,.85)'; g.font = '10px sans-serif';
  g.textAlign = 'left'; g.textBaseline = 'top';
  g.fillText('舊色階世代', x0 + 4, y0 + 4);
  g.restore();
}
function homoSlice(a) { return a.slice(NLEG); }
/* 取「最後 nDays 個日曆日」內的樣本（資料有缺漏，用筆數會跨越空窗） */
function lastDays(a, nDays) {
  var cut = DAYIDX[ND - 1] - nDays, out = [];
  for (var i = 0; i < ND; i++) if (DAYIDX[i] > cut) out.push(a[i]);
  return out;
}
function tOf(i) { return DAYIDX[i] / SPAN; }
/* 於資料缺口 (> 5 日) 斷開折線，避免把空窗畫成直線 */
function gapBreak(pts) {
  var out = [];
  for (var i = 0; i < pts.length; i++) {
    if (i > 0 && DAYIDX[i] - DAYIDX[i - 1] > 5) out.push(null);
    out.push(pts[i]);
  }
  return out;
}

/* 年份刻度：格線畫在年初，標籤置於該年資料重心，樣本過少或過於接近者不標 */
function yearTicks(minGap, minN) {
  minGap = minGap || 0.045; minN = minN || 3;
  var raw = [];
  for (var k = 0; k < M.years.length; k++) {
    var y = M.years[k], i0 = -1, i1 = -1, n = 0;
    for (var i = 0; i < ND; i++) if (YEAR[i] === y) { if (i0 < 0) i0 = i; i1 = i; n++; }
    if (i0 < 0) continue;
    raw.push({ g: tOf(i0), t: (tOf(i0) + tOf(i1)) / 2, l: String(y), n: n });
  }
  var last = -9;
  for (var j = 0; j < raw.length; j++) {
    if (raw[j].n < minN || raw[j].t - last < minGap) { raw[j].l = ''; }
    else last = raw[j].t;
  }
  return raw;
}

/* ---------- canvas helpers ---------- */
function prep(cv, h) {
  var dpr = Math.min(window.devicePixelRatio || 1, 2);
  var w = cv.clientWidth || cv.parentNode.clientWidth || 800;
  if (h) cv.style.height = h + 'px';
  var hh = h || cv.clientHeight || 300;
  cv.width = Math.max(1, Math.round(w * dpr));
  cv.height = Math.max(1, Math.round(hh * dpr));
  var g = cv.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, hh);
  return { g: g, w: w, h: hh };
}
function axes(g, x0, y0, w, h, xt, yt, opt) {
  opt = opt || {};
  g.strokeStyle = '#26405a'; g.lineWidth = 1;
  g.fillStyle = '#8ba6bd'; g.font = '11px ' + (opt.mono ? 'monospace' : 'sans-serif');
  g.textAlign = 'right'; g.textBaseline = 'middle';
  for (var i = 0; i < yt.length; i++) {
    var y = Math.round(y0 + h - yt[i].t * h) + 0.5;
    g.beginPath(); g.moveTo(x0, y); g.lineTo(x0 + w, y);
    g.strokeStyle = yt[i].zero ? '#40607e' : '#1e3448'; g.stroke();
    g.fillText(yt[i].l, x0 - 6, y);
  }
  g.textAlign = 'center'; g.textBaseline = 'top';
  for (var j = 0; j < xt.length; j++) {
    var x = Math.round(x0 + xt[j].t * w) + 0.5;
    var xg = Math.round(x0 + (xt[j].g !== undefined ? xt[j].g : xt[j].t) * w) + 0.5;
    if (xt[j].grid !== false) { g.beginPath(); g.moveTo(xg, y0); g.lineTo(xg, y0 + h); g.strokeStyle = '#1a2d3f'; g.stroke(); }
    if (xt[j].l) g.fillText(xt[j].l, x, y0 + h + 6);
  }
  g.strokeStyle = '#26405a';
  g.strokeRect(x0 + 0.5, y0 + 0.5, w, h);
}
function niceTicks(lo, hi, n) {
  var span = hi - lo; if (!(span > 0)) return [{ t: 0.5, l: fmt(lo, 1) }];
  var step = Math.pow(10, Math.floor(Math.log10(span / n)));
  var err = span / n / step;
  if (err >= 7.5) step *= 10; else if (err >= 3.5) step *= 5; else if (err >= 1.5) step *= 2;
  var out = [], s = Math.ceil(lo / step) * step;
  for (var v = s; v <= hi + 1e-9; v += step) {
    var dec = step < 0.05 ? 2 : step < 0.5 ? 1 : step < 1 ? 1 : 0;
    out.push({ t: (v - lo) / span, l: v.toFixed(dec), zero: Math.abs(v) < 1e-9 });
  }
  return out;
}
function polyline(g, pts, color, wd, dash) {
  g.save(); g.strokeStyle = color; g.lineWidth = wd || 1.5; g.lineJoin = 'round'; g.lineCap = 'round';
  if (dash) g.setLineDash(dash);
  g.beginPath(); var pen = false;
  for (var i = 0; i < pts.length; i++) {
    var p = pts[i];
    if (p === null) { pen = false; continue; }
    if (!pen) { g.moveTo(p[0], p[1]); pen = true; } else g.lineTo(p[0], p[1]);
  }
  g.stroke(); g.restore();
}

/* ---------- 陸地圖層 ----------
   land_mask.png 為 0.01° 網格 (1200×1400) 之二值遮罩，直接由原始水溫圖萃取，
   與資料共用同一組經緯度轉換，故海岸線與網格完全對齊。
   於載入時一次性上色為「陸地填色 + 海岸線」圖層，之後每次繪圖只需 drawImage。 */
var LAND_FILL = [52, 64, 78], LAND_COAST = [150, 176, 198];
function buildLand(mask, w, h) {
  var im = new ImageData(w, h), d = im.data;
  for (var y = 0; y < h; y++) {
    for (var x = 0; x < w; x++) {
      var i = y * w + x;
      if (mask[i] < 128) continue;                 // 海：透明
      // 海岸線＝與海相鄰的陸地像素
      var edge = (x === 0 || mask[i - 1] < 128) || (x === w - 1 || mask[i + 1] < 128) ||
                 (y === 0 || mask[i - w] < 128) || (y === h - 1 || mask[i + w] < 128);
      var c = edge ? LAND_COAST : LAND_FILL, o = i * 4;
      d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2]; d[o + 3] = 255;
    }
  }
  var cv = document.createElement('canvas');
  cv.width = w; cv.height = h;
  cv.getContext('2d').putImageData(im, 0, 0);
  return cv;
}

/* ---------- map rendering ---------- */
function fieldImage(getter, colf) {
  var im = new ImageData(NLON, NLAT);
  for (var j = 0; j < NLAT; j++) for (var i = 0; i < NLON; i++) {
    var c = j * NLON + i, v = getter(c), o = c * 4;
    if (!isFinite(v)) { im.data[o + 3] = 0; continue; }
    var col = colf(v);
    im.data[o] = col[0]; im.data[o + 1] = col[1]; im.data[o + 2] = col[2]; im.data[o + 3] = 255;
  }
  return im;
}
var _off = null;
function drawField(cv, getter, colf, opt) {
  opt = opt || {};
  var dpr = Math.min(window.devicePixelRatio || 1, 2);
  // 等距圓柱投影，於中緯度套用 cos(φ) 修正，使島形比例接近實際
  var ASPECT = (NLAT / NLON) / Math.cos(25 * Math.PI / 180);
  var W = cv.parentNode.clientWidth || cv.clientWidth || 640;
  var H = W * ASPECT;
  var maxH = opt.maxH || 640;
  if (H > maxH) { H = maxH; W = H / ASPECT; }
  W = Math.round(W); H = Math.round(H);
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
  cv.style.width = W + 'px'; cv.style.height = H + 'px'; cv.style.margin = '0 auto';
  var g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.fillStyle = '#0a1219'; g.fillRect(0, 0, W, H);

  if (!_off) { _off = document.createElement('canvas'); _off.width = NLON; _off.height = NLAT; }
  _off.getContext('2d').putImageData(fieldImage(getter, colf), 0, 0);
  g.imageSmoothingEnabled = opt.smooth !== false;
  g.imageSmoothingQuality = 'high';
  g.drawImage(_off, 0, 0, W, H);

  // 陸地疊在水溫場之上：遮住平滑造成的越岸暈染，並給出清晰海岸線
  if (LANDCV) {
    g.imageSmoothingEnabled = true;
    g.imageSmoothingQuality = 'high';
    g.drawImage(LANDCV, 0, 0, W, H);
  }

  // graticule
  g.save();
  g.strokeStyle = 'rgba(255,255,255,.14)'; g.lineWidth = 1;
  g.font = '10px monospace';
  g.shadowColor = 'rgba(0,0,0,.85)'; g.shadowBlur = 3;
  g.fillStyle = 'rgba(255,255,255,.9)';
  var gl = M.grid;
  for (var lon = 118; lon <= 126; lon += 2) {
    var x = (lon - gl.lon0) / (gl.nlon * gl.d) * W;
    g.beginPath(); g.moveTo(x, 0); g.lineTo(x, H); g.stroke();
    g.textAlign = 'left'; g.textBaseline = 'bottom'; g.fillText(lon + '°E', x + 3, H - 3);
  }
  for (var lat = 20; lat <= 30; lat += 2) {
    var y = (gl.lat1 - lat) / (gl.nlat * gl.d) * H;
    g.beginPath(); g.moveTo(0, y); g.lineTo(W, y); g.stroke();
    g.textAlign = 'left'; g.textBaseline = 'bottom'; g.fillText(lat + '°N', 3, y - 2);
  }
  g.restore();

  if (opt.iso) drawIso(g, getter, W, H, opt.isoStep || 1, opt.isoLabel);
  g.strokeStyle = 'rgba(139,166,189,.45)'; g.lineWidth = 1;
  g.strokeRect(0.5, 0.5, W - 1, H - 1);
  return { g: g, W: W, H: H };
}
/* marching squares 等溫線 */
function drawIso(g, getter, W, H, step, labelMode) {
  var vals = new Float32Array(NC);
  for (var c = 0; c < NC; c++) vals[c] = getter(c);
  var lo = Infinity, hi = -Infinity;
  for (var k = 0; k < NC; k++) if (isFinite(vals[k])) { if (vals[k] < lo) lo = vals[k]; if (vals[k] > hi) hi = vals[k]; }
  if (!isFinite(lo)) return;
  var cw = W / NLON, ch = H / NLAT;
  g.save(); g.lineWidth = 1; g.font = '10px monospace'; g.textAlign = 'center'; g.textBaseline = 'middle';
  for (var L = Math.ceil(lo / step) * step; L <= hi; L += step) {
    var major = Math.abs(L / (step * 5) - Math.round(L / (step * 5))) < 1e-6;
    g.strokeStyle = major ? 'rgba(255,255,255,.75)' : 'rgba(255,255,255,.3)';
    g.lineWidth = major ? 1.3 : 0.8;
    g.beginPath();
    var lbl = null;
    for (var j = 0; j < NLAT - 1; j++) for (var i = 0; i < NLON - 1; i++) {
      var a = vals[j * NLON + i], b = vals[j * NLON + i + 1],
        d = vals[(j + 1) * NLON + i], e = vals[(j + 1) * NLON + i + 1];
      if (!(isFinite(a) && isFinite(b) && isFinite(d) && isFinite(e))) continue;
      var idx = (a > L ? 8 : 0) | (b > L ? 4 : 0) | (e > L ? 2 : 0) | (d > L ? 1 : 0);
      if (idx === 0 || idx === 15) continue;
      var x0 = (i + 0.5) * cw, y0 = (j + 0.5) * ch;
      function T(p, q) { return (L - p) / (q - p); }
      var top = [x0 + T(a, b) * cw, y0], rgt = [x0 + cw, y0 + T(b, e) * ch],
        bot = [x0 + T(d, e) * cw, y0 + ch], lft = [x0, y0 + T(a, d) * ch];
      var segs = { 1: [lft, bot], 2: [bot, rgt], 3: [lft, rgt], 4: [top, rgt], 5: [top, lft, bot, rgt], 6: [top, bot], 7: [top, lft], 8: [top, lft], 9: [top, bot], 10: [top, rgt, bot, lft], 11: [top, rgt], 12: [lft, rgt], 13: [bot, rgt], 14: [lft, bot] }[idx];
      for (var s = 0; s < segs.length; s += 2) { g.moveTo(segs[s][0], segs[s][1]); g.lineTo(segs[s + 1][0], segs[s + 1][1]); }
      if (major && !lbl && i > 4 && i < NLON - 6) lbl = top;
    }
    g.stroke();
    if (labelMode !== false && lbl) {
      g.fillStyle = 'rgba(6,20,32,.85)';
      var t = String(Math.round(L * 10) / 10);
      var wd = g.measureText(t).width + 6;
      g.fillRect(lbl[0] - wd / 2, lbl[1] - 7, wd, 13);
      g.fillStyle = '#fff'; g.fillText(t, lbl[0], lbl[1]);
    }
  }
  g.restore();
}
function drawLegendBar(cv, colf) {
  var g = cv.getContext('2d');
  for (var x = 0; x < cv.width; x++) {
    g.fillStyle = css(colf(x / (cv.width - 1)));
    g.fillRect(x, 0, 1, cv.height);
  }
}

/* =================== BOOT =================== */
function loadPNG(url) {
  return new Promise(function (res, rej) {
    var im = new Image();
    im.onload = function () {
      var c = document.createElement('canvas'); c.width = im.width; c.height = im.height;
      var g = c.getContext('2d', { willReadFrequently: true });
      g.drawImage(im, 0, 0);
      var d = g.getImageData(0, 0, im.width, im.height).data;
      var out = new Uint8Array(im.width * im.height);
      for (var i = 0, n = out.length; i < n; i++) out[i] = d[i * 4];
      res({ data: out, w: im.width, h: im.height });
    };
    im.onerror = function () { rej(new Error('load ' + url)); };
    im.src = url;
  });
}
function j(u) { return fetch(u).then(function (r) { return r.json(); }); }

Promise.all([j('data/meta.json'), j('data/series.json'), j('data/monthly.json'),
  loadPNG('data/sst_cube.png'), loadPNG('data/clim_cube.png'),
  fetch('data/qc.json').then(function (r) { return r.ok ? r.json() : null; }).catch(function () { return null; }),
  loadPNG('data/land_mask.png').catch(function () { return null; }),
  fetch('data/insitu.json').then(function (r) { return r.ok ? r.json() : null; }).catch(function () { return null; })
]).then(function (r) {
  M = r[0]; SER = r[1]; MON = r[2]; SST = r[3].data; CLIM = r[4].data; QC = r[5];
  if (r[6]) LANDCV = buildLand(r[6].data, r[6].w, r[6].h);
  IS = r[7];
  NLON = M.grid.nlon; NLAT = M.grid.nlat; NC = NLON * NLAT;
  SCALE = M.scale; NODATA = M.nodata; ND = M.dates.length; NLEG = M.legacy_n || 0;
  DOY = new Int16Array(ND); YEAR = new Int16Array(ND); MONTHI = new Int8Array(ND);
  for (var i = 0; i < ND; i++) {
    var p = M.dates[i].split('-'), y = +p[0], mo = +p[1], dd = +p[2];
    YEAR[i] = y; MONTHI[i] = mo;
    var leap = (y % 4 === 0 && (y % 100 !== 0 || y % 400 === 0));
    var cum = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    var k = cum[mo - 1] + dd + (leap && mo > 2 ? 1 : 0);
    if (leap) { if (k === 60) k = 59; else if (k > 60) k -= 1; }
    DOY[i] = k;
  }
  // 真實時間軸（資料有缺漏，索引軸會把 2024 年的空窗壓成零寬度）
  DAYIDX = new Int32Array(ND);
  var t0 = Date.parse(M.dates[0] + 'T00:00:00Z');
  for (var q = 0; q < ND; q++)
    DAYIDX[q] = Math.round((Date.parse(M.dates[q] + 'T00:00:00Z') - t0) / 86400000);
  SPAN = Math.max(DAYIDX[ND - 1], 1);
  init();
  document.getElementById('loading').style.display = 'none';
}).catch(function (e) {
  document.getElementById('loading').innerHTML = '<div style="color:#ff7043">資料載入失敗：' + e.message + '</div>';
});

/* =================== INIT =================== */
var cur = 0, mapMode = 'sst', playing = null, pickCell = -1;

function init() {
  document.getElementById('hdr-range').textContent = M.dates[0] + ' ～ ' + M.dates[ND - 1];
  document.getElementById('hdr-n').textContent = ND.toLocaleString();
  document.getElementById('genTime').textContent = M.generated;
  var sp = document.getElementById('specPeriod');
  if (sp) sp.textContent = M.dates[0] + ' ～ ' + M.dates[ND - 1] + '，共 ' + ND.toLocaleString() +
    ' 日（同質基期 ' + (ND - NLEG).toLocaleString() + ' 日）';

  // tabs
  var btns = document.querySelectorAll('nav.tabs button');
  for (var i = 0; i < btns.length; i++) btns[i].onclick = function () {
    for (var k = 0; k < btns.length; k++) btns[k].classList.remove('on');
    this.classList.add('on');
    var ps = document.querySelectorAll('section.page');
    for (var q = 0; q < ps.length; q++) ps[q].classList.remove('on');
    document.getElementById('p-' + this.dataset.p).classList.add('on');
    render(this.dataset.p);
  };

  // region selects
  ['srRegion', 'mxRegion', 'mhRegion'].forEach(function (id) {
    var s = document.getElementById(id);
    M.regions.forEach(function (rg) {
      var o = document.createElement('option'); o.value = rg.id;
      o.textContent = rg.zh + '　' + rg.en; s.appendChild(o);
    });
    s.value = 'all';
  });
  document.getElementById('srRegion').value = 'east';
  document.getElementById('mhRegion').value = 'all';
  document.getElementById('srRegion').onchange = drawSeries;
  document.getElementById('mxRegion').onchange = drawMatrix;
  document.getElementById('mhRegion').onchange = drawMHW;
  segBind('srMode', drawSeries);
  segBind('hovMode', drawHov);
  segBind('mapMode', function (m) { mapMode = m; drawMap(); });

  // map controls
  var sl = document.getElementById('dslider');
  sl.max = ND - 1; sl.value = ND - 1; cur = ND - 1;
  sl.oninput = function () { cur = +this.value; drawMap(); };
  document.getElementById('btnPrev').onclick = function () { step(-1); };
  document.getElementById('btnNext').onclick = function () { step(1); };
  document.getElementById('btnPlay').onclick = togglePlay;
  document.getElementById('btnIso').onclick = function () { this.classList.toggle('on'); drawMap(); };
  document.getElementById('btnSmooth').onclick = function () { this.classList.toggle('on'); drawMap(); };
  var jy = document.getElementById('jumpYear');
  jy.appendChild(new Option('跳至年份…', ''));
  M.years.forEach(function (y) { jy.appendChild(new Option(y, y)); });
  jy.onchange = function () {
    if (!this.value) return;
    for (var i = 0; i < ND; i++) if (YEAR[i] === +this.value) { cur = i; break; }
    document.getElementById('dslider').value = cur; drawMap();
  };
  document.getElementById('mapcv').onmousemove = mapHover;
  document.getElementById('mapcv').onmouseleave = function () { document.getElementById('tip').style.display = 'none'; };
  document.getElementById('mapcv').onclick = mapClick;

  // hov lon bands
  var hl = document.getElementById('hovLon');
  [[116, 128, '全域 116–128°E'], [118, 121, '臺灣海峽 118–121°E'], [121, 124, '臺灣以東 121–124°E'], [124, 128, '西北太平洋 124–128°E']]
    .forEach(function (b) { hl.appendChild(new Option(b[2], b[0] + ',' + b[1])); });
  hl.onchange = drawHov;

  drawLegendBar(document.getElementById('mxLg'), function (t) { return ramp(ANO_STOPS, t); });

  if (IS) {
    var seaOpts = IS.meta.seasons.filter(function (s) { return s.tp !== 17 && s.tp !== 18; });
    [['vaMapSea', '15'], ['vtSea', '15']].forEach(function (cfg) {
      var el = document.getElementById(cfg[0]);
      seaOpts.forEach(function (s) { el.appendChild(new Option(s.name, s.tp)); });
      el.value = cfg[1];
    });
    var cu = document.getElementById('cuSea');
    [0, 13, 14, 15, 16].forEach(function (tp) {
      var nm = IS.meta.seasons.filter(function (s) { return s.tp === tp; })[0].name;
      cu.appendChild(new Option(nm, tp));
    });
    cu.value = '0';
    document.getElementById('vaMapSea').onchange = drawBiasMap;
    document.getElementById('vtSea').onchange = drawVert;
    cu.onchange = drawCurrents;
    document.getElementById('cuAxis').onclick = function () { this.classList.toggle('on'); drawCurrents(); };
    document.getElementById('cuSST').onclick = function () { this.classList.toggle('on'); drawCurrents(); };
    segBind('vaScMode', drawScatter);
    segBind('vtField', drawVert);
  } else {
    ['valid', 'vert'].forEach(function (p) {
      var b = document.querySelector('nav.tabs button[data-p="' + p + '"]');
      if (b) b.style.display = 'none';
    });
  }
  window.addEventListener('resize', debounce(function () {
    var on = document.querySelector('nav.tabs button.on');
    render(on.dataset.p);
  }, 250));

  render('overview');
}
function segBind(id, cb) {
  var el = document.getElementById(id), bs = el.querySelectorAll('button');
  for (var i = 0; i < bs.length; i++) bs[i].onclick = function () {
    for (var k = 0; k < bs.length; k++) bs[k].classList.remove('on');
    this.classList.add('on'); cb(this.dataset.m);
  };
}
function segVal(id) { return document.getElementById(id).querySelector('button.on').dataset.m; }
function debounce(f, ms) { var t; return function () { clearTimeout(t); t = setTimeout(f, ms); }; }

var drawn = {};
function render(p) {
  if (p === 'overview') drawOverview();
  else if (p === 'map') drawMap();
  else if (p === 'series') drawSeries();
  else if (p === 'hov') drawHov();
  else if (p === 'matrix') drawMatrix();
  else if (p === 'mhw') drawMHW();
  else if (p === 'valid') drawValid();
  else if (p === 'vert') { drawVert(); drawCurrents(); }
  else if (p === 'method') drawQC();
}

/* =================== OVERVIEW =================== */
function drawOverview() {
  var all = SER.all, mean = all.mean, anom = all.anom;
  var last = ND - 1;
  // recent 30d
  var r30 = lastDays(mean, 30), a30 = lastDays(anom, 30);
  var a365 = lastDays(anom, 365);
  var n30 = a30.filter(function (v) { return v !== null; }).length;
  var n365 = a365.filter(function (v) { return v !== null; }).length;
  var fit = linfit(homoSlice(anom), NLEG);     // 距平序列，已去除季節循環
  var fitRaw = linfit(homoSlice(mean), NLEG);  // 供繪圖用（原始尺度）
  var trendDecade = fit ? fit.m * 3652.5 : NaN;

  var nMHW = 0, dMHW = 0, maxEv = null;
  all.mhw.forEach(function (e) { nMHW++; dMHW += e.days; if (!maxEv || e.days > maxEv.days) maxEv = e; });

  var hottest = -Infinity, hotIdx = -1;
  for (var i = 0; i < ND; i++) if (mean[i] !== null && mean[i] > hottest) { hottest = mean[i]; hotIdx = i; }

  var k = [
    ['最新觀測日 Latest', M.dates[last], fmt(mean[last]) + ' °C　距平 ' + sgn(anom[last]), 'acc'],
    ['近 30 日均溫', fmt(meanOf(r30)) + ' °C', '距平 ' + sgn(meanOf(a30)) + '（' + n30 + ' 日樣本）', meanOf(a30) >= 0 ? 'warm' : 'cool'],
    ['近 1 年距平', sgn(meanOf(a365)) + ' °C', '相對同質基期氣候值（' + n365 + ' 日樣本）', meanOf(a365) >= 0 ? 'warm' : 'cool'],
    ['期間距平趨勢', sgn(trendDecade) + ' °C/10yr', '僅 ' + (fit ? fit.n : 0) + ' 日同質資料，非氣候趨勢', trendDecade >= 0 ? 'warm' : 'cool'],
    ['歷史最暖日', fmt(hottest) + ' °C', hotIdx >= 0 ? M.dates[hotIdx] : '—', 'warm'],
    ['海洋熱浪事件', nMHW + ' 次', '累計 ' + dMHW + ' 日' + (maxEv ? '，最長 ' + maxEv.days + ' 日' : ''), 'warm'],
    ['資料涵蓋', ND.toLocaleString() + ' 日', M.dates[0] + ' 起', 'acc'],
    ['空間解析度', '0.25°', NLON + ' × ' + NLAT + ' 格　116–128°E / 18–32°N', 'acc']
  ];
  document.getElementById('kpis').innerHTML = k.map(function (x) {
    return '<div class="kpi ' + x[3] + '"><div class="lab">' + x[0] + '</div><div class="val">' + x[1] + '</div><div class="sub">' + x[2] + '</div></div>';
  }).join('');

  // trend chart
  var P = prep(document.getElementById('ovTrend'), 300);
  var x0 = 46, y0 = 12, w = P.w - 66, h = P.h - 42;
  var lo = Infinity, hi = -Infinity;
  for (var q = 0; q < ND; q++) if (mean[q] !== null) { lo = Math.min(lo, mean[q]); hi = Math.max(hi, mean[q]); }
  var pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
  var xt = yearTicks();
  axes(P.g, x0, y0, w, h, xt, niceTicks(lo, hi, 6));
  shadeLegacy(P.g, x0, y0, w, h);
  function pt(i, v) { return v === null || !isFinite(v) ? null : [x0 + tOf(i) * w, y0 + h - (v - lo) / (hi - lo) * h]; }
  polyline(P.g, gapBreak(mean.map(function (v, i) { return pt(i, v); })), 'rgba(139,166,189,.45)', 1);
  polyline(P.g, gapBreak(movavg(mean, 31).map(function (v, i) { return pt(i, v); })), '#4fc3f7', 2);
  if (fitRaw) polyline(P.g, [pt(NLEG, fitRaw.b), pt(ND - 1, fitRaw.b + fitRaw.m * spanDays(NLEG, ND - 1))], '#ff7043', 1.8, [6, 4]);
  P.g.fillStyle = '#8ba6bd'; P.g.font = '11px sans-serif'; P.g.textAlign = 'left';
  P.g.fillText('°C', 6, y0 + 4);
  P.g.fillStyle = '#ff7043'; P.g.textAlign = 'right';
  P.g.fillText('距平趨勢 ' + sgn(trendDecade) + ' °C / 10yr', x0 + w - 6, y0 + 12);

  // region table
  var yrs = M.years;
  var html = '<thead><tr><th>海域 Region</th>' + yrs.map(function (y) { return '<th>' + y + '</th>'; }).join('') + '<th>全期均溫</th></tr></thead><tbody>';
  M.regions.forEach(function (rg) {
    var s = SER[rg.id]; if (!s) return;
    html += '<tr><td>' + rg.zh + '</td>';
    yrs.forEach(function (y) {
      var acc = 0, n = 0;
      for (var i = 0; i < ND; i++) if (YEAR[i] === y && s.anom[i] !== null) { acc += s.anom[i]; n++; }
      if (n < 3) { html += '<td class="mono" style="color:#4b6478">—</td>'; return; }
      var v = acc / n, thin = n < 30;
      html += '<td class="mono" title="' + n + ' 日樣本' + (thin ? '（樣本少，僅供參考）' : '') + '" style="background:' +
        css(anoColor(v, 1.2)) + ';color:' + (Math.abs(v) > 0.7 ? '#fff' : '#0d1620') +
        (thin ? ';opacity:.45;font-style:italic' : '') + '">' + sgn(v, 2) + (thin ? '*' : '') + '</td>';
    });
    html += '<td class="mono">' + fmt(meanOf(s.mean)) + '</td></tr>';
  });
  document.getElementById('ovRegion').innerHTML = html + '</tbody>';

  // latest field
  document.getElementById('ovMapDate').textContent = M.dates[ND - 1] + '　實測 SST (°C)';
  var rng = fieldRange(ND - 1);
  drawField(document.getElementById('ovMap'), function (c) { return sstAt(ND - 1, c); },
    function (v) { return sstColor(v, rng[0], rng[1]); }, { smooth: true, iso: true, isoStep: 1, maxH: 380 });
}

/* =================== MAP =================== */
function fieldRange(d) {
  var lo = Infinity, hi = -Infinity;
  for (var c = 0; c < NC; c++) { var v = sstAt(d, c); if (isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } }
  if (!isFinite(lo)) return [20, 30];
  // 低對比日（如盛夏）值域很窄，直接鋪滿整條色階會把最冷格點壓到最深的藍，
  // 看起來像無資料；兩端各留 0.5 °C 緩衝。
  var pad = Math.max((hi - lo) * 0.08, 0.5);
  return [Math.floor(lo - pad), Math.ceil(hi + pad)];
}
function step(k) {
  cur = clamp(cur + k, 0, ND - 1);
  document.getElementById('dslider').value = cur; drawMap();
}
function togglePlay() {
  var b = document.getElementById('btnPlay');
  if (playing) { clearInterval(playing); playing = null; b.textContent = '▶ 播放'; b.classList.remove('on'); return; }
  b.textContent = '❚❚ 暫停'; b.classList.add('on');
  var ms = +document.getElementById('spd').value;
  playing = setInterval(function () {
    cur = (cur + 1) % ND;
    document.getElementById('dslider').value = cur;
    drawMap();
  }, ms);
}
var mapGeom = null;
function drawMap() {
  var iso = document.getElementById('btnIso').classList.contains('on');
  var sm = document.getElementById('btnSmooth').classList.contains('on');
  document.getElementById('mapDate').textContent = M.dates[cur];
  var getter, colf, lgf, lo, hi;
  if (mapMode === 'sst') {
    var r = fieldRange(cur); lo = r[0]; hi = r[1];
    getter = function (c) { return sstAt(cur, c); };
    colf = function (v) { return sstColor(v, lo, hi); };
    lgf = function (t) { return ramp(SST_STOPS, t); };
    document.getElementById('lgMin').textContent = lo; document.getElementById('lgMax').textContent = hi;
  } else {
    lo = -3; hi = 3;
    getter = function (c) { return anomAt(cur, c); };
    colf = function (v) { return anoColor(v, 3); };
    lgf = function (t) { return ramp(ANO_STOPS, t); };
    document.getElementById('lgMin').textContent = '-3'; document.getElementById('lgMax').textContent = '+3';
  }
  mapGeom = drawField(document.getElementById('mapcv'), getter, colf,
    { smooth: sm, iso: iso, isoStep: mapMode === 'sst' ? 1 : 1 });
  drawLegendBar(document.getElementById('lgcv'), lgf);

  // daily stats
  var vals = [], anos = [];
  for (var c = 0; c < NC; c++) { var v = sstAt(cur, c); if (isFinite(v)) vals.push(v); var a = anomAt(cur, c); if (isFinite(a)) anos.push(a); }
  vals.sort(function (a, b) { return a - b; });
  var st = [
    ['日期 Date', M.dates[cur] + (cur < NLEG ? ' <span class="pill" style="color:#ffb300;border-color:#6b4d12">舊色階</span>' : '')],
    ['有效格點 Valid', vals.length + ' / ' + NC],
    ['全域平均 Mean', fmt(meanOf(vals)) + ' °C'],
    ['最低 Min', fmt(vals[0]) + ' °C'],
    ['中位數 Median', fmt(vals[vals.length >> 1]) + ' °C'],
    ['最高 Max', fmt(vals[vals.length - 1]) + ' °C'],
    ['平均距平 Anomaly', sgn(meanOf(anos)) + ' °C'],
    ['距平 &gt; +1 °C 面積', (100 * anos.filter(function (x) { return x > 1; }).length / Math.max(anos.length, 1)).toFixed(1) + ' %']
  ];
  document.getElementById('mapStats').innerHTML = st.map(function (x) {
    return '<tr><td>' + x[0] + '</td><td class="mono">' + x[1] + '</td></tr>';
  }).join('');
  if (pickCell >= 0) drawPoint(pickCell);
}
function cellFromEvent(ev) {
  var cv = document.getElementById('mapcv'), r = cv.getBoundingClientRect();
  var i = Math.floor((ev.clientX - r.left) / r.width * NLON);
  var jj = Math.floor((ev.clientY - r.top) / r.height * NLAT);
  if (i < 0 || i >= NLON || jj < 0 || jj >= NLAT) return -1;
  return jj * NLON + i;
}
function cellLL(c) {
  var i = c % NLON, jj = (c / NLON) | 0;
  return [M.grid.lon0 + (i + 0.5) * M.grid.d, M.grid.lat1 - (jj + 0.5) * M.grid.d];
}
function mapHover(ev) {
  var c = cellFromEvent(ev), tip = document.getElementById('tip');
  if (c < 0) { tip.style.display = 'none'; return; }
  var v = sstAt(cur, c), a = anomAt(cur, c), ll = cellLL(c);
  if (!isFinite(v)) { tip.style.display = 'none'; return; }
  tip.style.display = 'block';
  tip.innerHTML = ll[0].toFixed(2) + '°E ' + ll[1].toFixed(2) + '°N<br>' + fmt(v) + ' °C　(' + sgn(a) + ')';
  var w = document.getElementById('mapwrap').getBoundingClientRect();
  tip.style.left = Math.min(ev.clientX - w.left + 12, w.width - 160) + 'px';
  tip.style.top = (ev.clientY - w.top + 12) + 'px';
}
function mapClick(ev) { var c = cellFromEvent(ev); if (c >= 0 && isFinite(sstAt(cur, c))) { pickCell = c; drawPoint(c); } }
function drawPoint(c) {
  var ll = cellLL(c);
  document.getElementById('ptHint').innerHTML = '<b class="mono">' + ll[0].toFixed(2) + '°E, ' + ll[1].toFixed(2) + '°N</b>　逐日海溫（灰）與 31 日平滑（藍）';
  var ser = new Array(ND);
  for (var i = 0; i < ND; i++) { var v = sstAt(i, c); ser[i] = isFinite(v) ? v : null; }
  var P = prep(document.getElementById('ptCv'), 220);
  var x0 = 34, y0 = 10, w = P.w - 44, h = P.h - 34;
  var lo = Infinity, hi = -Infinity;
  for (var q = 0; q < ND; q++) if (ser[q] !== null) { lo = Math.min(lo, ser[q]); hi = Math.max(hi, ser[q]); }
  if (!isFinite(lo)) return;
  var xt = yearTicks(0.09).map(function (t) { return { t: t.t, g: t.g, l: t.l ? t.l.slice(2) : '' }; });
  axes(P.g, x0, y0, w, h, xt, niceTicks(lo, hi, 4));
  function pt(i, v) { return v === null ? null : [x0 + tOf(i) * w, y0 + h - (v - lo) / (hi - lo) * h]; }
  polyline(P.g, gapBreak(ser.map(function (v, i) { return pt(i, v); })), 'rgba(139,166,189,.45)', 1);
  polyline(P.g, gapBreak(movavg(ser, 31).map(function (v, i) { return pt(i, v); })), '#4fc3f7', 1.8);
  var m = pt(cur, ser[cur]);
  if (m) { P.g.fillStyle = '#ff7043'; P.g.beginPath(); P.g.arc(m[0], m[1], 3.5, 0, 7); P.g.fill(); }
}

/* =================== SERIES =================== */
var YCOL = ['#4fc3f7', '#66bb6a', '#ffb300', '#ff7043', '#ab47bc', '#26c6da', '#ec407a', '#8d6e63', '#9ccc65', '#5c6bc0'];
function drawSeries() {
  var rid = document.getElementById('srRegion').value, mode = segVal('srMode');
  var rg = M.regions.filter(function (r) { return r.id === rid; })[0], s = SER[rid];
  document.getElementById('srTitle').textContent = rg.zh + '　' + rg.en + '　(' + rg.box[0] + '–' + rg.box[1] + '°E, ' + rg.box[2] + '–' + rg.box[3] + '°N)';
  var P = prep(document.getElementById('srCv'), 380);
  var x0 = 48, y0 = 14, w = P.w - 68, h = P.h - 44;
  var g = P.g, leg = '';

  if (mode === 'year') {
    document.getElementById('srNote').textContent = '以日序（day-of-year）疊圖比較各年季節循環；白色粗線為同質基期之氣候平均。';
    var lo = Infinity, hi = -Infinity;
    for (var i = 0; i < ND; i++) if (s.mean[i] !== null) { lo = Math.min(lo, s.mean[i]); hi = Math.max(hi, s.mean[i]); }
    var mLbl = ['1月', '2月', '3月', '4月', '5月', '6月', '7月', '8月', '9月', '10月', '11月', '12月'];
    var cum = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    var xt = cum.map(function (c, i) { return { t: c / 365, l: mLbl[i] }; });
    axes(g, x0, y0, w, h, xt, niceTicks(lo, hi, 6));
    M.years.forEach(function (y, yi) {
      var pts = [];
      for (var i = 0; i < ND; i++) if (YEAR[i] === y) {
        pts.push(s.mean[i] === null ? null : [x0 + (DOY[i] - 1) / 364 * w, y0 + h - (s.mean[i] - lo) / (hi - lo) * h]);
      }
      polyline(g, pts, YCOL[yi % YCOL.length], 1.4);
      leg += '<span style="color:' + YCOL[yi % YCOL.length] + '">■</span> ' + y + '　';
    });
    var cp = [];
    for (var k = 0; k < 365; k++) cp.push(s.clim[k] === null ? null : [x0 + k / 364 * w, y0 + h - (s.clim[k] - lo) / (hi - lo) * h]);
    polyline(g, cp, 'rgba(255,255,255,.9)', 2.6);
    leg += '<span style="color:#fff">■</span> 氣候平均 Climatology';
  } else if (mode === 'full') {
    document.getElementById('srNote').textContent = '逐日海域平均海溫（灰）與 31 日移動平均（藍），紅虛線為線性趨勢。';
    var lo2 = Infinity, hi2 = -Infinity;
    for (var q = 0; q < ND; q++) if (s.mean[q] !== null) { lo2 = Math.min(lo2, s.mean[q]); hi2 = Math.max(hi2, s.mean[q]); }
    var xt2 = yearTicks();
    axes(g, x0, y0, w, h, xt2, niceTicks(lo2, hi2, 6));
    shadeLegacy(g, x0, y0, w, h);
    var pf = function (i, v) { return v === null ? null : [x0 + tOf(i) * w, y0 + h - (v - lo2) / (hi2 - lo2) * h]; };
    polyline(g, gapBreak(s.mean.map(function (v, i) { return pf(i, v); })), 'rgba(139,166,189,.45)', 1);
    polyline(g, gapBreak(movavg(s.mean, 31).map(function (v, i) { return pf(i, v); })), '#4fc3f7', 2);
    var f = linfit(homoSlice(s.anom), NLEG), fr = linfit(homoSlice(s.mean), NLEG);
    if (f && fr) {
      polyline(g, [pf(NLEG, fr.b), pf(ND - 1, fr.b + fr.m * spanDays(NLEG, ND - 1))], '#ff7043', 1.8, [6, 4]);
      leg = '線性趨勢 <b style="color:#ff7043">' + sgn(f.m * 3652.5) + ' °C / 10 年</b>（僅以 ' + M.baseline[0] + ' 起之同質資料迴歸；記錄僅約 6 年，不足以代表氣候趨勢）';
    }
  } else {
    document.getElementById('srNote').textContent = '逐日距平（相對本期日序氣候值），紅＝偏暖、藍＝偏冷；黑線為 31 日移動平均。';
    var r = 0;
    for (var z = 0; z < ND; z++) if (s.anom[z] !== null) r = Math.max(r, Math.abs(s.anom[z]));
    r = Math.ceil(r * 2) / 2;
    var xt3 = yearTicks();
    axes(g, x0, y0, w, h, xt3, niceTicks(-r, r, 6));
    shadeLegacy(g, x0, y0, w, h);
    var zy = y0 + h / 2;
    for (var i2 = 0; i2 < ND; i2++) {
      var v = s.anom[i2]; if (v === null) continue;
      var x = x0 + tOf(i2) * w, yy = y0 + h - (v + r) / (2 * r) * h;
      g.strokeStyle = v >= 0 ? 'rgba(255,112,67,.75)' : 'rgba(66,165,245,.75)';
      g.lineWidth = Math.max(w / ND, 0.8);
      g.beginPath(); g.moveTo(x, zy); g.lineTo(x, yy); g.stroke();
    }
    polyline(g, gapBreak(movavg(s.anom, 31).map(function (v, i) { return v === null ? null : [x0 + tOf(i) * w, y0 + h - (v + r) / (2 * r) * h]; })), '#ffffff', 1.8);
    leg = '<span style="color:#ff7043">■</span> 正距平　<span style="color:#42a5f5">■</span> 負距平　<span style="color:#fff">■</span> 31 日平均';
  }
  g.fillStyle = '#8ba6bd'; g.font = '11px sans-serif'; g.textAlign = 'left'; g.fillText('°C', 6, y0 + 4);
  document.getElementById('srLegend').innerHTML = leg;

  // monthly table
  var mLab = ['1月', '2月', '3月', '4月', '5月', '6月', '7月', '8月', '9月', '10月', '11月', '12月'];
  var rows = [['月份 Month'], ['平均 Mean'], ['最低 Min'], ['最高 Max'], ['樣本 N']];
  for (var mo = 1; mo <= 12; mo++) {
    var vv = [];
    for (var i3 = 0; i3 < ND; i3++) if (MONTHI[i3] === mo && s.mean[i3] !== null) vv.push(s.mean[i3]);
    rows[0].push(mLab[mo - 1]);
    rows[1].push(vv.length ? fmt(meanOf(vv)) : '—');
    rows[2].push(vv.length ? fmt(Math.min.apply(null, vv)) : '—');
    rows[3].push(vv.length ? fmt(Math.max.apply(null, vv)) : '—');
    rows[4].push(vv.length);
  }
  var t = '<thead><tr>' + rows[0].map(function (c, i) { return '<th>' + (i ? c : '　') + '</th>'; }).join('') + '</tr></thead><tbody>';
  for (var k2 = 1; k2 < rows.length; k2++)
    t += '<tr><td>' + rows[k2][0] + '</td>' + rows[k2].slice(1).map(function (c) { return '<td class="mono">' + c + '</td>'; }).join('') + '</tr>';
  document.getElementById('srMonthly').innerHTML = t + '</tbody>';
}

/* =================== HOVMOLLER =================== */
function drawHov() {
  var mode = segVal('hovMode');
  var band = document.getElementById('hovLon').value.split(',').map(Number);
  var i0 = Math.max(0, Math.round((band[0] - M.grid.lon0) / M.grid.d));
  var i1 = Math.min(NLON, Math.round((band[1] - M.grid.lon0) / M.grid.d));
  hovRender(document.getElementById('hovCv'), 420, mode, i0, i1, 'lat');
  hovRender(document.getElementById('hovCv2'), 340, mode, 0, NLON, 'lon');
}
function hovRender(cv, hgt, mode, i0, i1, axis) {
  var n = axis === 'lat' ? NLAT : NLON;
  var arr = new Float32Array(ND * n); arr.fill(NaN);
  for (var d = 0; d < ND; d++) {
    for (var k = 0; k < n; k++) {
      var s = 0, c = 0;
      if (axis === 'lat') {
        for (var i = i0; i < i1; i++) { var v = mode === 'sst' ? sstAt(d, k * NLON + i) : anomAt(d, k * NLON + i); if (isFinite(v)) { s += v; c++; } }
      } else {
        for (var jj = 0; jj < NLAT; jj++) { var v2 = mode === 'sst' ? sstAt(d, jj * NLON + k) : anomAt(d, jj * NLON + k); if (isFinite(v2)) { s += v2; c++; } }
      }
      if (c >= 3) arr[d * n + k] = s / c;
    }
  }
  var lo = Infinity, hi = -Infinity;
  for (var q = 0; q < arr.length; q++) if (isFinite(arr[q])) { if (arr[q] < lo) lo = arr[q]; if (arr[q] > hi) hi = arr[q]; }
  if (mode === 'anom') { var r = Math.max(Math.abs(lo), Math.abs(hi)); lo = -r; hi = r; }
  var colf = mode === 'sst' ? function (v) { return sstColor(v, lo, hi); } : function (v) { return anoColor(v, hi); };

  // 以真實日曆日為橫軸；無資料之日期留白，使 2024 年的空窗如實呈現
  var WCOL = SPAN + 1;
  var im = new ImageData(WCOL, n);
  for (var q = 0; q < im.data.length; q += 4) {
    im.data[q] = 26; im.data[q + 1] = 38; im.data[q + 2] = 50; im.data[q + 3] = 255;
  }
  for (var d2 = 0; d2 < ND; d2++) {
    var xc = DAYIDX[d2];
    for (var k2 = 0; k2 < n; k2++) {
      var v3 = arr[d2 * n + k2], o = (k2 * WCOL + xc) * 4;
      if (!isFinite(v3)) { im.data[o] = 12; im.data[o + 1] = 22; im.data[o + 2] = 31; im.data[o + 3] = 255; continue; }
      var col = colf(v3);
      im.data[o] = col[0]; im.data[o + 1] = col[1]; im.data[o + 2] = col[2]; im.data[o + 3] = 255;
    }
  }
  var off = document.createElement('canvas'); off.width = WCOL; off.height = n;
  off.getContext('2d').putImageData(im, 0, 0);

  var P = prep(cv, hgt);
  var x0 = 48, y0 = 10, w = P.w - 62, h = P.h - 40;
  P.g.imageSmoothingEnabled = false;
  P.g.drawImage(off, x0, y0, w, h);
  P.g.imageSmoothingEnabled = true;
  P.g.strokeStyle = '#26405a'; P.g.strokeRect(x0 + 0.5, y0 + 0.5, w, h);
  // axes
  P.g.fillStyle = '#8ba6bd'; P.g.font = '11px monospace';
  P.g.textAlign = 'right'; P.g.textBaseline = 'middle';
  if (axis === 'lat') {
    for (var lat = 20; lat <= 30; lat += 2) {
      var y = y0 + (M.grid.lat1 - lat) / (NLAT * M.grid.d) * h;
      P.g.fillText(lat + '°N', x0 - 6, y);
      P.g.strokeStyle = 'rgba(255,255,255,.15)'; P.g.beginPath(); P.g.moveTo(x0, y); P.g.lineTo(x0 + w, y); P.g.stroke();
    }
  } else {
    for (var lon = 118; lon <= 126; lon += 2) {
      var y2 = y0 + (lon - M.grid.lon0) / (NLON * M.grid.d) * h;
      P.g.fillText(lon + '°E', x0 - 6, y2);
      P.g.strokeStyle = 'rgba(255,255,255,.15)'; P.g.beginPath(); P.g.moveTo(x0, y2); P.g.lineTo(x0 + w, y2); P.g.stroke();
    }
  }
  P.g.textAlign = 'center'; P.g.textBaseline = 'top';
  yearTicks().forEach(function (tk) {
    var xg = x0 + tk.g * w;
    P.g.strokeStyle = 'rgba(255,255,255,.28)'; P.g.beginPath(); P.g.moveTo(xg, y0); P.g.lineTo(xg, y0 + h); P.g.stroke();
    if (tk.l) { P.g.fillStyle = '#8ba6bd'; P.g.fillText(tk.l, x0 + tk.t * w, y0 + h + 6); }
  });
  if (NLEG) {
    var xl = x0 + (tOf(NLEG - 1) + tOf(NLEG)) / 2 * w;
    P.g.save(); P.g.strokeStyle = '#ffb300'; P.g.setLineDash([5, 3]); P.g.lineWidth = 1.4;
    P.g.beginPath(); P.g.moveTo(xl, y0); P.g.lineTo(xl, y0 + h); P.g.stroke(); P.g.restore();
  }
  if (cv.id === 'hovCv') {
    document.getElementById('hovMin').textContent = fmt(lo, 1);
    document.getElementById('hovMax').textContent = fmt(hi, 1);
    drawLegendBar(document.getElementById('hovLg'), function (t) { return ramp(mode === 'sst' ? SST_STOPS : ANO_STOPS, t); });
  }
}

/* =================== MATRIX =================== */
function drawMatrix() {
  var rid = document.getElementById('mxRegion').value, mm = MON[rid];
  var mLab = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'];
  var t = '<thead><tr><th>年 Year</th>' + mLab.map(function (m) { return '<th>' + m + '月</th>'; }).join('') + '<th>年平均</th></tr></thead><tbody>';
  M.years.forEach(function (y, yi) {
    t += '<tr><td class="mono">' + y + '</td>';
    var acc = [], row = mm.v[yi], cnt = mm.n[yi];
    for (var m = 0; m < 12; m++) {
      var v = row[m], n = cnt[m];
      if (v === null || n < 2) { t += '<td style="background:#16232f;color:#3d5468">·</td>'; continue; }
      acc.push(v);
      var thin = n < 10;
      t += '<td class="mono" title="' + n + ' 日樣本' + (thin ? '（樣本少，僅供參考）' : '') + '" style="background:' +
        css(anoColor(v, 2)) + ';color:' + (Math.abs(v) > 1.1 ? '#fff' : '#0d1620') +
        (thin ? ';opacity:.45;font-style:italic' : '') + '">' + sgn(v, 1) + (thin ? '*' : '') + '</td>';
    }
    t += '<td class="mono">' + (acc.length >= 3 ? sgn(meanOf(acc), 2) : '—') + '</td></tr>';
  });
  document.getElementById('mxTable').innerHTML = t + '</tbody>';

  // cross-region last 36 months
  var P = prep(document.getElementById('mxCv'), 330);
  var x0 = 46, y0 = 12, w = P.w - 60, h = P.h - 44;
  var keys = [], labels = [];
  for (var yi2 = 0; yi2 < M.years.length; yi2++) for (var m2 = 0; m2 < 12; m2++) {
    if (MON.all.n[yi2][m2] >= 5) { keys.push([yi2, m2]); labels.push(M.years[yi2] + '-' + String(m2 + 1).padStart(2, '0')); }
  }
  keys = keys.slice(-36); labels = labels.slice(-36);
  var lo = 0, hi = 0;
  M.regions.forEach(function (rg) {
    keys.forEach(function (k) { var v = MON[rg.id].v[k[0]][k[1]]; if (v !== null) { lo = Math.min(lo, v); hi = Math.max(hi, v); } });
  });
  var r = Math.max(Math.abs(lo), Math.abs(hi)) * 1.1 || 1;
  var xt = labels.map(function (l, i) { return { t: keys.length < 2 ? 0 : i / (keys.length - 1), l: i % 6 === 0 ? l.slice(2) : '', grid: i % 6 === 0 }; });
  axes(P.g, x0, y0, w, h, xt, niceTicks(-r, r, 5));
  var leg = '';
  M.regions.forEach(function (rg, ri) {
    var pts = keys.map(function (k, i) {
      var v = MON[rg.id].v[k[0]][k[1]];
      return v === null ? null : [x0 + (keys.length < 2 ? 0 : i / (keys.length - 1)) * w, y0 + h - (v + r) / (2 * r) * h];
    });
    var col = rg.id === 'all' ? '#ffffff' : YCOL[ri % YCOL.length];
    polyline(P.g, pts, col, rg.id === 'all' ? 2.4 : 1.2);
    leg += '<span style="color:' + col + '">■</span> ' + rg.zh + '　';
  });
  document.getElementById('mxLegend').innerHTML = leg;
}

/* =================== MHW =================== */
function drawMHW() {
  var rid = document.getElementById('mhRegion').value, s = SER[rid];
  var ev = s.mhw, tot = 0, mx = null, worst = null;
  ev.forEach(function (e) { tot += e.days; if (!mx || e.days > mx.days) mx = e; if (!worst || e.imax > worst.imax) worst = e; });
  var homoDays = ND - NLEG;
  var homoYrs = Math.max(spanDays(NLEG, ND - 1) / 365.25, 1e-6);
  var k = [
    ['事件總數 Events', ev.length, '平均 ' + (ev.length / homoYrs).toFixed(1) + ' 次/年（基期 ' + homoYrs.toFixed(1) + ' 年）', 'warm'],
    ['累計熱浪日數', tot + ' 日', '占同質期間 ' + (100 * tot / homoDays).toFixed(1) + ' %', 'warm'],
    ['最長事件 Longest', mx ? mx.days + ' 日' : '—', mx ? mx.start + ' ～ ' + mx.end : '', 'warm'],
    ['最強事件 Peak', worst ? '+' + fmt(worst.imax) + ' °C' : '—', worst ? worst.start + '（最高 ' + fmt(worst.tmax) + ' °C）' : '', 'warm']
  ];
  document.getElementById('mhKpi').innerHTML = k.map(function (x) {
    return '<div class="kpi ' + x[3] + '"><div class="lab">' + x[0] + '</div><div class="val">' + x[1] + '</div><div class="sub">' + x[2] + '</div></div>';
  }).join('');

  var P = prep(document.getElementById('mhCv'), 330);
  var x0 = 46, y0 = 12, w = P.w - 60, h = P.h - 42;
  var lo = Infinity, hi = -Infinity;
  for (var i = 0; i < ND; i++) if (s.mean[i] !== null) { lo = Math.min(lo, s.mean[i]); hi = Math.max(hi, s.mean[i]); }
  var pad = (hi - lo) * .05; lo -= pad; hi += pad;
  var xt = yearTicks();
  axes(P.g, x0, y0, w, h, xt, niceTicks(lo, hi, 6));
  shadeLegacy(P.g, x0, y0, w, h);
  function px(i) { return x0 + tOf(i) * w; }
  function py(v) { return y0 + h - (v - lo) / (hi - lo) * h; }
  // fill above threshold
  P.g.save(); P.g.fillStyle = 'rgba(255,112,67,.5)';
  for (var d = NLEG; d < ND; d++) {
    var th = s.p90[DOY[d] - 1], v = s.mean[d];
    if (v === null || th === null || v <= th) continue;
    P.g.fillRect(px(d) - Math.max(w / ND / 2, .5), py(v), Math.max(w / ND, 1.1), py(th) - py(v));
  }
  P.g.restore();
  polyline(P.g, gapBreak(s.mean.map(function (v, i) { return v === null ? null : [px(i), py(v)]; })), 'rgba(139,166,189,.6)', 1);
  polyline(P.g, gapBreak(Array.from({ length: ND }, function (_, i) { var v = s.clim[DOY[i] - 1]; return v === null || v === undefined ? null : [px(i), py(v)]; })), '#66bb6a', 1.4);
  polyline(P.g, gapBreak(Array.from({ length: ND }, function (_, i) { var v = s.p90[DOY[i] - 1]; return v === null || v === undefined ? null : [px(i), py(v)]; })), '#ffb300', 1.2, [4, 3]);

  var t = '<thead><tr><th>#</th><th>起始 Start</th><th>結束 End</th><th>日數 Days</th><th>平均強度 i̅ (°C)</th><th>最大強度 i_max (°C)</th><th>期間最高溫 (°C)</th><th>等級</th></tr></thead><tbody>';
  ev.slice().sort(function (a, b) { return b.days - a.days; }).forEach(function (e, i) {
    var cat = e.imax >= 2 ? '嚴重 Severe' : e.imax >= 1 ? '強 Strong' : '中度 Moderate';
    t += '<tr><td class="mono">' + (i + 1) + '</td><td class="mono">' + e.start + '</td><td class="mono">' + e.end + '</td>' +
      '<td class="mono">' + e.days + '</td><td class="mono">+' + fmt(e.imean) + '</td><td class="mono">+' + fmt(e.imax) + '</td>' +
      '<td class="mono">' + fmt(e.tmax) + '</td><td><span class="pill hot">' + cat + '</span></td></tr>';
  });
  document.getElementById('mhTable').innerHTML = t + '</tbody>';
}

/* =================== QC =================== */
function drawQC() {
  var el = document.getElementById('qcTable');
  if (!QC) { el.innerHTML = '<tbody><tr><td>品管統計檔未提供</td></tr></tbody>'; return; }
  var t = '<thead><tr><th>年 Year</th><th>影像數</th><th>成功反演</th><th>色階 −2~30</th><th>色階 −2~35</th>' +
          '<th>自動判定率</th><th>平均有效格點</th><th>年均海溫 (°C)</th></tr></thead><tbody>';
  QC.years.forEach(function (r) {
    t += '<tr><td class="mono">' + r.y + '</td><td class="mono">' + r.n + '</td><td class="mono">' + r.ok + '</td>' +
      '<td class="mono">' + (r.v30 || '—') + '</td><td class="mono">' + (r.v35 || '—') + '</td>' +
      '<td class="mono">' + (100 * r.auto / Math.max(r.ok, 1)).toFixed(0) + ' %</td>' +
      '<td class="mono">' + r.nvalid.toFixed(0) + '</td><td class="mono">' + fmt(r.tmean) + '</td></tr>';
  });
  var T = QC.total;
  t += '<tr style="font-weight:700;background:#17293a"><td>合計 Total</td><td class="mono">' + T.n + '</td><td class="mono">' + T.ok + '</td>' +
    '<td class="mono">' + T.v30 + '</td><td class="mono">' + T.v35 + '</td><td class="mono">' +
    (100 * T.auto / Math.max(T.ok, 1)).toFixed(1) + ' %</td><td class="mono">' + T.nvalid.toFixed(0) + ' / ' + NC +
    '</td><td class="mono">' + fmt(T.tmean) + '</td></tr>';
  el.innerHTML = t + '</tbody>';

  var d = document.getElementById('qcNotes');
  if (!d) return;
  var rows = [
    ['資料期間 Period', M.dates[0] + ' ～ ' + M.dates[ND - 1] + '（' + ND.toLocaleString() + ' 日）'],
    ['同質基期 Baseline', M.baseline[0] + ' ～ ' + M.baseline[1] + '（' + (ND - NLEG).toLocaleString() + ' 日）'],
    ['舊色階世代 Legacy', M.dates[0] + ' ～ ' + M.legacy_end + '（' + NLEG + ' 日，' +
      '相對基期平均偏差 <b style="color:#ffb300">' + sgn(M.legacy_offset) + ' °C</b>，不納入基期）'],
    ['色階自動判定率', (100 * T.auto / Math.max(T.ok, 1)).toFixed(1) + ' %（其餘以時間鄰近日之判定結果沿用）'],
    ['剔除影像 Rejected', (QC.dropped || 0) + ' 張' + (QC.dropped_list && QC.dropped_list.length ?
      '：' + QC.dropped_list.map(function (x) { return x.date + '（' + x.why + '）'; }).join('、') : '')],
    ['低結構影像 Low-structure', (QC.lowsd || 0) + ' 日（空間標準差 &lt; 0.6 °C）。經目視檢查為原始檔案庫中<b>較粗糙的重繪版本</b>：等溫線為黑色而非紅色、細部渦流結構遺失。此類影像仍保留，但其空間細節不可信，域平均可能偏低。'],
    ['有效格點門檻', '每日 ≥ ' + QC.min_valid + ' 格（實際平均 ' + T.nvalid.toFixed(0) + ' / ' + NC + '）'],
    ['物理值域檢查', QC.trange[0] + ' ~ ' + QC.trange[1] + ' °C，超出者視為疊加物污染剔除']
  ];
  d.innerHTML = '<table>' + rows.map(function (r) {
    return '<tr><td style="width:190px;color:#8ba6bd">' + r[0] + '</td><td>' + r[1] + '</td></tr>';
  }).join('') + '</table>';
}

/* =================== 現場觀測整合 (CTD / SADCP 15 弧分氣候圖集) =================== */
var IS = null;                       // insitu.json
var VT_INFO = {
  mld: ['混合層深度 (m)', '溫度較 10 m 低 0.5 °C 之深度。淺＝表層熱量被侷限在薄層內，升溫快、對表面擾動敏感；深＝熱量向下混合，表水溫變化遲鈍。', 'm', true],
  dT:  ['層化強度 ΔT(10–100 m) (°C)', '10 m 與 100 m 之溫差。值大代表溫躍層強、上下水體交換受阻。', '°C', false],
  swi: ['表層增溫指標 SWI (°C)', '10 m 溫度減去 0–100 m 平均溫度。值大代表熱量集中於表層，衛星看到的高溫「不深」。', '°C', false],
  t10: ['CTD 10 m 溫度 (°C)', '圖集之現場實測 10 dbar 溫度氣候平均。', '°C', false]
};
var GRAD_SEQ = [[0.00,12,24,92],[0.2,20,82,178],[0.42,32,150,200],[0.6,150,214,110],[0.78,238,220,80],[0.9,246,160,48],[1,190,50,35]];

function isGrid() { return { lon: IS.meta.glon, lat: IS.meta.glat }; }

/* 圖集網格之場繪製：格點落在 0.25° 交點，故影像需外推半格以正確套疊 */
function drawAtlas(cv, arr, colf, opt) {
  opt = opt || {};
  var g0 = isGrid(), NX = g0.lon.length, NY = g0.lat.length;
  var dpr = Math.min(window.devicePixelRatio || 1, 2);
  var ASPECT = (NLAT / NLON) / Math.cos(25 * Math.PI / 180);
  var W = cv.parentNode.clientWidth || 640, H = W * ASPECT;
  var maxH = opt.maxH || 620;
  if (H > maxH) { H = maxH; W = H / ASPECT; }
  W = Math.round(W); H = Math.round(H);
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
  cv.style.width = W + 'px'; cv.style.height = H + 'px'; cv.style.margin = '0 auto';
  var g = cv.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.fillStyle = '#0a1219'; g.fillRect(0, 0, W, H);

  if (arr) {
    var im = new ImageData(NX, NY), d = im.data;
    for (var j = 0; j < NY; j++) for (var i = 0; i < NX; i++) {
      var v = arr[NY - 1 - j][i], o = (j * NX + i) * 4;   // 陣列緯度由南而北，影像由北而南
      if (v === null || v === undefined) { d[o + 3] = 0; continue; }
      var c = colf(v);
      d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2]; d[o + 3] = 255;
    }
    var off = document.createElement('canvas'); off.width = NX; off.height = NY;
    off.getContext('2d').putImageData(im, 0, 0);
    // 半格外推：資料涵蓋 115.875–128.125 / 17.875–32.125
    var sx = W / 12, sy = H / 14;
    g.imageSmoothingEnabled = opt.smooth !== false; g.imageSmoothingQuality = 'high';
    g.drawImage(off, -0.125 * sx, -0.125 * sy, 12.25 * sx, 14.25 * sy);
  }
  if (LANDCV) { g.imageSmoothingEnabled = true; g.drawImage(LANDCV, 0, 0, W, H); }

  g.save();
  g.strokeStyle = 'rgba(255,255,255,.14)'; g.lineWidth = 1; g.font = '10px monospace';
  g.shadowColor = 'rgba(0,0,0,.85)'; g.shadowBlur = 3; g.fillStyle = 'rgba(255,255,255,.9)';
  for (var lon = 118; lon <= 126; lon += 2) {
    var x = (lon - 116) / 12 * W;
    g.beginPath(); g.moveTo(x, 0); g.lineTo(x, H); g.stroke();
    g.textAlign = 'left'; g.textBaseline = 'bottom'; g.fillText(lon + '°E', x + 3, H - 3);
  }
  for (var lat = 20; lat <= 30; lat += 2) {
    var y = (32 - lat) / 14 * H;
    g.beginPath(); g.moveTo(0, y); g.lineTo(W, y); g.stroke();
    g.textAlign = 'left'; g.textBaseline = 'bottom'; g.fillText(lat + '°N', 3, y - 2);
  }
  g.restore();
  g.strokeStyle = 'rgba(139,166,189,.45)'; g.lineWidth = 1; g.strokeRect(0.5, 0.5, W - 1, H - 1);
  return { g: g, W: W, H: H };
}

function arrRange(arr, plo, phi) {
  var v = [];
  for (var j = 0; j < arr.length; j++) for (var i = 0; i < arr[j].length; i++)
    if (arr[j][i] !== null && arr[j][i] !== undefined) v.push(arr[j][i]);
  if (!v.length) return [0, 1];
  v.sort(function (a, b) { return a - b; });
  var q = function (p) { return v[Math.min(v.length - 1, Math.max(0, Math.round(p * (v.length - 1))))]; };
  return [q(plo === undefined ? 0.02 : plo), q(phi === undefined ? 0.98 : phi)];
}

/* ---------- 分頁：現場驗證 ---------- */
function drawValid() {
  if (!IS) return;
  var S = IS.meta.seasons, byTp = {};
  S.forEach(function (r) { byTp[r.tp] = r; });
  document.getElementById('vaBase').textContent = IS.meta.baseline[0] + ' ～ ' + IS.meta.baseline[1];

  var w = byTp[13], su = byTp[15], an = byTp[0];
  var k = [
    ['冬季偏差 DJF bias', sgn(w.bias) + ' °C', 'RMSE ' + fmt(w.rmse) + '　R = ' + w.r.toFixed(3), Math.abs(w.bias) < 0.2 ? 'acc' : 'warm'],
    ['夏季偏差 JJA bias', sgn(su.bias) + ' °C', 'RMSE ' + fmt(su.rmse) + '　R = ' + su.r.toFixed(3), 'warm'],
    ['冬季混合層 DJF MLD', fmt(w.mld, 0) + ' m', '層化 ΔT ' + fmt(w.dT) + ' °C（弱）', 'cool'],
    ['夏季混合層 JJA MLD', fmt(su.mld, 0) + ' m', '層化 ΔT ' + fmt(su.dT) + ' °C（強）', 'warm'],
    ['比對格點數', an.n.toLocaleString(), '0.25° 圖集格點（年平均）', 'acc'],
    ['CTD 資料列', IS.meta.nctd.toLocaleString(), '15 弧分網格氣候圖集', 'acc'],
    ['SADCP 資料列', IS.meta.nadcp.toLocaleString(), '0–500 m 流速剖面', 'acc']
  ];
  document.getElementById('vaKpi').innerHTML = k.map(function (x) {
    return '<div class="kpi ' + x[3] + '"><div class="lab">' + x[0] + '</div><div class="val">' + x[1] + '</div><div class="sub">' + x[2] + '</div></div>';
  }).join('');

  var t = '<thead><tr><th>期間 Period</th><th>格點數 n</th><th>偏差 bias (°C)</th><th>RMSE (°C)</th>'
        + '<th>中位偏差</th><th>標準差 SD</th><th>相關 R</th><th>MLD 中位 (m)</th><th>ΔT 中位 (°C)</th><th>衛星日數</th></tr></thead><tbody>';
  S.forEach(function (r) {
    var hot = Math.abs(r.bias) >= 0.4;
    t += '<tr><td>' + r.name + '</td><td class="mono">' + r.n + '</td>'
      + '<td class="mono" style="color:' + (hot ? '#ff7043' : '#66bb6a') + '">' + sgn(r.bias) + '</td>'
      + '<td class="mono">' + fmt(r.rmse) + '</td><td class="mono">' + sgn(r.median) + '</td>'
      + '<td class="mono">' + fmt(r.sd) + '</td><td class="mono">' + r.r.toFixed(3) + '</td>'
      + '<td class="mono">' + fmt(r.mld, 0) + '</td><td class="mono">' + fmt(r.dT) + '</td>'
      + '<td class="mono">' + r.satdays + '</td></tr>';
  });
  document.getElementById('vaTable').innerHTML = t + '</tbody>';

  drawScatter(); drawBinned(); drawBiasMap();
}

function drawScatter() {
  var tp = segVal('vaScMode'), pts = IS.scatter[tp] || [];
  var P = prep(document.getElementById('vaScCv'), 360);
  var x0 = 46, y0 = 12, w = P.w - 60, h = P.h - 44;
  var all = [];
  pts.forEach(function (p) { all.push(p[0]); all.push(p[1]); });
  var lo = Math.floor(Math.min.apply(null, all) - 0.5), hi = Math.ceil(Math.max.apply(null, all) + 0.5);
  var tk = niceTicks(lo, hi, 6);
  axes(P.g, x0, y0, w, h, tk, tk);
  var X = function (v) { return x0 + (v - lo) / (hi - lo) * w; };
  var Y = function (v) { return y0 + h - (v - lo) / (hi - lo) * h; };
  polyline(P.g, [[X(lo), Y(lo)], [X(hi), Y(hi)]], 'rgba(139,166,189,.6)', 1.2, [5, 4]);
  P.g.fillStyle = 'rgba(79,195,247,.55)';
  pts.forEach(function (p) {
    P.g.beginPath(); P.g.arc(X(p[0]), Y(p[1]), 2.2, 0, 7); P.g.fill();
  });
  P.g.fillStyle = '#8ba6bd'; P.g.font = '11px sans-serif';
  P.g.textAlign = 'center'; P.g.fillText('CTD 10 m (°C)', x0 + w / 2, y0 + h + 24);
  P.g.save(); P.g.translate(12, y0 + h / 2); P.g.rotate(-Math.PI / 2);
  P.g.fillText('衛星反演 SST (°C)', 0, 0); P.g.restore();
  var r = IS.meta.seasons.filter(function (s) { return String(s.tp) === String(tp); })[0];
  P.g.textAlign = 'left'; P.g.fillStyle = '#e8f1f8';
  P.g.fillText('n=' + r.n + '　bias ' + sgn(r.bias) + '　RMSE ' + fmt(r.rmse) + '　R ' + r.r.toFixed(3), x0 + 8, y0 + 14);
}

function drawBinned() {
  var P = prep(document.getElementById('vaBinCv'), 300);
  var x0 = 52, y0 = 14, w = P.w - 70, h = P.h - 52;
  var byTp = {}; IS.meta.seasons.forEach(function (r) { byTp[r.tp] = r; });
  var W13 = byTp[13].binned, W15 = byTp[15].binned;
  var vals = W13.concat(W15).map(function (b) { return b.bias; }).filter(function (v) { return v !== null; });
  var hi = Math.max(1.2, Math.ceil(Math.max.apply(null, vals) * 10) / 10 + 0.1);
  var lo = Math.min(-0.4, Math.floor(Math.min.apply(null, vals) * 10) / 10 - 0.1);
  var n = W15.length;
  var xt = W15.map(function (b, i) { return { t: (i + 0.5) / n, l: b.lab + ' m', grid: false }; });
  axes(P.g, x0, y0, w, h, xt, niceTicks(lo, hi, 6));
  var zy = y0 + h - (0 - lo) / (hi - lo) * h;
  P.g.strokeStyle = '#40607e'; P.g.beginPath(); P.g.moveTo(x0, zy); P.g.lineTo(x0 + w, zy); P.g.stroke();
  var bw = w / n * 0.34;
  [[W13, '#4fc3f7', -1], [W15, '#ff7043', 1]].forEach(function (S) {
    S[0].forEach(function (b, i) {
      if (b.bias === null) return;
      var cx = x0 + (i + 0.5) / n * w + S[2] * bw * 0.55;
      var y = y0 + h - (b.bias - lo) / (hi - lo) * h;
      P.g.fillStyle = S[1];
      P.g.fillRect(cx - bw / 2, Math.min(y, zy), bw, Math.abs(zy - y));
      P.g.fillStyle = '#8ba6bd'; P.g.font = '10px monospace'; P.g.textAlign = 'center';
      P.g.fillText(b.n, cx, (b.bias >= 0 ? y - 4 : y + 12));
    });
  });
  P.g.fillStyle = '#8ba6bd'; P.g.font = '11px sans-serif'; P.g.textAlign = 'left';
  P.g.fillText('偏差 (°C)', 6, y0 + 4);
}

function drawBiasMap() {
  var sel = document.getElementById('vaMapSea');
  var tp = sel.value || '15';
  var arr = IS.fields[tp].bias;
  drawAtlas(document.getElementById('vaMapCv'), arr,
    function (v) { return anoColor(v, 2); }, { smooth: false, maxH: 420 });
  drawLegendBar(document.getElementById('vaMapLg'), function (t) { return ramp(ANO_STOPS, t); });
}

/* ---------- 分頁：垂直結構與海流 ---------- */
function drawVert() {
  if (!IS) return;
  var tp = document.getElementById('vtSea').value || '15';
  var f = segVal('vtField');
  var arr = IS.fields[tp][f], info = VT_INFO[f];
  var rg = arrRange(arr, 0.03, 0.97);
  var inv = info[3];                     // MLD：淺為暖色（風險高）
  var colf = function (v) {
    var t = (v - rg[0]) / (rg[1] - rg[0]);
    return ramp(GRAD_SEQ, inv ? 1 - t : t);
  };
  drawAtlas(document.getElementById('vtMapCv'), arr, colf, { smooth: false, maxH: 560 });
  drawLegendBar(document.getElementById('vtLg'), function (t) { return ramp(GRAD_SEQ, inv ? 1 - t : t); });
  document.getElementById('vtLo').textContent = fmt(rg[0], f === 'mld' ? 0 : 1);
  document.getElementById('vtHi').textContent = fmt(rg[1], f === 'mld' ? 0 : 1);
  document.getElementById('vtUnit').textContent = info[2];
  document.getElementById('vtNote').innerHTML = '<b>' + info[0] + '</b><br>' + info[1];

  var v = [];
  for (var j = 0; j < arr.length; j++) for (var i = 0; i < arr[j].length; i++)
    if (arr[j][i] !== null) v.push(arr[j][i]);
  v.sort(function (a, b) { return a - b; });
  var q = function (p) { return v[Math.round(p * (v.length - 1))]; };
  var dec = f === 'mld' ? 0 : 2;
  var rows = [['有效格點', v.length], ['最小值', fmt(q(0), dec)], ['10 百分位', fmt(q(0.1), dec)],
              ['中位數', fmt(q(0.5), dec)], ['90 百分位', fmt(q(0.9), dec)], ['最大值', fmt(q(1), dec)]];
  document.getElementById('vtStats').innerHTML = rows.map(function (r) {
    return '<tr><td>' + r[0] + '</td><td class="mono">' + r[1] + '</td></tr>';
  }).join('');
}

function drawCurrents() {
  var tp = document.getElementById('cuSea').value || '0';
  var C0 = IS.cur[tp], ax = IS.axis[tp] || [];
  var showAxis = document.getElementById('cuAxis').classList.contains('on');
  var showSST = document.getElementById('cuSST').classList.contains('on');
  var base = null, rg = null;
  if (showSST) {
    base = IS.fields[tp].sat || IS.fields[tp].t10; rg = arrRange(base, 0.02, 0.98);
  }
  var R = drawAtlas(document.getElementById('cuMapCv'), base,
    base ? function (v) { return ramp(SST_STOPS, (v - rg[0]) / (rg[1] - rg[0])); } : null,
    { smooth: true, maxH: 560 });
  var g = R.g, W = R.W, H = R.H, G = isGrid();
  var X = function (lon) { return (lon - 116) / 12 * W; };
  var Y = function (lat) { return (32 - lat) / 14 * H; };

  var step = W < 420 ? 3 : 2, sc = Math.min(W / 12, 60) * 0.55;
  g.save(); g.lineWidth = 1.1; g.strokeStyle = 'rgba(232,241,248,.88)'; g.fillStyle = 'rgba(232,241,248,.88)';
  for (var j = 0; j < G.lat.length; j += step) for (var i = 0; i < G.lon.length; i += step) {
    var u = C0.u[j][i], v = C0.v[j][i];
    if (u === null || v === null) continue;
    var sp = Math.hypot(u, v); if (sp < 0.02) continue;
    var x = X(G.lon[i]), y = Y(G.lat[j]);
    var dx = u * sc, dy = -v * sc;
    g.beginPath(); g.moveTo(x, y); g.lineTo(x + dx, y + dy); g.stroke();
    var a = Math.atan2(dy, dx), hl = Math.min(4.5, 2 + sp * 4);
    g.beginPath(); g.moveTo(x + dx, y + dy);
    g.lineTo(x + dx - hl * Math.cos(a - 0.4), y + dy - hl * Math.sin(a - 0.4));
    g.lineTo(x + dx - hl * Math.cos(a + 0.4), y + dy - hl * Math.sin(a + 0.4));
    g.closePath(); g.fill();
  }
  g.restore();

  if (showAxis && ax.length) {
    polyline(g, ax.map(function (a) { return [X(a.lon), Y(a.lat)]; }), '#ffb300', 2.4);
    g.fillStyle = '#ffb300';
    ax.forEach(function (a) { g.beginPath(); g.arc(X(a.lon), Y(a.lat), 2.6, 0, 7); g.fill(); });
  }

  var t = '<thead><tr><th>緯度</th><th>軸心經度</th><th>北向流 v (m/s)</th><th>流速 |V|</th><th>軸心 SST</th><th>西側 1° SST</th><th>ΔSST</th></tr></thead><tbody>';
  ax.forEach(function (a) {
    var d = (a.sst !== null && a.sstW !== null) ? a.sst - a.sstW : null;
    t += '<tr><td class="mono">' + a.lat.toFixed(2) + '°N</td><td class="mono">' + a.lon.toFixed(2) + '°E</td>'
      + '<td class="mono">' + fmt(a.v, 3) + '</td><td class="mono">' + fmt(a.spd, 3) + '</td>'
      + '<td class="mono">' + fmt(a.sst) + '</td><td class="mono">' + fmt(a.sstW) + '</td>'
      + '<td class="mono" style="color:' + (d > 0 ? '#ff7043' : '#42a5f5') + '">' + (d === null ? '—' : sgn(d)) + '</td></tr>';
  });
  document.getElementById('cuTable').innerHTML = t + '</tbody>';
}
