#!/usr/bin/env python3
"""
sgcal - ColorChecker SG based correction for FADGI / Open DICE workflows.

  fit    : read a photo of the ColorChecker SG, find the chart automatically,
           fit  (1) per-channel tone / gray-balance curves  and
                (2) a neutral-preserving root-polynomial colour transform,
           write a corrected copy of the SG shot, a reusable correction.json,
           a pass/fail report and an overlay image showing where patches were read.
  apply  : apply a saved correction.json to any number of images that were shot
           in the same session (same light, lens, exposure, same colour space).

The measurements mimic Open DICE (D50 Lab from the embedded RGB profile). They are
very close to Open DICE's numbers but not identical, so the report warns when a
result is close to a limit.
"""
import argparse
import csv
import datetime
import json
import os
import struct
import sys

import numpy as np
import tifffile
from scipy import ndimage
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares

VERSION = "1.0"
COLS = "ABCDEFGHIJKLMN"
NCOLS, NROWS = 14, 10
PATCH_NAMES = [f"{c}{r}" for r in range(1, NROWS + 1) for c in COLS]

# FADGI defaults (General materials) per star level: tone band, white balance, mean dE, 90th dE
DEFAULT_LIMITS = {
    4: dict(tone=1.5, wb=2.0, mean=2.0, p90=4.0),
    3: dict(tone=3.0, wb=4.0, mean=3.5, p90=7.0),
    2: dict(tone=4.5, wb=6.0, mean=5.0, p90=10.0),
    1: dict(tone=6.0, wb=8.0, mean=6.5, p90=13.0),
}


class SgcalError(Exception):
    pass


# --------------------------------------------------------------------------- colour spaces
_BUILTIN = {
    # RGB -> XYZ (D50 adapted), gamma / curve
    "adobergb": dict(
        name="Adobe RGB (1998)",
        M=[[0.5767309, 0.1855540, 0.1881852], [0.2973769, 0.6273491, 0.0752741], [0.0270343, 0.0706872, 0.9911085]],
        curve=("gamma", 2.19921875)),
    "srgb": dict(
        name="sRGB",
        M=[[0.4360747, 0.3850649, 0.1430804], [0.2225045, 0.7168786, 0.0606169], [0.0139322, 0.0971045, 0.7141733]],
        curve=("srgb", 0)),
    "prophoto": dict(
        name="ProPhoto RGB",
        M=[[0.7976749, 0.1351917, 0.0313534], [0.2880402, 0.7118741, 0.0000857], [0.0, 0.0, 0.8252100]],
        curve=("gamma", 1.8)),
}


def _srgb_decode(x):
    x = np.asarray(x, dtype=np.float64)
    return np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)


def _parse_icc(data):
    n = struct.unpack(">I", data[128:132])[0]
    tags = {}
    for i in range(n):
        sig, off, size = struct.unpack(">4sII", data[132 + 12 * i:144 + 12 * i])
        tags[sig.decode("latin1")] = data[off:off + size]
    return tags


def _icc_xyz(b):
    return np.array(struct.unpack(">3i", b[8:20])) / 65536.0


def _icc_curve(b):
    t = b[:4]
    if t == b"curv":
        n = struct.unpack(">I", b[8:12])[0]
        if n == 0:
            return lambda x: np.asarray(x, dtype=np.float64)
        if n == 1:
            g = struct.unpack(">H", b[12:14])[0] / 256.0
            return lambda x: np.power(np.clip(x, 0, 1), g)
        tab = np.array(struct.unpack(">%dH" % n, b[12:12 + 2 * n])) / 65535.0
        xs = np.linspace(0, 1, n)
        return lambda x: np.interp(x, xs, tab)
    if t == b"para":
        ft = struct.unpack(">H", b[8:10])[0]
        npar = {0: 1, 1: 3, 2: 4, 3: 5, 4: 7}[ft]
        p = np.array(struct.unpack(">%di" % npar, b[12:12 + 4 * npar])) / 65536.0
        g = p[0]
        if ft == 0:
            return lambda x: np.power(np.clip(x, 0, 1), g)
        a, bb = p[1], p[2]
        if ft == 1:
            return lambda x: np.where(x >= -bb / a, np.power(np.maximum(a * x + bb, 0), g), 0.0)
        c = p[3]
        if ft == 2:
            return lambda x: np.where(x >= -bb / a, np.power(np.maximum(a * x + bb, 0), g) + c, c)
        d = p[4]
        if ft == 3:
            return lambda x: np.where(x >= d, np.power(np.maximum(a * x + bb, 0), g), c * x)
        e, f = p[5], p[6]
        return lambda x: np.where(x >= d, np.power(np.maximum(a * x + bb, 0), g) + e, c * x + f)
    raise SgcalError("Unsupported ICC curve type in the image profile.")


class ColorSpace:
    """Matrix/TRC RGB space. Linear RGB <-> XYZ(D50) <-> Lab(D50)."""

    def __init__(self, name, M, curves):
        self.name = name
        self.M = np.array(M, dtype=np.float64)
        self.Minv = np.linalg.inv(self.M)
        self.curves = curves  # 3 callables: encoded[0..1] -> linear
        self.Wn = self.M @ np.ones(3)  # neutral axis == white point of the space
        self._dec = {}

    @classmethod
    def from_builtin(cls, key):
        d = _BUILTIN[key]
        kind, g = d["curve"]
        fn = _srgb_decode if kind == "srgb" else (lambda x, g=g: np.power(np.clip(x, 0, 1), g))
        return cls(d["name"], d["M"], [fn] * 3)

    @classmethod
    def from_icc(cls, data):
        tags = _parse_icc(data)
        if not all(k in tags for k in ("rXYZ", "gXYZ", "bXYZ", "rTRC", "gTRC", "bTRC")):
            raise SgcalError("The embedded ICC profile is not a matrix/TRC RGB profile. "
                             "Re-export the TIFF as Adobe RGB, sRGB or ProPhoto, or pass --colorspace.")
        M = np.stack([_icc_xyz(tags["rXYZ"]), _icc_xyz(tags["gXYZ"]), _icc_xyz(tags["bXYZ"])], axis=1)
        curves = [_icc_curve(tags["rTRC"]), _icc_curve(tags["gTRC"]), _icc_curve(tags["bTRC"])]
        name = "ICC profile"
        if "desc" in tags:
            b = tags["desc"]
            try:
                if b[:4] == b"desc":
                    n = struct.unpack(">I", b[8:12])[0]
                    name = b[12:12 + n].split(b"\0")[0].decode("latin1")
                elif b[:4] == b"mluc":
                    cnt = struct.unpack(">I", b[8:12])[0]
                    ln, off = struct.unpack(">II", b[20:28])
                    name = b[off:off + ln].decode("utf-16-be")
            except Exception:
                pass
        return cls(name, M, curves)

    # code (0..N-1) -> linear, per channel, as float32 tables
    def decode_lut(self, bits):
        if bits not in self._dec:
            n = 2 ** bits
            x = np.arange(n, dtype=np.float64) / (n - 1)
            self._dec[bits] = [np.asarray(c(x), dtype=np.float64) for c in self.curves]
        return self._dec[bits]

    def encode(self, lin, ch, bits):
        """linear[0..1] -> integer code, for channel ch, via an inverse of the decode table."""
        dec16 = self._dec16(ch)
        n = 2 ** bits - 1
        code = np.interp(lin, dec16, np.linspace(0.0, 1.0, dec16.size))
        return np.rint(code * n)

    def _dec16(self, ch):
        key = ("d16", ch)
        if key not in self._dec:
            x = np.linspace(0, 1, 65536)
            d = np.asarray(self.curves[ch](x), dtype=np.float64)
            self._dec[key] = np.maximum.accumulate(d)  # guarantee monotone for np.interp
        return self._dec[key]

    def lab(self, lin):
        """lin: (N,3) linear RGB -> Lab (N,3) D50 relative to this space's white."""
        xyz = np.clip(lin, 0, 1) @ self.M.T / self.Wn
        f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16 / 116)
        return np.stack([116 * f[:, 1] - 16, 500 * (f[:, 0] - f[:, 1]), 200 * (f[:, 1] - f[:, 2])], 1)

    def lab2lin(self, L, a, b):
        fy = (L + 16) / 116
        f = np.array([fy + a / 500, fy, fy - b / 200])
        xyz = np.where(f ** 3 > 0.008856, f ** 3, (f - 16 / 116) / 7.787) * self.Wn
        return self.Minv @ xyz

    def same_as(self, other):
        return np.allclose(self.M, other.M, atol=2e-3)


def Lstar_from_density(d):
    y = 10.0 ** (-np.asarray(d, dtype=np.float64))
    return np.where(y > 0.008856, 116 * np.cbrt(y) - 16, 903.3 * y)


# --------------------------------------------------------------------------- CIEDE2000
def deltaE2000(lab1, lab2):
    lab1 = np.atleast_2d(np.asarray(lab1, float))
    lab2 = np.atleast_2d(np.asarray(lab2, float))
    L1, a1, b1 = lab1[:, 0], lab1[:, 1], lab1[:, 2]
    L2, a2, b2 = lab2[:, 0], lab2[:, 1], lab2[:, 2]
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2)
    Cm = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cm ** 7 / (Cm ** 7 + 25.0 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360
    dLp = L2 - L1
    dCp = C2p - C1p
    dh = h2p - h1p
    dh = np.where(dh > 180, dh - 360, np.where(dh < -180, dh + 360, dh))
    dh = np.where(C1p * C2p == 0, 0, dh)
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dh / 2))
    Lpm = (L1 + L2) / 2
    Cpm = (C1p + C2p) / 2
    hs = h1p + h2p
    hpm = np.where(np.abs(h1p - h2p) <= 180, hs / 2, np.where(hs < 360, (hs + 360) / 2, (hs - 360) / 2))
    hpm = np.where(C1p * C2p == 0, hs, hpm)
    T = (1 - 0.17 * np.cos(np.radians(hpm - 30)) + 0.24 * np.cos(np.radians(2 * hpm))
         + 0.32 * np.cos(np.radians(3 * hpm + 6)) - 0.20 * np.cos(np.radians(4 * hpm - 63)))
    dth = 30 * np.exp(-(((hpm - 275) / 25) ** 2))
    Rc = 2 * np.sqrt(Cpm ** 7 / (Cpm ** 7 + 25.0 ** 7))
    Sl = 1 + 0.015 * (Lpm - 50) ** 2 / np.sqrt(20 + (Lpm - 50) ** 2)
    Sc = 1 + 0.045 * Cpm
    Sh = 1 + 0.015 * Cpm * T
    Rt = -np.sin(np.radians(2 * dth)) * Rc
    return np.sqrt((dLp / Sl) ** 2 + (dCp / Sc) ** 2 + (dHp / Sh) ** 2 + Rt * (dCp / Sc) * (dHp / Sh))


# --------------------------------------------------------------------------- reference data
def read_reference(xlsx_path):
    """Profile_ColorCheckerSG xlsx: col A name, B L*, C a*, D b*, E density (tone patches only)."""
    import openpyxl
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb.worksheets[0]
    ref, dens = {}, {}
    for row in ws.iter_rows(values_only=True):
        if not row or not isinstance(row[0], str) or row[0] not in PATCH_NAMES:
            continue
        if None in row[1:4]:
            continue
        ref[row[0]] = np.array(row[1:4], dtype=float)
        if len(row) > 4 and isinstance(row[4], (int, float)):
            dens[row[0]] = float(row[4])
    missing = [n for n in PATCH_NAMES if n not in ref]
    if missing:
        raise SgcalError("Profile sheet is missing patches: " + ", ".join(missing[:10]))
    if len(dens) < 6:
        raise SgcalError("Profile sheet has no density column (column E) for the tone patches.")
    return ref, dens


def read_limits(config_xlsx, material, level):
    lim = dict(DEFAULT_LIMITS[level])
    src = "built-in FADGI defaults"
    if config_xlsx:
        import openpyxl
        wb = openpyxl.load_workbook(config_xlsx, data_only=True)
        if material not in wb.sheetnames:
            raise SgcalError("Material '%s' not in config file. Available: %s" % (material, ", ".join(wb.sheetnames)))
        ws = wb[material]
        rows = [r for r in ws.iter_rows(values_only=True) if r and r[0]]
        header = next(r for r in rows if str(r[0]).strip() == "FADGI Criteria")
        want = "*" * level
        heads = [str(c).strip() if c is not None else "" for c in header]
        if want not in heads:
            raise SgcalError("Material '%s' has no %d-star column." % (material, level))
        col = heads.index(want)
        key = {"OECF band": "tone", "White Balance": "wb", "Mean DeltaE": "mean", "90th DeltaE": "p90"}
        for r in rows:
            k = key.get(str(r[0]).strip())
            if k and r[col] is not None:
                lim[k] = float(r[col])
        src = "%s / %s / %d-star" % (os.path.basename(config_xlsx), material, level)
    lim["source"] = src
    return lim


def list_materials(config_xlsx):
    import openpyxl
    return openpyxl.load_workbook(config_xlsx, data_only=True).sheetnames


# --------------------------------------------------------------------------- image IO
class Img:
    def __init__(self, path, colorspace_override=None, lazy=False):
        with tifffile.TiffFile(path) as t:
            pg = t.pages[0]
            try:
                self.a = pg.asarray(out="memmap") if lazy else pg.asarray()
            except Exception:
                try:
                    self.a = pg.asarray()
                except ValueError as e:
                    if "imagecodecs" in str(e):
                        raise SgcalError("This TIFF is compressed (%s). Install the decoder once with:  py -3 -m pip install imagecodecs"
                                         % getattr(pg.compression, "name", pg.compression))
                    raise
            tg = pg.tags
            self.icc = tg["InterColorProfile"].value if "InterColorProfile" in tg else None
            xr = tg["XResolution"].value if "XResolution" in tg else None
            yr = tg["YResolution"].value if "YResolution" in tg else None
            self.res = (xr[0] / xr[1] if xr else 300.0, yr[0] / yr[1] if yr else 300.0)
            ru = tg["ResolutionUnit"].value if "ResolutionUnit" in tg else 2
            self.resunit = {1: "none", 2: "inch", 3: "centimeter"}.get(int(ru), "inch")
            self.software = str(tg["Software"].value) if "Software" in tg else ""
        if self.a.ndim != 3 or self.a.shape[2] < 3:
            raise SgcalError("Need an RGB TIFF (got shape %s)." % (self.a.shape,))
        self.a = self.a[..., :3]
        if self.a.dtype == np.uint8:
            self.bits = 8
        elif self.a.dtype == np.uint16:
            self.bits = 16
        else:
            raise SgcalError("Only 8-bit and 16-bit integer TIFFs are supported (got %s)." % self.a.dtype)
        if colorspace_override:
            self.cs = ColorSpace.from_builtin(colorspace_override)
        elif self.icc:
            self.cs = ColorSpace.from_icc(bytes(self.icc))
        else:
            self.cs = ColorSpace.from_builtin("srgb")
            print("  note: no embedded ICC profile, assuming sRGB (use --colorspace to change).")


def write_tiff(path, arr, src, compress=False):
    extra = []
    if src.icc:
        b = bytes(src.icc)
        extra.append((34675, "B", len(b), b, True))
    tifffile.imwrite(path, arr, photometric="rgb", resolution=src.res, resolutionunit=src.resunit,
                     extratags=extra, compression=("zlib" if compress else None), predictor=bool(compress),
                     bigtiff=arr.nbytes > 3.5e9,
                     software="sgcal %s (corrected)" % VERSION)


def check_not_corrected(img, path, force=False):
    if img.software.startswith("sgcal") and not force:
        raise SgcalError("%s was already produced by sgcal. Applying the correction twice would over-correct it. "
                         "Use the original (uncorrected) file, or pass --force." % os.path.basename(path))


# --------------------------------------------------------------------------- chart detection
def _gray_small(img):
    a = img.a
    H, W = a.shape[:2]
    f = max(1, int(max(H, W) / 1600))
    h, w = H // f * f, W // f * f
    out = np.empty((h // f, w // f), np.float32)
    step = f * 64
    for y in range(0, h, step):
        y1 = min(y + step, h)
        blk = a[y:y1, :w].astype(np.float32)
        out[y // f:y1 // f] = blk.reshape((y1 - y) // f, f, w // f, f, 3).mean(axis=(1, 3, 4))
    return out / float(2 ** img.bits - 1), f


def _homography(src, dst):
    A, b = [], []
    for (x, y), (u, v) in zip(src, dst):
        A.append([x, y, 1, 0, 0, 0, -u * x, -u * y]); b.append(u)
        A.append([0, 0, 0, x, y, 1, -v * x, -v * y]); b.append(v)
    h = np.linalg.solve(np.array(A, float), np.array(b, float))
    return np.append(h, 1).reshape(3, 3)


def _apply_h(Hm, pts):
    p = np.c_[pts, np.ones(len(pts))] @ Hm.T
    return p[:, :2] / p[:, 2:3]


def find_corner_patches(img):
    """Return 4 centres (full-res xy) of the white corner patches A1,N1,A10,N10 as TL,TR,BL,BR."""
    g, f = _gray_small(img)
    hi = np.percentile(g, 99.5)
    thr = 0.6 * hi
    lab, n = ndimage.label(g > thr)
    if n == 0:
        raise SgcalError("Could not find any bright patches. Is the chart in the frame and exposed?")
    objs = ndimage.find_objects(lab)
    areas = ndimage.sum(np.ones_like(g), lab, index=np.arange(1, n + 1))
    cand = []
    for i, sl in enumerate(objs):
        bh, bw = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        if bw < 5 or bh < 5:
            continue
        fill = areas[i] / (bh * bw)
        asp = bw / bh
        if fill > 0.75 and 0.75 < asp < 1.33:
            cand.append((i + 1, areas[i]))
    if len(cand) < 8:
        raise SgcalError("Found too few square patches (%d) - chart not detected." % len(cand))
    med = np.median([a for _, a in cand])
    cand = [(i, a) for i, a in cand if 0.5 * med <= a <= 1.8 * med]
    if len(cand) < 6:
        raise SgcalError("Chart patches not consistent in size - chart not detected.")
    cen = np.array(ndimage.center_of_mass(np.ones_like(g), lab, [i for i, _ in cand]))[:, ::-1]  # x,y
    s, d = cen.sum(1), cen[:, 0] - cen[:, 1]
    tl, br = cen[np.argmin(s)], cen[np.argmax(s)]
    tr, bl = cen[np.argmax(d)], cen[np.argmin(d)]
    pts = np.array([tl, tr, bl, br]) * f + f / 2.0
    return pts, float(np.sqrt(med)) * f


def grid_centres(corners):
    """corners: TL,TR,BL,BR centres -> (140,2) patch centres in name order, plus pitch."""
    hm = _homography(np.array([[0, 0], [13, 0], [0, 9], [13, 9]], float), corners)
    gx, gy = np.meshgrid(np.arange(NCOLS), np.arange(NROWS))
    pts = _apply_h(hm, np.c_[gx.ravel(), gy.ravel()].astype(float))
    pitch = np.mean([np.linalg.norm(corners[1] - corners[0]) / 13, np.linalg.norm(corners[3] - corners[2]) / 13,
                     np.linalg.norm(corners[2] - corners[0]) / 9, np.linalg.norm(corners[3] - corners[1]) / 9])
    return pts, float(pitch)


def read_patches(a, centres, pitch, frac=0.15):
    """Median code per channel inside a square of half-width frac*pitch. Returns (N,3) codes (float)."""
    h = max(2, int(round(frac * pitch)))
    H, W = a.shape[:2]
    out = np.zeros((len(centres), 3))
    for i, (x, y) in enumerate(centres):
        x, y = int(round(x)), int(round(y))
        x0, x1, y0, y1 = max(0, x - h), min(W, x + h + 1), max(0, y - h), min(H, y + h + 1)
        if x1 - x0 < 3 or y1 - y0 < 3:
            raise SgcalError("Patch %s falls outside the image." % PATCH_NAMES[i])
        out[i] = np.median(a[y0:y1, x0:x1].reshape(-1, 3), axis=0)
    return out


def detect_and_read(img, ref, corners=None, log=print):
    """Find the chart (or use given corners), choose the orientation that best matches the reference."""
    if corners is None:
        corners, _ = find_corner_patches(img)
    else:
        corners = np.array(corners, float).reshape(4, 2)
    refL = np.array([ref[n][0] for n in PATCH_NAMES])
    best = None
    for rot in (0, 180):
        c = corners if rot == 0 else corners[::-1]  # 180 deg: TL<->BR, TR<->BL
        centres, pitch = grid_centres(c)
        codes = read_patches(img.a, centres, pitch)
        lin = np.stack([img.cs.decode_lut(img.bits)[k][codes[:, k].astype(int)] for k in range(3)], 1)
        L = img.cs.lab(lin)[:, 0]
        r = np.corrcoef(L, refL)[0, 1]
        if best is None or r > best[0]:
            best = (r, rot, centres, pitch, codes, lin)
    r, rot, centres, pitch, codes, lin = best
    log("  chart found: pitch %.1f px, orientation %d deg, match to reference r=%.3f" % (pitch, rot, r))
    if r < 0.97:
        raise SgcalError("Chart detected but patch values do not match the reference (r=%.2f). "
                         "Check exposure/glare, or give the four white corner patch centres with --corners "
                         "x1,y1,x2,y2,x3,y3,x4,y4 (top-left, top-right, bottom-left, bottom-right)." % r)
    return dict(centres=centres, pitch=pitch, codes=codes, lin=lin, corr=r, rot=rot)


def write_overlay(img, det, path):
    from PIL import Image, ImageDraw
    a = img.a
    H, W = a.shape[:2]
    f = max(1, int(max(H, W) / 1400))
    small = a[::f, ::f].astype(np.float64) / (2 ** img.bits - 1)
    small = (np.clip(small, 0, 1) * 255).astype(np.uint8)
    im = Image.fromarray(small)
    d = ImageDraw.Draw(im)
    h = max(2, 0.15 * det["pitch"]) / f
    for (x, y), n in zip(det["centres"], PATCH_NAMES):
        x, y = x / f, y / f
        col = (255, 0, 0) if n in det.get("tone_names", ()) else (255, 255, 0)
        d.rectangle([x - h, y - h, x + h, y + h], outline=col, width=1)
        if det["pitch"] / f > 28:
            d.text((x - h, y - h - 9), n, fill=col)
    im.save(path)


# --------------------------------------------------------------------------- correction model
def _feats(X):
    X = np.clip(X, 0, 1)
    r, g, b = X[:, 0], X[:, 1], X[:, 2]
    return np.c_[r, g, b, np.sqrt(r * g), np.sqrt(r * b), np.sqrt(g * b)]


def _full(q):
    Q = q.reshape(5, 3)
    return np.vstack([Q, 1 - Q.sum(0)])  # each output column sums to 1 -> neutral axis untouched


def _fit_matrix(X, labref, cs):
    base = np.zeros((6, 3)); base[:3] = np.eye(3)
    q0 = base[:-1].ravel()
    F = _feats(X)

    def res(q):
        return deltaE2000(cs.lab(F @ _full(q)), labref)

    sol = least_squares(res, q0, method="trf", max_nfev=300)
    return _full(sol.x)


def fit_correction(img, det, ref, dens, limits, margin_frac=0.6, log=print):
    cs = img.cs
    lin = det["lin"]
    names = PATCH_NAMES
    tone = [n for n in dens if n in names]
    idx = [names.index(n) for n in tone]
    aim = Lstar_from_density([dens[n] for n in tone])
    refL = np.array([ref[n][0] for n in tone])
    margin = margin_frac * limits["tone"]
    tL = aim + np.clip(refL - aim, -margin, margin)
    dst = np.array([cs.lab2lin(tL[i], ref[n][1], ref[n][2]) for i, n in enumerate(tone)])
    src = lin[idx]
    nodes_x, nodes_y = [], []
    for c in range(3):
        o = np.argsort(src[:, c])
        x, y = src[o, c], dst[o, c]
        keep = np.r_[True, np.diff(x) > 1e-6] & (x > 1e-6) & (x < 1 - 1e-6)  # drop duplicate / end-point inputs
        x, y = x[keep], y[keep]
        y = np.maximum.accumulate(y) + np.arange(len(y)) * 1e-7  # strictly increasing
        nodes_x.append(np.r_[0.0, x, 1.0].tolist())
        nodes_y.append(np.r_[0.0, y, 1.0].tolist())
    stage1 = [PchipInterpolator(nodes_x[c], nodes_y[c]) for c in range(3)]
    X = np.stack([np.clip(stage1[c](np.clip(lin[:, c], 0, 1)), 0, 1) for c in range(3)], 1)
    labref = np.array([ref[n] for n in names])
    P = _fit_matrix(X, labref, cs)
    model = dict(nodes_x=nodes_x, nodes_y=nodes_y, P=P.tolist())
    log("  fitted: %d tone/gray patches (target %.1f L* inside the %.1f band), %d colour patches"
        % (len(tone), margin, limits["tone"], len(names)))
    return model, tone


class Corrector:
    def __init__(self, model, cs):
        self.cs = cs
        self.stage1 = [PchipInterpolator(model["nodes_x"][c], model["nodes_y"][c]) for c in range(3)]
        self.P = np.array(model["P"])

    def luts(self, bits):
        dec = self.cs.decode_lut(bits)
        return [np.clip(self.stage1[c](np.clip(dec[c], 0, 1)), 0, 1).astype(np.float32) for c in range(3)]

    def apply_linear(self, X):
        """X (N,3) stage-1 linear values -> final linear."""
        return np.clip(_feats(X) @ self.P, 0, 1)

    def apply_codes(self, codes, bits):
        """(N,3) integer codes -> (N,3) corrected codes (float)."""
        lu = self.luts(bits)
        X = np.stack([lu[c][codes[:, c].astype(int)] for c in range(3)], 1).astype(np.float64)
        Y = self.apply_linear(X)
        return np.stack([self.cs.encode(Y[:, c], c, bits) for c in range(3)], 1)

    def apply_image(self, a, bits, strip=256, progress=None):
        lu = self.luts(bits)
        dt = np.uint8 if bits == 8 else np.uint16
        out = np.empty(a.shape, dtype=dt)
        H = a.shape[0]
        for y0 in range(0, H, strip):
            blk = a[y0:y0 + strip]
            sh = blk.shape
            X = np.stack([lu[c][blk[..., c]] for c in range(3)], -1).reshape(-1, 3)
            Y = np.clip(_feats(X.astype(np.float32)) @ self.P.astype(np.float32), 0, 1)
            for c in range(3):
                out[y0:y0 + strip, ..., c] = self.cs.encode(Y[:, c], c, bits).reshape(sh[:2]).astype(dt)
            if progress:
                progress(min(1.0, (y0 + strip) / H))
        return out


# --------------------------------------------------------------------------- evaluation / report
def evaluate(codes, img, ref, dens, limits):
    cs = img.cs
    lin = np.stack([cs.decode_lut(img.bits)[k][codes[:, k].astype(int)] for k in range(3)], 1)
    lab = cs.lab(lin)
    labref = np.array([ref[n] for n in PATCH_NAMES])
    de = deltaE2000(lab, labref)
    tone = [n for n in dens if n in PATCH_NAMES]
    rows = []
    for n in tone:
        i = PATCH_NAMES.index(n)
        aim = float(Lstar_from_density(dens[n]))
        wb = float(np.hypot(lab[i, 1] - ref[n][1], lab[i, 2] - ref[n][2]))
        rows.append(dict(patch=n, aim=aim, L=float(lab[i, 0]), tone_err=float(lab[i, 0] - aim), wb=wb))
    res = dict(
        tone_rows=rows,
        tone_worst=max(abs(r["tone_err"]) for r in rows),
        wb_worst=max(r["wb"] for r in rows),
        mean=float(de.mean()), p90=float(np.percentile(de, 90)), max=float(de.max()),
        de=de, lab=lab)
    res["pass"] = dict(tone=res["tone_worst"] <= limits["tone"], wb=res["wb_worst"] <= limits["wb"],
                       mean=res["mean"] <= limits["mean"], p90=res["p90"] <= limits["p90"])
    return res


def format_report(title, res, limits, level):
    def flag(ok, val, lim):
        if not ok:
            return "FAIL"
        return "pass" if val <= 0.88 * lim else "pass (close to limit - Open DICE may differ slightly)"
    L = []
    L.append("%s   [%d-star, limits from %s]" % (title, level, limits["source"]))
    L.append("  Tone  worst |L* - aim| : %5.2f   limit %.2f   %s" % (res["tone_worst"], limits["tone"], flag(res["pass"]["tone"], res["tone_worst"], limits["tone"])))
    L.append("  White balance worst dEab: %5.2f   limit %.2f   %s" % (res["wb_worst"], limits["wb"], flag(res["pass"]["wb"], res["wb_worst"], limits["wb"])))
    L.append("  Colour mean   dE2000    : %5.2f   limit %.2f   %s" % (res["mean"], limits["mean"], flag(res["pass"]["mean"], res["mean"], limits["mean"])))
    L.append("  Colour 90th   dE2000    : %5.2f   limit %.2f   %s" % (res["p90"], limits["p90"], flag(res["pass"]["p90"], res["p90"], limits["p90"])))
    L.append("  Colour max    dE2000    : %5.2f" % res["max"])
    worst = np.argsort(-res["de"])[:5]
    L.append("  Worst colour patches    : " + ", ".join("%s %.1f" % (PATCH_NAMES[i], res["de"][i]) for i in worst))
    return "\n".join(L)


def write_patch_csv(path, before, after, ref, dens):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["patch", "ref_L", "ref_a", "ref_b", "before_L", "before_a", "before_b", "before_dE2000",
                    "after_L", "after_a", "after_b", "after_dE2000", "tone_aim_L", "after_tone_err"])
        tone = {r["patch"]: r for r in after["tone_rows"]}
        for i, n in enumerate(PATCH_NAMES):
            t = tone.get(n)
            w.writerow([n] + ["%.2f" % v for v in ref[n]] + ["%.2f" % v for v in before["lab"][i]] +
                       ["%.2f" % before["de"][i]] + ["%.2f" % v for v in after["lab"][i]] + ["%.2f" % after["de"][i]] +
                       (["%.2f" % t["aim"], "%.2f" % t["tone_err"]] if t else ["", ""]))


# --------------------------------------------------------------------------- commands
def cmd_fit(a, log=print, progress=None):
    level = a.level
    if not os.path.isfile(a.sg):
        raise SgcalError("SG image not found: " + a.sg)
    if not os.path.isfile(a.profile):
        raise SgcalError("Profile sheet not found: " + a.profile)
    log("Reading %s" % os.path.basename(a.sg))
    img = Img(a.sg, a.colorspace)
    check_not_corrected(img, a.sg, getattr(a, "force", False))
    log("  %dx%d, %d-bit, colour space: %s" % (img.a.shape[1], img.a.shape[0], img.bits, img.cs.name))
    ref, dens = read_reference(a.profile)
    limits = read_limits(a.config, a.material, level)
    corners = [float(v) for v in a.corners.split(",")] if a.corners else None
    det = detect_and_read(img, ref, corners, log)
    det["tone_names"] = set(dens)
    before = evaluate(det["codes"], img, ref, dens, limits)
    model, tone = fit_correction(img, det, ref, dens, limits, a.margin, log)
    cor = Corrector(model, img.cs)

    outdir = a.out or os.path.join(os.path.dirname(os.path.abspath(a.sg)),
                                   os.path.splitext(os.path.basename(a.sg))[0] + "_sgcal")
    os.makedirs(outdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(a.sg))[0]
    log("Applying correction to the SG image ...")
    out = cor.apply_image(img.a, img.bits, progress=progress)
    out_path = os.path.join(outdir, stem + "_corrected.tif")
    write_tiff(out_path, out, img)

    codes2 = read_patches(out, det["centres"], det["pitch"])
    after = evaluate(codes2, img, ref, dens, limits)

    corr_path = os.path.join(outdir, "correction.json")
    with open(corr_path, "w") as fh:
        json.dump(dict(version=VERSION, created=datetime.datetime.now().isoformat(timespec="seconds"),
                       source=os.path.basename(a.sg), colorspace=img.cs.name, M=img.cs.M.tolist(),
                       level=level, limits=limits, model=model,
                       result=dict(tone_worst=after["tone_worst"], wb_worst=after["wb_worst"], mean=after["mean"],
                                   p90=after["p90"], max=after["max"])), fh, indent=1)
    write_patch_csv(os.path.join(outdir, "patches.csv"), before, after, ref, dens)
    try:
        write_overlay(img, det, os.path.join(outdir, "overlay.png"))
    except Exception as e:  # overlay is a convenience only
        log("  (overlay not written: %s)" % e)
    rep = format_report("BEFORE correction", before, limits, level) + "\n\n" + \
          format_report("AFTER correction (measured on the saved file)", after, limits, level)
    with open(os.path.join(outdir, "report.txt"), "w") as fh:
        fh.write(rep + "\n")
    log("\n" + rep)
    log("\nWrote to %s:\n  %s\n  correction.json   (apply to the rest of the session)\n  report.txt, patches.csv, overlay.png"
        % (outdir, os.path.basename(out_path)))
    if not all(after["pass"].values()):
        log("\nWARNING: at least one test is predicted to fail at %d-star." % level)
    return dict(outdir=outdir, corrected=out_path, correction=corr_path, after=after, before=before, limits=limits)


def cmd_apply(a, log=print, progress=None):
    with open(a.correction) as fh:
        c = json.load(fh)
    outdir = a.out
    done, skipped = [], []
    n = len(a.images)
    for k, p in enumerate(a.images):
        log("Applying to %s" % os.path.basename(p))
        try:
            img = Img(p, a.colorspace)
            check_not_corrected(img, p, getattr(a, "force", False))
            if not np.allclose(np.array(c["M"]), img.cs.M, atol=2e-3):
                raise SgcalError("%s uses a different colour space (%s) than the SG shot the correction was fitted on (%s)."
                                 % (os.path.basename(p), img.cs.name, c["colorspace"]))
            cor = Corrector(c["model"], img.cs)
            out = cor.apply_image(img.a, img.bits, progress=(lambda f, k=k: progress((k + f) / n)) if progress else None)
            od = outdir or os.path.dirname(os.path.abspath(p))
            os.makedirs(od, exist_ok=True)
            op = os.path.join(od, os.path.splitext(os.path.basename(p))[0] + "_corrected.tif")
            write_tiff(op, out, img)
            log("  -> " + op)
            done.append(op)
        except Exception as e:  # one bad file must not stop the batch
            msg = str(e) if isinstance(e, SgcalError) else "%s: %s" % (type(e).__name__, e)
            log("  SKIPPED %s: %s" % (os.path.basename(p), msg))
            skipped.append((p, msg))
    if skipped:
        log("\n%d file(s) skipped:" % len(skipped))
        for p, m in skipped:
            log("  %s - %s" % (os.path.basename(p), m))
    if not done and skipped:
        raise SgcalError(skipped[0][1])
    return done


def find_sg_file(paths, ref, log=print):
    """Scan the selected TIFFs and return [(match, path)] for those that contain the SG chart (best first)."""
    hits = []
    for k, p in enumerate(paths):
        try:
            img = Img(p, lazy=True)
            if img.software.startswith("sgcal"):
                continue
            det = detect_and_read(img, ref, None, lambda *x: None)
            hits.append((det["corr"], p))
            log("  SG chart found in %s (match %.3f)" % (os.path.basename(p), det["corr"]))
        except SgcalError:
            pass
        except Exception as e:
            log("  skipped %s (%s)" % (os.path.basename(p), e))
    return sorted(hits, reverse=True)


def cmd_session(a, log=print, progress=None):
    files = [p for p in a.images if os.path.isfile(p)]
    if not files:
        raise SgcalError("No image files given.")
    ref, _ = read_reference(a.profile)
    log("Looking for the ColorChecker SG among %d file(s) ..." % len(files))
    hits = find_sg_file(files, ref, log)
    if not hits:
        raise SgcalError("None of the %d selected files contains the ColorChecker SG chart. "
                         "Include the SG shot in the selection (or use 'fit' with --corners)." % len(files))
    sg = hits[0][1]
    if len(hits) > 1:
        log("  note: the chart appears in %d files; using %s. Re-run with only one SG shot to choose another."
            % (len(hits), os.path.basename(sg)))
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(sg)), "sgcal_output")
    fa = argparse.Namespace(sg=sg, profile=a.profile, config=a.config, material=a.material, level=a.level,
                            margin=a.margin, corners=None, colorspace=a.colorspace, out=out, force=False)
    res = cmd_fit(fa, log=log, progress=(lambda f: progress(0.25 * f)) if progress else None)
    others = [p for p in files if os.path.abspath(p) != os.path.abspath(sg)]
    done = [res["corrected"]]
    if others:
        log("\nApplying the correction to %d other image(s) ..." % len(others))
        aa = argparse.Namespace(correction=res["correction"], images=others, out=out,
                                colorspace=a.colorspace, force=False)
        done += cmd_apply(aa, log=log, progress=(lambda f: progress(0.25 + 0.75 * f)) if progress else None)
    res["applied"] = done
    log("\nDone: %d corrected file(s) in %s" % (len(done), out))
    return res


def build_parser():
    ap = argparse.ArgumentParser(prog="sgcal", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    f = sp.add_parser("fit", help="fit a correction from a ColorChecker SG shot")
    f.add_argument("sg", help="TIFF of the ColorChecker SG (8 or 16-bit, same exposure/light as the session)")
    f.add_argument("--profile", required=True, help="Profile_ColorCheckerSG xlsx (Lab references + densities)")
    f.add_argument("--config", help="Open DICE Config_materials xlsx (limits). Optional.")
    f.add_argument("--material", default="BoundVolumeGeneral", help="sheet name in the config xlsx")
    f.add_argument("--level", type=int, default=4, choices=[1, 2, 3, 4], help="FADGI star level to aim for (default 4)")
    f.add_argument("--margin", type=float, default=0.6,
                   help="how far (fraction of the tone band) neutrals may move toward the SG reference L*; default 0.6")
    f.add_argument("--corners", help="manual chart corners: x1,y1,...,x4,y4 = centres of A1,N1,A10,N10 (TL,TR,BL,BR)")
    f.add_argument("--colorspace", choices=list(_BUILTIN), help="override when the TIFF has no/unsupported ICC profile")
    f.add_argument("--out", help="output folder (default: <sg>_sgcal next to the image)")
    f.add_argument("--force", action="store_true", help="process a file that sgcal already corrected")
    h = sp.add_parser("session", help="give it ALL the session TIFFs: finds the SG shot, fits, corrects everything")
    h.add_argument("images", nargs="+", help="every TIFF of the session, SG shot included")
    h.add_argument("--profile", required=True)
    h.add_argument("--config")
    h.add_argument("--material", default="BoundVolumeGeneral")
    h.add_argument("--level", type=int, default=4, choices=[1, 2, 3, 4])
    h.add_argument("--margin", type=float, default=0.6)
    h.add_argument("--colorspace", choices=list(_BUILTIN))
    h.add_argument("--out", help="output folder (default: sgcal_output next to the SG shot)")
    c = sp.add_parser("chart", help="correct EACH image from the reference chart that is in its own frame (handles exposure differences)")
    c.add_argument("images", nargs="+")
    c.add_argument("--chart", default="golden_thread_DL1_1", help="chart definition name in the charts folder")
    c.add_argument("--sg-correction", help="correction.json from an SG fit: supplies the colour transform (optional)")
    c.add_argument("--out", help="output folder (default: chart_corrected next to each image)")
    c.add_argument("--compress", action="store_true", help="write compressed (deflate) TIFFs - smaller but slower")
    c.add_argument("--force", action="store_true", help="process a file that sgcal already corrected")
    g = sp.add_parser("apply", help="apply a saved correction.json to images")
    g.add_argument("correction")
    g.add_argument("images", nargs="+")
    g.add_argument("--out", help="output folder (default: next to each image)")
    g.add_argument("--colorspace", choices=list(_BUILTIN))
    g.add_argument("--force", action="store_true", help="process a file that sgcal already corrected")
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    try:
        if a.cmd == "chart":
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            import sgchart
            sgchart.cmd_chart(a)
        else:
            {"fit": cmd_fit, "apply": cmd_apply, "session": cmd_session}[a.cmd](a)
    except SgcalError as e:
        print("\nERROR: %s" % e)
        return 1
    except Exception as e:  # unreadable file, missing sheet, ...: friendly message, full trace with SGCAL_DEBUG=1
        if os.environ.get("SGCAL_DEBUG"):
            raise
        print("\nERROR: %s: %s\n(set SGCAL_DEBUG=1 for the full traceback)" % (type(e).__name__, e))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
