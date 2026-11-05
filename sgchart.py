#!/usr/bin/env python3
"""
sgchart - per-image correction from a reference chart that is IN the frame (e.g. Golden Thread).

Every image gets its own tone / gray-balance curves from the chart's neutral patches, so differences in
exposure, shutter, aperture, ISO or Phocus rendering between shots are absorbed per image. Colour comes
from the shared SG transform (optional --sg-correction) which does not depend on exposure.

The chart is located by feature matching against a reference scan of the printed target, so it can be any
size/position/rotation in the frame. Chart definitions live in the `charts` folder (JSON + reference scan).
"""
import datetime
import glob
import json
import os

import cv2
import numpy as np
from scipy.interpolate import PchipInterpolator

import sgcal
from sgcal import (ColorSpace, Corrector, Img, SgcalError, check_not_corrected, deltaE2000, write_tiff)

HERE = os.path.dirname(os.path.abspath(__file__))
CHART_DIR = os.path.join(HERE, "charts")


# --------------------------------------------------------------------------- chart definitions
def list_charts():
    return sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(CHART_DIR, "*.json")))


class Chart:
    def __init__(self, key):
        path = key if os.path.isfile(key) else os.path.join(CHART_DIR, key + ".json")
        if not os.path.isfile(path):
            raise SgcalError("Unknown chart '%s'. Available: %s" % (key, ", ".join(list_charts()) or "none"))
        d = json.load(open(path))
        self.key, self.name = key, d["name"]
        self.ref_path = os.path.join(os.path.dirname(path), d["reference_image"])
        if not os.path.isfile(self.ref_path):
            raise SgcalError("Reference scan missing: " + self.ref_path)
        self.ids = sorted(int(k) for k in d["patches"])
        self.lab = {int(k): np.array([v["L"], v["a"], v["b"]], float) for k, v in d["patches"].items()}
        self.density = {int(k): v.get("density") for k, v in d["patches"].items()}
        self.gray = list(d["gray_patches"])
        self.colour = list(d["colour_patches"])
        self.boxes = {}
        for bar in d["bars"]:
            x0, y0, x1, y1 = bar["box"]
            n = bar["last"] - bar["first"] + 1
            w = (x1 - x0 + 1) / n
            for i in range(n):
                self.boxes[bar["first"] + i] = (x0 + i * w, y0, x0 + (i + 1) * w, y1 + 1)
        self._ref = None

    def ref_gray(self):
        if self._ref is None:
            im = cv2.imread(self.ref_path, cv2.IMREAD_GRAYSCALE)
            if im is None:
                raise SgcalError("Cannot read reference scan " + self.ref_path)
            self._ref = im
        return self._ref


# --------------------------------------------------------------------------- registration
def _gray8(a, bits, maxside=3000):
    H, W = a.shape[:2]
    f = max(1, int(np.ceil(max(H, W) / maxside)))
    h, w = H // f * f, W // f * f
    g = np.empty((h // f, w // f), np.float32)
    step = f * 64
    for y in range(0, h, step):
        y1 = min(y + step, h)
        g[y // f:y1 // f] = a[y:y1, :w].astype(np.float32).reshape((y1 - y) // f, f, w // f, f, 3).mean(axis=(1, 3, 4))
    g = g / float(2 ** bits - 1) * 255
    return np.clip(g, 0, 255).astype(np.uint8), f


def register(img, chart, log=print):
    """Homography mapping reference-scan pixels -> full-resolution frame pixels."""
    g, f = _gray8(img.a, img.bits)
    clahe = cv2.createCLAHE(3.0, (8, 8))
    sift = cv2.SIFT_create(6000)
    fk, fd = sift.detectAndCompute(clahe.apply(g), None)
    if fd is None or len(fk) < 50:
        raise SgcalError("Not enough image detail to find the chart.")
    best = None
    ref = chart.ref_gray()
    for s in (1.0, 0.5, 0.25):  # chart may be much smaller in the frame than in the scan
        r = ref if s == 1.0 else cv2.resize(ref, None, fx=s, fy=s, interpolation=cv2.INTER_AREA)
        rk, rd = sift.detectAndCompute(clahe.apply(r), None)
        if rd is None or len(rk) < 20:
            continue
        m = cv2.BFMatcher().knnMatch(rd, fd, k=2)
        good = [p[0] for p in m if len(p) == 2 and p[0].distance < 0.75 * p[1].distance]
        if len(good) < 12:
            continue
        src = np.float32([rk[q.queryIdx].pt for q in good]) / s
        dst = np.float32([fk[q.trainIdx].pt for q in good])
        Hm, inl = cv2.findHomography(src, dst, cv2.RANSAC, 4.0)
        if Hm is None:
            continue
        n = int(inl.sum())
        if best is None or n > best[0]:
            best = (n, Hm)
        if n >= 300:
            break
    if best is None or best[0] < 25:
        raise SgcalError("Could not find the %s in this image (too few matching features: %s)."
                         % (chart.key, 0 if best is None else best[0]))
    n, Hm = best
    H = np.diag([f, f, 1.0]) @ Hm
    sc = float(np.sqrt(abs(np.linalg.det(H[:2, :2]))))
    rot = float(np.degrees(np.arctan2(H[1, 0], H[0, 0])))
    log("  chart found: %d feature matches, scale %.2fx, rotation %.1f deg" % (n, sc, rot))
    return H, n


def quads(chart, H, frac=0.5):
    out = {}
    for n, (x0, y0, x1, y1) in chart.boxes.items():
        cx, cy, hw, hh = (x0 + x1) / 2, (y0 + y1) / 2, (x1 - x0) * frac / 2, (y1 - y0) * frac / 2
        pts = np.float32([[cx - hw, cy - hh], [cx + hw, cy - hh], [cx + hw, cy + hh], [cx - hw, cy + hh]]).reshape(-1, 1, 2)
        out[n] = cv2.perspectiveTransform(pts, H).reshape(-1, 2)
    return out


def measure(a, chart, H, frac=0.5):
    """Median code per channel in the inner part of every patch."""
    h, w = a.shape[:2]
    codes = []
    for n, q in sorted(quads(chart, H, frac).items()):
        x0, x1 = int(np.floor(q[:, 0].min())), int(np.ceil(q[:, 0].max()))
        y0, y1 = int(np.floor(q[:, 1].min())), int(np.ceil(q[:, 1].max()))
        if x0 < 0 or y0 < 0 or x1 > w or y1 > h or x1 - x0 < 4 or y1 - y0 < 4:
            raise SgcalError("Patch %d falls outside the image - is the whole chart in the frame?" % n)
        mask = np.zeros((y1 - y0, x1 - x0), np.uint8)
        cv2.fillConvexPoly(mask, np.round(q - [x0, y0]).astype(np.int32), 1)
        codes.append(np.median(a[y0:y1, x0:x1][mask > 0].reshape(-1, 3), axis=0))
    return np.array(codes)


# --------------------------------------------------------------------------- fitting
def _to_lin(codes, cs, bits):
    dec = cs.decode_lut(bits)
    return np.stack([dec[c][codes[:, c].astype(int)] for c in range(3)], 1)


def fit_tone(lin, chart, cs, idx):
    """Per-channel monotone curves: measured linear -> the chart's neutral (L,a,b) as linear RGB."""
    nx, ny = [], []
    for c in range(3):
        x = np.array([lin[n - 1, c] for n in idx])
        y = np.array([cs.lab2lin(*chart.lab[n])[c] for n in idx])
        o = np.argsort(x)
        x, y = x[o], y[o]
        k = np.r_[True, np.diff(x) > 1e-6] & (x > 1e-6) & (x < 1 - 1e-6)
        x, y = x[k], y[k]
        y = np.maximum.accumulate(y) + np.arange(len(y)) * 1e-7
        nx.append(np.r_[0.0, x, 1.0].tolist())
        ny.append(np.r_[0.0, y, 1.0].tolist())
    return nx, ny


def _identity_P():
    P = np.zeros((6, 3))
    P[:3] = np.eye(3)
    return P


def evaluate(codes, chart, cs, bits):
    lin = _to_lin(codes, cs, bits)
    lab = cs.lab(lin)
    ref = np.array([chart.lab[n] for n in chart.ids])
    de = deltaE2000(lab, ref)
    g = [n - 1 for n in chart.gray]
    c = [n - 1 for n in chart.colour]
    return dict(
        lab=lab, de=de,
        gray_dL=float(np.abs(lab[g, 0] - ref[g, 0]).max()),
        gray_dL_signed=(lab[g, 0] - ref[g, 0]).tolist(),
        gray_wb=float(np.hypot(lab[g, 1] - ref[g, 1], lab[g, 2] - ref[g, 2]).max()),
        col_mean=float(de[c].mean()), col_p90=float(np.percentile(de[c], 90)), col_max=float(de[c].max()),
        all_mean=float(de.mean()), all_p90=float(np.percentile(de, 90)))


def exposure_estimate(lin, chart, cs):
    """Approximate exposure error in stops from the mid/light grays (log2 of reference / measured luminance)."""
    use = [n for n in chart.gray if 40 < chart.lab[n][0] < 90]
    yr = np.array([cs.lab2lin(*chart.lab[n]) @ cs.M[1] for n in use])
    ym = np.array([lin[n - 1] @ cs.M[1] for n in use])
    return float(np.median(np.log2(yr / np.maximum(ym, 1e-6))))


def warnings_for(codes, chart, bits, ev, curves_nodes, n_match):
    w = []
    mx = 2 ** bits - 1
    white = codes[chart.gray[0] - 1]
    if white.max() >= 0.985 * mx:
        w.append("Highlights are clipped or nearly clipped (white patch %s of %d): detail above that cannot be recovered."
                 % (white.astype(int).tolist(), mx))
    black = codes[chart.gray[-1] - 1]
    if black.max() < 0.02 * mx:
        w.append("Darkest patch is very close to black (%s): shadow detail will be noisy after correction." % black.astype(int).tolist())
    if abs(ev) > 0.7:
        w.append("Exposure is about %+.1f stop(s) off the chart. Large corrections amplify noise - consider re-shooting." % ev)
    if n_match < 80:
        w.append("Chart registration is weak (%d feature matches): check overlay.png." % n_match)
    # steepest part of the correction (gain)
    nx, ny = curves_nodes
    gains = [np.max(np.diff(ny[c]) / np.maximum(np.diff(nx[c]), 1e-9)) for c in range(3)]
    if max(gains) > 6:
        w.append("The tone curve has a very steep section (gain %.1fx), which can show banding/noise in shadows." % max(gains))
    return w


def _loo_gray(lin, chart, cs):
    err = []
    for n in chart.gray:
        idx = [m for m in chart.gray if m != n]
        nx, ny = fit_tone(lin, chart, cs, idx)
        cur = [PchipInterpolator(nx[c], ny[c]) for c in range(3)]
        x = np.array([[float(np.clip(cur[c](np.clip(lin[n - 1, c], 0, 1)), 0, 1)) for c in range(3)]])
        err.append(float(cs.lab(x)[0, 0] - chart.lab[n][0]))
    return err


def correct_one(path, chart, out_dir, P_sg=None, sg_cs_M=None, log=print, progress=None, compress=False, force=False):
    name = os.path.basename(path)
    log("\n== %s" % name)
    img = Img(path)
    check_not_corrected(img, path, force)
    cs = img.cs
    if sg_cs_M is not None and not np.allclose(sg_cs_M, cs.M, atol=2e-3):
        raise SgcalError("%s uses a different colour space (%s) than the SG correction was fitted in." % (name, cs.name))
    log("  %dx%d, %d-bit, %s" % (img.a.shape[1], img.a.shape[0], img.bits, cs.name))
    H, n_match = register(img, chart, log)
    codes = measure(img.a, chart, H)
    lin = _to_lin(codes, cs, img.bits)
    before = evaluate(codes, chart, cs, img.bits)
    ev = exposure_estimate(lin, chart, cs)
    nx, ny = fit_tone(lin, chart, cs, chart.gray)
    P = P_sg if P_sg is not None else _identity_P()
    model = dict(nodes_x=nx, nodes_y=ny, P=P.tolist())
    warns = warnings_for(codes, chart, img.bits, ev, (nx, ny), n_match)
    loo = _loo_gray(lin, chart, cs)
    cor = Corrector(model, cs)
    log("  applying correction ...")
    out = cor.apply_image(img.a, img.bits, progress=progress)
    codes2 = measure(out, chart, H)
    after = evaluate(codes2, chart, cs, img.bits)
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(name)[0]
    op = os.path.join(out_dir, stem + "_corrected.tif")
    write_tiff(op, out, img, compress=compress)
    with open(os.path.join(out_dir, stem + "_chartcal.json"), "w") as fh:
        json.dump(dict(version=sgcal.VERSION, created=datetime.datetime.now().isoformat(timespec="seconds"),
                       source=name, chart=chart.key, colorspace=cs.name, M=cs.M.tolist(),
                       exposure_estimate_stops=ev, colour_matrix="SG" if P_sg is not None else "none",
                       model=model, warnings=warns,
                       before={k: v for k, v in before.items() if k not in ("lab", "de")},
                       after={k: v for k, v in after.items() if k not in ("lab", "de")}), fh, indent=1)
    write_chart_csv(os.path.join(out_dir, stem + "_chart_patches.csv"), chart, before, after)
    try:
        write_overlay(img, chart, H, os.path.join(out_dir, stem + "_chart_overlay.png"))
    except Exception as e:
        log("  (overlay not written: %s)" % e)
    rep = format_report(name, chart, before, after, ev, loo, warns, n_match, P_sg is not None)
    with open(os.path.join(out_dir, stem + "_chart_report.txt"), "w") as fh:
        fh.write(rep + "\n")
    log(rep)
    log("  -> " + op)
    return dict(source=path, corrected=op, before=before, after=after, ev=ev, warnings=warns, loo=loo)


# --------------------------------------------------------------------------- outputs
def format_report(name, chart, before, after, ev, loo, warns, n_match, with_sg):
    L = ["  CHART: %s | est. exposure error %+.2f stop | colour from %s" % (chart.name.split(" (")[0], ev, "SG transform" if with_sg else "tone fit only")]
    L.append("                          before     after")
    L.append("  grays worst |dL*|     : %6.2f   %6.2f   (fitted on these patches)" % (before["gray_dL"], after["gray_dL"]))
    L.append("  gray-patch check      : leave-one-out worst |dL*| %.2f  (how well the curve interpolates between patches)" % max(abs(e) for e in loo))
    L.append("  grays worst dEab (WB) : %6.2f   %6.2f" % (before["gray_wb"], after["gray_wb"]))
    L.append("  colours mean dE2000   : %6.2f   %6.2f   <- independent of the tone fit" % (before["col_mean"], after["col_mean"]))
    L.append("  colours 90th / max    : %4.2f/%4.2f -> %4.2f/%4.2f" % (before["col_p90"], before["col_max"], after["col_p90"], after["col_max"]))
    for w in warns:
        L.append("  WARNING: " + w)
    return "\n".join(L)


def write_chart_csv(path, chart, before, after):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["patch", "ref_L", "ref_a", "ref_b", "before_L", "before_a", "before_b", "before_dE2000",
                    "after_L", "after_a", "after_b", "after_dE2000", "density"])
        for i, n in enumerate(chart.ids):
            w.writerow([n] + ["%.2f" % v for v in chart.lab[n]] + ["%.2f" % v for v in before["lab"][i]] + ["%.2f" % before["de"][i]]
                       + ["%.2f" % v for v in after["lab"][i]] + ["%.2f" % after["de"][i]] + [chart.density[n] if chart.density[n] is not None else ""])


def write_overlay(img, chart, H, path):
    a = img.a
    f = max(1, int(max(a.shape[:2]) / 1600))
    sm = (np.clip(a[::f, ::f, ::-1].astype(np.float32) / (2 ** img.bits - 1), 0, 1) * 255).astype(np.uint8).copy()
    for n, q in quads(chart, H).items():
        pts = np.round(q / f).astype(np.int32).reshape(-1, 1, 2)
        col = (0, 0, 255) if n in chart.gray else (0, 255, 255)
        cv2.polylines(sm, [pts], True, col, 2)
        cv2.putText(sm, str(n), tuple(int(v) for v in pts[0, 0] + [0, -4]), cv2.FONT_HERSHEY_SIMPLEX, 0.5, col, 1)
    cv2.imwrite(path, sm)


# --------------------------------------------------------------------------- command
def cmd_chart(a, log=print, progress=None):
    chart = Chart(a.chart)
    P_sg, M_sg = None, None
    if a.sg_correction:
        c = json.load(open(a.sg_correction))
        P_sg, M_sg = np.array(c["model"]["P"]), np.array(c["M"])
        log("Using the colour transform from %s (fitted on %s)" % (os.path.basename(a.sg_correction), c.get("source", "?")))
    done, skipped = [], []
    n = len(a.images)
    for k, p in enumerate(a.images):
        out_dir = a.out or os.path.join(os.path.dirname(os.path.abspath(p)), "chart_corrected")
        try:
            done.append(correct_one(p, chart, out_dir, P_sg, M_sg, log,
                                    (lambda f, k=k: progress((k + f) / n)) if progress else None,
                                    compress=getattr(a, "compress", False), force=getattr(a, "force", False)))
        except Exception as e:
            msg = str(e) if isinstance(e, SgcalError) else "%s: %s" % (type(e).__name__, e)
            log("  SKIPPED %s: %s" % (os.path.basename(p), msg))
            skipped.append((p, msg))
    if len(done) > 1:
        log("\nSummary (colour mean dE2000 on the 18 colour patches, before -> after):")
        for r in done:
            log("  %-34s %.2f -> %.2f   exposure %+.2f st%s" % (os.path.basename(r["source"]), r["before"]["col_mean"], r["after"]["col_mean"],
                                                               r["ev"], "   (warnings)" if r["warnings"] else ""))
    if skipped:
        log("\n%d file(s) skipped:" % len(skipped))
        for p, m in skipped:
            log("  %s - %s" % (os.path.basename(p), m))
    if not done and skipped:
        raise SgcalError(skipped[0][1])
    return done
