#!/usr/bin/env python3
"""SG Calibrator - simple window around sgcal.py (fit from a ColorChecker SG shot, apply to a session)."""
import argparse
import json
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import sgcal
try:
    import sgchart
except Exception:  # OpenCV missing: chart mode disabled, the rest still works
    sgchart = None

SETTINGS = os.path.join(os.path.expanduser("~"), ".sgcal_gui.json")


def load_settings():
    try:
        with open(SETTINGS) as fh:
            return json.load(fh)
    except Exception:
        return {}


def save_settings(d):
    try:
        with open(SETTINGS, "w") as fh:
            json.dump(d, fh, indent=1)
    except Exception:
        pass


def open_path(p):
    try:
        if sys.platform.startswith("win"):
            os.startfile(p)  # noqa
        elif sys.platform == "darwin":
            subprocess.Popen(["open", p])
        else:
            subprocess.Popen(["xdg-open", p])
    except Exception as e:
        messagebox.showerror("Open", str(e))


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("SG Calibrator - FADGI tone / white balance / colour")
        self.geometry("860x700")
        self.minsize(760, 600)
        self.q = queue.Queue()
        self.busy = False
        self.last_out = None
        s = load_settings()
        self.v_sg = tk.StringVar(value=s.get("sg", ""))
        self.v_profile = tk.StringVar(value=s.get("profile", ""))
        self.v_config = tk.StringVar(value=s.get("config", ""))
        self.v_material = tk.StringVar(value=s.get("material", "BoundVolumeGeneral"))
        self.v_level = tk.StringVar(value=s.get("level", "4"))
        self.v_corr = tk.StringVar(value=s.get("corr", ""))
        self.v_outdir = tk.StringVar(value=s.get("outdir", ""))
        self._build()
        self._refresh_materials()
        self.after(100, self._poll)

    # ------------------------------------------------------------------ layout
    def _row(self, parent, r, label, var, browse):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(parent, textvariable=var).grid(row=r, column=1, sticky="ew", pady=3)
        ttk.Button(parent, text="Browse...", command=browse).grid(row=r, column=2, padx=(8, 0), pady=3)

    def _build(self):
        pad = dict(padx=12, pady=6)
        f1 = ttk.LabelFrame(self, text="1. Fit a correction from the ColorChecker SG shot (do this once per session)")
        f1.pack(fill="x", **pad)
        f1.columnconfigure(1, weight=1)
        self._row(f1, 0, "SG shot (TIFF):", self.v_sg, lambda: self._pick(self.v_sg, "SG TIFF", [("TIFF", "*.tif *.tiff")]))
        self._row(f1, 1, "Profile sheet (xlsx):", self.v_profile, lambda: self._pick(self.v_profile, "Profile_ColorCheckerSG xlsx", [("Excel", "*.xlsx")]))
        self._row(f1, 2, "Config sheet (xlsx, optional):", self.v_config, self._pick_config)
        opt = ttk.Frame(f1)
        opt.grid(row=3, column=0, columnspan=3, sticky="w", pady=3)
        ttk.Label(opt, text="Material:").pack(side="left")
        self.cb_mat = ttk.Combobox(opt, textvariable=self.v_material, width=34)
        self.cb_mat.pack(side="left", padx=(6, 18))
        ttk.Label(opt, text="Aim for:").pack(side="left")
        ttk.Combobox(opt, textvariable=self.v_level, values=["4", "3", "2", "1"], width=3, state="readonly").pack(side="left", padx=6)
        ttk.Label(opt, text="star").pack(side="left")
        bf = ttk.Frame(f1)
        bf.grid(row=4, column=0, columnspan=3, sticky="e", pady=(6, 4))
        self.btn_session = ttk.Button(bf, text="One step: select ALL session images (finds the SG shot itself)...", command=self.on_session)
        self.btn_session.pack(side="left", padx=(0, 8))
        self.btn_fit = ttk.Button(bf, text="Fit only (from the SG shot above)", command=self.on_fit)
        self.btn_fit.pack(side="left")

        f2 = ttk.LabelFrame(self, text="2. Apply the correction to the session's images (same light, lens, exposure and colour space)")
        f2.pack(fill="x", **pad)
        f2.columnconfigure(1, weight=1)
        self._row(f2, 0, "Correction file:", self.v_corr, lambda: self._pick(self.v_corr, "correction.json", [("JSON", "*.json")]))
        self._row(f2, 1, "Output folder (blank = next to each image):", self.v_outdir, self._pick_dir)
        self.btn_apply = ttk.Button(f2, text="Choose images and apply...", command=self.on_apply)
        self.btn_apply.grid(row=2, column=0, columnspan=3, sticky="e", pady=(6, 4))

        fc = ttk.LabelFrame(self, text="3. Chart in every frame - corrects each image from its own chart (handles different exposures)")
        fc.pack(fill="x", **pad)
        row = ttk.Frame(fc)
        row.pack(fill="x", pady=3)
        ttk.Label(row, text="Chart:").pack(side="left")
        charts = sgchart.list_charts() if sgchart else []
        self.v_chart = tk.StringVar(value=(charts[0] if charts else "(OpenCV not installed)"))
        ttk.Combobox(row, textvariable=self.v_chart, values=charts, width=28, state="readonly").pack(side="left", padx=(6, 14))
        self.v_usesg = tk.BooleanVar(value=True)
        ttk.Checkbutton(row, text="use colour transform from the correction file above", variable=self.v_usesg).pack(side="left", padx=(0, 14))
        self.v_compress = tk.BooleanVar(value=False)
        ttk.Checkbutton(row, text="compress output TIFFs", variable=self.v_compress).pack(side="left")
        self.btn_chart = ttk.Button(fc, text="Select images and correct each from its own chart...", command=self.on_chart,
                                    state=("normal" if sgchart else "disabled"))
        self.btn_chart.pack(anchor="e", pady=(2, 4))

        f3 = ttk.Frame(self)
        f3.pack(fill="both", expand=True, **pad)
        self.result = tk.Label(f3, text="", anchor="w", justify="left", font=("Segoe UI", 11, "bold"))
        self.result.pack(fill="x")
        self.prog = ttk.Progressbar(f3, mode="determinate", maximum=1.0)
        self.prog.pack(fill="x", pady=4)
        self.log = tk.Text(f3, height=14, wrap="word", font=("Consolas", 9))
        sb = ttk.Scrollbar(f3, command=self.log.yview)
        self.log.configure(yscrollcommand=sb.set)
        self.log.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        fb = ttk.Frame(self)
        fb.pack(fill="x", padx=12, pady=(0, 10))
        self.btn_folder = ttk.Button(fb, text="Open output folder", command=lambda: self.last_out and open_path(self.last_out), state="disabled")
        self.btn_folder.pack(side="left")
        self.btn_overlay = ttk.Button(fb, text="Show where patches were read", command=lambda: self.last_out and open_path(os.path.join(self.last_out, "overlay.png")), state="disabled")
        self.btn_overlay.pack(side="left", padx=8)

    # ------------------------------------------------------------------ pickers
    def _pick(self, var, title, types):
        init = os.path.dirname(var.get()) if var.get() else ""
        p = filedialog.askopenfilename(title=title, filetypes=types, initialdir=init or None)
        if p:
            var.set(p)

    def _pick_config(self):
        self._pick(self.v_config, "Open DICE Config_materials xlsx", [("Excel", "*.xlsx")])
        self._refresh_materials()

    def _pick_dir(self):
        p = filedialog.askdirectory(title="Output folder")
        if p:
            self.v_outdir.set(p)

    def _refresh_materials(self):
        mats = ["BoundVolumeGeneral"]
        if self.v_config.get() and os.path.isfile(self.v_config.get()):
            try:
                mats = sgcal.list_materials(self.v_config.get())
            except Exception:
                pass
        self.cb_mat["values"] = mats
        if self.v_material.get() not in mats:
            self.v_material.set(mats[0])

    # ------------------------------------------------------------------ worker plumbing
    def _log(self, msg):
        self.q.put(("log", msg))

    def _poll(self):
        try:
            while True:
                kind, val = self.q.get_nowait()
                if kind == "log":
                    self.log.insert("end", val + "\n")
                    self.log.see("end")
                elif kind == "progress":
                    self.prog["value"] = val
                elif kind == "done":
                    self._finish(val)
        except queue.Empty:
            pass
        self.after(100, self._poll)

    def _set_busy(self, b):
        self.busy = b
        st = "disabled" if b else "normal"
        self.btn_fit.configure(state=st)
        self.btn_session.configure(state=st)
        self.btn_apply.configure(state=st)
        self.btn_chart.configure(state=("disabled" if (b or not sgchart) else "normal"))

    def _run(self, fn):
        if self.busy:
            return
        self._set_busy(True)
        self.prog["value"] = 0
        self.log.delete("1.0", "end")
        self.result.configure(text="Working...", fg="#444444")

        def work():
            try:
                res = fn()
                self.q.put(("done", ("ok", res)))
            except sgcal.SgcalError as e:
                self.q.put(("done", ("err", str(e))))
            except Exception as e:
                self.q.put(("done", ("err", "%s: %s" % (type(e).__name__, e))))
        threading.Thread(target=work, daemon=True).start()

    def _finish(self, val):
        self._set_busy(False)
        self.prog["value"] = 1.0
        status, res = val
        if status == "err":
            self.result.configure(text="Error - see log", fg="#b00020")
            self._log("\nERROR: " + res)
            return
        if isinstance(res, dict) and "after" in res:
            ok = all(res["after"]["pass"].values())
            lim = res["limits"]
            txt = ("Predicted PASS at %d-star" if ok else "Predicted FAIL at %d-star - see report") % int(self.v_level.get())
            if "applied" in res:
                txt += "  |  %d file(s) corrected" % len(res["applied"])
            self.result.configure(text=txt, fg="#0a7d28" if ok else "#b00020")
            self.v_corr.set(res["correction"])
            self.last_out = res["outdir"]
            self.btn_folder.configure(state="normal")
            self.btn_overlay.configure(state="normal")
        elif res and isinstance(res[0], dict):   # chart mode
            after = [r["after"]["col_mean"] for r in res]
            warn = sum(1 for r in res if r["warnings"])
            txt = "%d image(s) corrected from their own chart - colour dE2000 after: mean %.2f, worst image %.2f" % (len(res), sum(after) / len(after), max(after))
            if warn:
                txt += "  |  %d with warnings (see log)" % warn
            self.result.configure(text=txt, fg="#b36b00" if warn else "#0a7d28")
            self.last_out = os.path.dirname(res[0]["corrected"])
            self.btn_folder.configure(state="normal")
        else:
            self.result.configure(text="Done - %d image(s) written" % len(res), fg="#0a7d28")
            if res:
                self.last_out = os.path.dirname(res[0])
                self.btn_folder.configure(state="normal")
        self._save()

    def _save(self):
        save_settings(dict(sg=self.v_sg.get(), profile=self.v_profile.get(), config=self.v_config.get(),
                           material=self.v_material.get(), level=self.v_level.get(), corr=self.v_corr.get(),
                           outdir=self.v_outdir.get()))

    # ------------------------------------------------------------------ actions
    def on_fit(self):
        if not self.v_sg.get() or not self.v_profile.get():
            messagebox.showinfo("Fit", "Choose the SG TIFF and the Profile_ColorCheckerSG xlsx first.")
            return
        ns = argparse.Namespace(sg=self.v_sg.get(), profile=self.v_profile.get(), config=self.v_config.get() or None,
                                material=self.v_material.get(), level=int(self.v_level.get()), margin=0.6,
                                corners=None, colorspace=None, out=None, force=False)
        self._save()
        self._run(lambda: sgcal.cmd_fit(ns, log=self._log, progress=lambda f: self.q.put(("progress", f))))

    def on_session(self):
        if not self.v_profile.get():
            messagebox.showinfo("Session", "Choose the Profile_ColorCheckerSG xlsx first.")
            return
        init = os.path.dirname(self.v_sg.get()) if self.v_sg.get() else None
        files = filedialog.askopenfilenames(title="Select ALL TIFFs of the session (include the SG shot)",
                                            filetypes=[("TIFF", "*.tif *.tiff")], initialdir=init)
        if not files:
            return
        ns = argparse.Namespace(images=list(files), profile=self.v_profile.get(), config=self.v_config.get() or None,
                                material=self.v_material.get(), level=int(self.v_level.get()), margin=0.6,
                                colorspace=None, out=self.v_outdir.get() or None)
        self._save()
        self._run(lambda: sgcal.cmd_session(ns, log=self._log, progress=lambda f: self.q.put(("progress", f))))

    def on_chart(self):
        files = filedialog.askopenfilenames(title="Images with the reference chart in the frame", filetypes=[("TIFF", "*.tif *.tiff")])
        if not files:
            return
        corr = self.v_corr.get() if (self.v_usesg.get() and self.v_corr.get() and os.path.isfile(self.v_corr.get())) else None
        ns = argparse.Namespace(images=list(files), chart=self.v_chart.get(), sg_correction=corr, out=self.v_outdir.get() or None,
                                compress=self.v_compress.get(), force=False)
        self._save()
        self._run(lambda: sgchart.cmd_chart(ns, log=self._log, progress=lambda f: self.q.put(("progress", f))))

    def on_apply(self):
        if not self.v_corr.get() or not os.path.isfile(self.v_corr.get()):
            messagebox.showinfo("Apply", "Choose a correction.json first (it is filled in automatically after a fit).")
            return
        files = filedialog.askopenfilenames(title="Images to correct", filetypes=[("TIFF", "*.tif *.tiff")])
        if not files:
            return
        ns = argparse.Namespace(correction=self.v_corr.get(), images=list(files), out=self.v_outdir.get() or None,
                                colorspace=None, force=False)
        self._save()
        self._run(lambda: sgcal.cmd_apply(ns, log=self._log, progress=lambda f: self.q.put(("progress", f))))


def main():
    App().mainloop()


if __name__ == "__main__":
    main()
