# SG Calibrator

Fits a tone, white-balance and colour correction from a photo of the **ColorChecker SG** and checks it against
**FADGI** limits the way Open DICE does. Built for Hasselblad / Phocus capture workflows, but works on any TIFF with an embedded ICC profile.

## Set up (once)
1. Install Python 3 from python.org (tick *Add python.exe to PATH* and keep *tcl/tk* selected).
2. Double-click `install.bat`.
3. Double-click `run_gui.bat`.

## When exposure differs between shots: chart in every frame (recommended for the Golden Thread workflow)
Each image is corrected from the reference chart **inside that image**, so different shutter / aperture / ISO / Phocus
rendering between shots is absorbed per image.
1. Put the chart in the frame (the whole chart must be visible; any size, position or rotation). Chart definitions live in
   the `charts` folder - see `charts/README.md` for how to add your target (the Golden Thread files are not bundled).
2. In the window (section 3) pick the chart, tick *use colour transform from the correction file above* if you have fitted an SG
   correction (optional - colour is good from the chart's grays alone, the SG matrix improves it slightly), press
   **Select images and correct each from its own chart...** and select the TIFFs.
3. Per image you get `*_corrected.tif`, a report, per-patch CSV, `*_chartcal.json` (the exact adjustment applied - keep it for
   the preservation record) and an overlay showing where each patch was read.
   Command line: `py -3 sgcal.py chart img1.tif img2.tif --sg-correction SG_sgcal\correction.json`
Reading the report: the 12 gray patches are what the tone curve is *fitted on*, so they match by construction. The honest
checks are the **colour patches** (not used in the tone fit) and the **leave-one-out** gray figure. Warnings flag clipped
highlights, near-black shadows, exposure more than ~0.7 stop off, and weak chart detection.
Limits: clipped highlights cannot be recovered; large corrections amplify noise; tests so far used simulated exposure changes plus a few real frames; validate on your own bracketed shots.

## Each session - fastest way (one step)
1. Export all of the session's TIFFs from Phocus into one folder (SG shot included, same light/lens/exposure/LCC).
2. In the window choose the `Profile_ColorCheckerSG` xlsx (and the Open DICE config sheet, material and star level),
   then press **One step: select ALL session images...** and select every TIFF.
3. It looks through the files, finds the one containing the ColorChecker SG, fits the correction from it, and corrects
   every selected image. Everything lands in a `sgcal_output` folder next to the SG shot: `*_corrected.tif` for each image,
   plus `correction.json`, `report.txt`, `patches.csv` and `overlay.png`.
   Files already corrected by this tool, or unreadable files, are skipped and listed - the rest carry on.
   If the chart appears in more than one file, the best match is used and the others are listed in the log.
   Command line: `py -3 sgcal.py session *.tif --profile Profile_ColorCheckerSG.xlsx --config Config_materials.xlsx`

## Each session - step by step
1. Shoot the SG with the **same light, lens, aperture, exposure and LCC** as the session. Export a TIFF from Phocus
   (8 or 16-bit, Adobe RGB / sRGB / ProPhoto, **no extra tone curve or sharpening**). Cropping is optional - the chart
   is found automatically, even rotated or upside down.
2. In the window: choose the SG TIFF, `Profile_ColorCheckerSG` xlsx and (optionally) the Open DICE `Config_materials` xlsx,
   pick the material and star level, press **Fit correction**.
3. Read the result. *Predicted PASS* means tone, white balance, mean and 90th-percentile colour are inside the limits.
   `overlay.png` shows exactly where each patch was read (red = the 12 tone patches).
4. Press **Choose images and apply...** and pick the session's TIFFs. Corrected copies are saved as `*_corrected.tif`
   (same bit depth, ICC profile and resolution; originals are never overwritten).
5. Run the corrected SG TIFF through Open DICE to confirm.

Output folder (`<sg name>_sgcal`): corrected SG TIFF, `correction.json` (reuse it for the session), `report.txt`,
`patches.csv` (per-patch before/after Lab and dE2000), `overlay.png`.

## Command line
    py -3 sgcal.py fit  SG.tif --profile Profile_ColorCheckerSG.xlsx --config Config_materials.xlsx --material BoundVolumeGeneral --level 4
    py -3 sgcal.py apply SG_sgcal\correction.json img1.tif img2.tif --out corrected
If the chart is not found: `--corners x1,y1,x2,y2,x3,y3,x4,y4` = pixel centres of the four white corner patches
A1, N1, A10, N10 (top-left, top-right, bottom-left, bottom-right). `--margin` (default 0.6) sets how far neutrals may
move toward the SG's own reference L* inside the tone band: raise it for better colour, lower it for more tone headroom.

## What it does
* Reads the colour matrix and curve from the TIFF's embedded ICC profile and converts to D50 Lab, like Open DICE.
* Tone aim per patch = L* of density (10^-D), from the density column of the profile sheet.
* Per-channel curves put the 12 gray patches on the aim (and on their reference a*b*), then a 3x3-class root-polynomial
  transform that leaves the gray axis alone is fitted on all 140 patches by minimising dE2000.
* Applies in strips, so a 100 MP 16-bit file takes roughly half a minute and ~1.5 GB of RAM.

## Limits - read these
* It corrects the *colour and tone* tests. Sharpness, noise, uniformity and registration come from the capture
  (focus, light, LCC) and are not changed.
* The correction is valid for the lighting, lens and exposure of the SG shot. Re-shoot the SG when any of them changes.
* Results are fitted on the SG and judged on the SG; Open DICE does the same, but the correction is not proof of
  accuracy on other targets. Keep the original capture settings close to correct - the closer Phocus is to right,
  the gentler (and safer) the correction.
* Predictions matched Open DICE to about 0.02 on the mean and 90th-percentile dE of a test file; keep an eye on any
  result flagged "close to limit".
* The tool refuses files it already corrected (use `--force` to override).


## Status and disclaimer
Independent project, not affiliated with or endorsed by FADGI, Open DICE, Hasselblad, X-Rite or Image Science Associates. Predicted pass/fail mirrors Open DICE's calculations closely but is **not** a certification: always confirm with Open DICE. Corrections change pixel values: keep originals and the `*_chartcal.json` / `correction.json` records alongside the corrected files.
