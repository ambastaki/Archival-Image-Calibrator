# Chart definitions

Chart mode (`sgcal.py chart ...`) needs one folder entry per target: a JSON file plus a reference scan of the target.
The files for the Golden Thread / Object-Level Target are **not included in this repository** because the target's
reference values and artwork belong to its makers. If you own the target and may use its values, create your own:

    charts/<key>.json            # see format below
    charts/<key>_ref.jpg         # a clean scan or photo of the target; box coordinates are in its pixels

Format:

    {
     "name": "My target",
     "reference_image": "<key>_ref.jpg",
     "bars": [ {"first": 1, "last": 9, "box": [x0, y0, x1, y1]} ],   // equal-width patches left to right inside box
     "gray_patches": [10, 11, 12],
     "colour_patches": [1, 2, 3],
     "patches": { "1": {"L": 38.9, "a": 13.3, "b": 14.3, "density": null} }   // D50 Lab, density for grays
    }
