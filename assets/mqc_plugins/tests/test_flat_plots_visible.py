#!/usr/bin/env python

"""Regression test for the static (flat) rendering of the custom bar plot modules.

MultiQC switches from interactive plots to static images once a report holds more
samples than config.plots_flat_numseries (100 by default), which is the code path
users hit on large runs. Install the plugins first, then run this file directly or
under pytest:

    pip install -e assets/mqc_plugins/
    python assets/mqc_plugins/tests/test_flat_plots_visible.py
"""

import os
import re
import subprocess
import sys
import tempfile

# MultiQC prefixes the image div id with "mqc_" plus the pconfig id of each module.
PLOT_IDS = {
    "composition_barplots": "Composition_barplots",
    "interestgroup_comp": "interestgroup_composition_bargraph",
}
SAMPLES = ["SAMPLE_{:02d}".format(i) for i in range(4)]


def _write_csv(path, index_name, rows):
    with open(path, "w") as fh:
        fh.write(",".join([index_name] + SAMPLES) + "\n")
        for name, values in rows:
            fh.write(",".join([name] + ["{:.3f}".format(v) for v in values]) + "\n")


def _write_inputs(indir):
    for level in range(1, 8):
        rows = [
            ("taxon_{}_{}".format(level, i), [40.0, 30.0, 20.0, 10.0])
            for i in range(1, 3)
        ]
        rows.append(("taxon_{}_3".format(level), [20.0, 40.0, 60.0, 80.0]))
        _write_csv(os.path.join(indir, "level-{}.csv".format(level)), "#OTU ID", rows)
    for group in ("GroupA", "GroupB"):
        rows = [("{}_taxon_{}".format(group, i), [1.0, 2.0, 3.0, 4.0]) for i in range(2)]
        _write_csv(
            os.path.join(indir, "{}_groupinterest_comp.csv".format(group)), "taxon", rows
        )


def _report_html(workdir):
    indir = os.path.join(workdir, "input")
    outdir = os.path.join(workdir, "output")
    os.makedirs(indir)
    _write_inputs(indir)
    cmd = [sys.executable, "-m", "multiqc", "--force", "--flat", "-o", outdir]
    for module in PLOT_IDS:
        cmd += ["-m", module]
    subprocess.check_call(cmd + [indir], stdout=subprocess.DEVNULL)
    with open(os.path.join(outdir, "multiqc_report.html"), encoding="utf-8") as fh:
        return fh.read()


def test_flat_plots_are_not_hidden():
    with tempfile.TemporaryDirectory() as workdir:
        html = _report_html(workdir)
    for module, plot_id in PLOT_IDS.items():
        images = re.findall(
            r'<div class="mqc_mplplot" id="mqc_{}[^"]*"([^>]*)>'.format(plot_id), html
        )
        assert images, "{}: no static image was rendered".format(module)
        visible = [attrs for attrs in images if "display:none" not in attrs]
        assert visible, "{}: all {} static images were rendered hidden".format(
            module, len(images)
        )


if __name__ == "__main__":
    test_flat_plots_are_not_hidden()
    print("OK: static bar plot images are rendered visible")
