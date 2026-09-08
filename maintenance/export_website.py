"""Export executed notebooks, scientific previews, and searchable website metadata."""

from __future__ import annotations
import argparse
import base64
import hashlib
import html
import io
import json
from pathlib import Path
import re
import shutil
import nbformat
from nbconvert import HTMLExporter
from PIL import Image, ImageOps
from catalog_data import catalog

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", type=Path, required=True)
    parser.add_argument("--executed", type=Path, required=True)
    parser.add_argument(
        "--allow-stale",
        action="store_true",
        help="Development preview only; never use for release",
    )
    args = parser.parse_args()
    site = args.site.resolve()
    exporter = HTMLExporter(
        template_name="classic", exclude_input_prompt=True, exclude_output_prompt=True
    )
    exporter.mathjax_url = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.9/MathJax.js?config=TeX-AMS_CHTML-full"
    data = []
    for slug, info in catalog().items():
        if "canonical" in info:
            continue
        source = ROOT / "python" / f"{slug}.ipynb"
        executed = args.executed / source.name
        report = json.loads((args.executed / f"{slug}.json").read_text())
        if report["status"] != "passed":
            raise RuntimeError(f"{slug} did not pass execution")
        if (
            not args.allow_stale
            and report["sha256"] != hashlib.sha256(source.read_bytes()).hexdigest()
        ):
            raise RuntimeError(f"{slug} has changed since execution")
        notebook = nbformat.read(executed, as_version=4)
        # Preview work may use fresh prose while a final execution is in progress.
        if args.allow_stale:
            prose = {
                c.id: c.source
                for c in nbformat.read(source, as_version=4).cells
                if c.cell_type == "markdown"
            }
            for c in notebook.cells:
                if c.cell_type == "markdown" and c.id in prose:
                    c.source = prose[c.id]
        images = []
        for ci, c in enumerate(notebook.cells):
            for output in c.get("outputs", []):
                png = output.get("data", {}).get("image/png")
                if png:
                    image = Image.open(io.BytesIO(base64.b64decode(png))).convert("RGB")
                    # Exclude empty figures and prefer substantive, landscape results.
                    if image.width >= 280 and image.height >= 150:
                        images.append((ci, image))
        if not images:
            raise RuntimeError(f"No scientific figure available for {slug}")
        # Choose the last substantive result unless a teaching-specific preview is supplied.
        preferred = {
            "ml_1_pca_nn": "cell-026",
            "shapes_7_isomap": "cell-012",
            "ml_8_deep_texture_synthesis": "cell-060",
            "ml_12_diffusion_models": "cell-064",
            "ml_12_diffusion_model_jax": "cell-060",
            "ml_11_conformal_prediction": "cell-021",
        }
        chosen = next(
            (
                im
                for ci, im in images
                if notebook.cells[ci].get("id") == preferred.get(slug)
            ),
            images[-1][1],
        )
        preview = ImageOps.pad(
            chosen, (720, 420), color="#f5f8fb", method=Image.Resampling.LANCZOS
        )
        image_path = site / "assets" / "tours" / f"{slug}.webp"
        image_path.parent.mkdir(parents=True, exist_ok=True)
        preview.save(image_path, "WEBP", quality=88)
        notebook.metadata["title"] = info["title"]
        body, _ = exporter.from_notebook_node(notebook)
        title = html.escape(info["title"])
        colab = f"https://colab.research.google.com/github/gpeyre/numerical-tours/blob/master/python/{slug}.ipynb"
        css = """<link rel="stylesheet" href="/assets/css/tours.css"><style>
body{background:#fff;padding:0;font-family:system-ui,-apple-system,sans-serif;font-size:16px;color:#192d43}#notebook{padding:25px 20px 70px}#notebook-container{max-width:1000px;width:100%;padding:20px 35px;box-shadow:none}div.text_cell_render{font-family:inherit;font-size:16px;line-height:1.75}div.text_cell_render h1{font-family:Georgia,serif;font-size:2.7rem;font-weight:400;line-height:1.18;color:#12243a;margin-top:10px}div.text_cell_render h2{font-family:Georgia,serif;font-size:1.85rem;font-weight:400;border-top:1px solid #dce4ec;padding-top:25px;margin-top:30px}div.text_cell_render h3{font-size:1.25rem}div.input_area{border:1px solid #dce4ec;background:#f5f8fb;border-radius:5px;padding:10px}div.input_area pre{font-size:13px;line-height:1.55}div.output_area pre{font-size:13px}div.output_subarea{max-width:100%}div.output_png img{max-width:100%;height:auto}a{color:#185bce}div.cell{padding:8px 0}div.output_stderr{background:#fff6e0;border-left:3px solid #bc8100}.celltoolbar{display:none}@media(max-width:650px){#notebook{padding:15px 8px}#notebook-container{padding:10px 8px}div.text_cell_render h1{font-size:2rem}div.text_cell_render{padding:8px 5px}div.input_area pre{font-size:12px}.rendered_html table{display:block;overflow-x:auto}}
</style>"""
        body = body.replace(
            "</head>",
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            + css
            + "</head>",
        )
        toolbar = f'''<a class="skip-link" href="#main">Skip to notebook</a><header class="site-header"><div class="header-inner"><a class="brand" href="/"><span class="brand-mark" aria-hidden="true">∿</span>Numerical Tours</a><nav aria-label="Main navigation"><a href="/python/">Python tours</a><a href="/installation_python/">Getting started</a><a href="/archive/">Other languages</a></nav></div></header><div class="notebook-toolbar"><a href="/python/">← All Python tours</a><div class="notebook-actions"><a href="/downloads/{slug}.ipynb" download>Download notebook</a><a class="colab-link" href="{colab}"><span class="colab-icon" aria-hidden="true">co</span> Open in Colab</a></div></div><main id="main" aria-label="{title}">'''
        body = re.sub(r"(<body[^>]*>)", lambda m: m[1] + toolbar, body, count=1)
        body = body.replace(
            "</body>",
            '</main><footer class="site-footer"><a href="/python/">Explore another tour →</a><a href="/about/">Gabriel Peyré and contributors</a></footer></body>',
        )
        destination = site / "python" / slug
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "index.html").write_text(
            "\n".join(line.rstrip() for line in body.splitlines()) + "\n"
        )
        (site / "downloads").mkdir(exist_ok=True)
        shutil.copy2(executed, site / "downloads" / source.name)
        data.append(
            dict(
                slug=slug,
                title=info["title"],
                category=info["category"],
                summary=info["intro"].split(". ")[0] + ".",
                search=" ".join(
                    [
                        info["title"],
                        info["category"],
                        info["intro"],
                        " ".join(info["keywords"]),
                    ]
                ),
                image=f"/assets/tours/{slug}.webp",
                image_alt=f"Computed results from {info['title']}",
                url=f"/python/{slug}/",
                colab=colab,
            )
        )
    (site / "_data" / "python_catalog.json").write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n"
    )
    print(f"Exported {len(data)} reading pages, notebooks, and figure previews.")


if __name__ == "__main__":
    main()
