# Maintaining the Python tours and website

## Install the validation environment

Use Python 3.12. From the content repository root:

```sh
python -m pip install -r maintenance/requirements.txt
```

`python/requirements*.txt` describe supported dependency ranges. `reports/environment.txt` records the exact versions used for the saved execution report.

## Execute the collection

```sh
python maintenance/validate_notebooks.py --output reports/executed --jobs 3
```

Every notebook runs in a fresh Jupyter kernel, in its own directory. No cells are skipped, exceptions stop execution, and the runner returns a failing exit code if any tour fails. It records elapsed times, code-cell counts, and the SHA-256 of each source notebook. Executed copies remain in the output directory; sources are not modified by the runner. Omit notebook names to run all 62 main notebooks, or supply names to rerun a subset.

The 128 unfinished conversions under `python/todo/`, checkpoint files, and notebooks stored with historical solution files are outside the maintained collection. They are not counted as passing tours.

Run the shared numerical invariants and structural checks with:

```sh
python -m pytest tests
ruff check python/*.ipynb --select F821,W605
```

The numerical checks cover adjoints, transform inversion and energy, interpolation boundaries, convolution, and linear-program feasibility. End-to-end execution alone is not a proof that every statement or experiment is mathematically correct.

## Editorial metadata

`catalog_data.py` contains each tour’s title, introduction, topic, keywords, and reading list. The notebooks remain the source of mathematical exposition and implementations. `apply_editorial.py` applies the shared introductory material and portable setup. After editing notebook content or running the editorial tool, rerun formatting and execution before exporting.

## Export to the website branch

With the `gh-pages` branch checked out alongside this repository:

```sh
python maintenance/export_website.py \
  --site ../numerical-tours-site \
  --executed reports/executed
```

The exporter requires a passing report matching the exact source hash. It produces 58 HTML reading pages, downloadable executed notebooks, figure previews extracted from actual outputs, and a searchable catalogue. Four historical variants retain their notebook URLs but are omitted from the main catalogue.

Build the website using its existing Jekyll workflow. The Colab buttons target `master`, so publish the tested content commit before publishing the corresponding website changes. A local Colab button cannot reference unpublished working-tree edits.

Do not commit editor checkpoints, caches, environments, temporary exports, or generated Jekyll `_site` output. Keep the compact execution report and exact environment record under `reports/`.
