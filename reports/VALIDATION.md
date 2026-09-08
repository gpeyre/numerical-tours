# Python collection validation

Validated on 8 September 2026 with Python 3.12.14 on macOS, using CPU execution.

## Results

- **62/62 main Python notebooks passed**, each in a fresh Jupyter kernel.
- **1,892 code cells executed**; no cells were skipped and notebook exceptions were not ignored.
- **71 tests passed**: numerical invariants, notebook structure, complete inline implementations, and 5–10 bibliography entries per notebook.
- Undefined-name and invalid-string-escape checks passed for every main notebook.
- Every saved execution report matches the SHA-256 of its current source notebook.
- Editorial regeneration was checked on a temporary copy of all 62 notebooks.

The machine-readable results are in [execution.json](execution.json). The exact installed package versions are in [environment.txt](environment.txt). Reproduce execution and export using the instructions in [maintenance/README.md](../maintenance/README.md).

## Website

The website export contains 58 distinct Python tours with complete executed reading pages, downloadable notebooks, previews extracted from their computed figures, and Colab links. Four historical notebook variants remain in the content repository but are omitted from the main catalogue.

Search tests cover multiword queries, accents, topic filters, empty results, reset, URL state, and literal input. The Jekyll build and generated-link checks are recorded in the website repository's `reports/VALIDATION.md`.

Source notebooks retain clean outputs; the downloadable website notebooks contain the executed outputs. Colab links target the published `master` branch, so the content changes must be published before the website changes.

## Scope and limits

The 128 unfinished conversions under `python/todo/`, editor checkpoints, and historical notebooks stored alongside solution files were excluded from execution. They are drafts, not passing members of this collection. The main notebooks no longer depend on hidden exercise-solution files.

Light editorial corrections were made to 16 legacy notebooks (8 MATLAB, 4 R, 4 Julia); these languages were not execution-tested. The website keeps their existing routes under a secondary language archive and repairs the old MATLAB pages' MathJax endpoint.

Validation checks code execution and selected mathematical invariants; it is not an independent proof of every mathematical claim. Browser visual testing and remote Colab execution were not performed. The deep texture tour requires an initial pretrained-weight download of approximately 550 MB; all experiments use practical CPU defaults.

Publish the validated content commit on `master` before the corresponding `gh-pages` website commit so that Colab and website downloads refer to the same revision of the teaching material.
