# Python Numerical Tours

The main collection contains **62 notebooks**, including four historical variants retained for existing links. The website presents **58 distinct tours** with introductions, complete worked examples, figures, and topic-specific references.

## Start in Colab

Open a tour on [the website](https://www.numerical-tours.com/python/) and select **Open in Colab**. Run the setup cell, then run all cells in order. The setup locates or downloads the companion toolbox and data and installs missing dependencies. There are no separate solution files to run.

## Run locally

Use Python 3.12 and a separate environment:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r python/requirements.txt jupyterlab
jupyter lab python
```

On Windows, activate with `.venv\Scripts\activate`. Install `python/requirements-torch.txt` for PyTorch tours or `python/requirements-jax.txt` for the JAX/Flax tour. The setup cells also check for these dependencies.

Default examples use bundled or simulated data. The deep-texture tour downloads pretrained VGG-19 ImageNet weights once (about 550 MB). Its default image size runs on a CPU; larger experiments benefit from a GPU. Diffusion examples use modest particle counts and networks to remain practical on a laptop.

Random seeds are set in the notebooks. Restart the kernel before running a tour from beginning to end. Floating-point and GPU differences can produce small variations across platforms.

## Collection boundaries

Only notebooks immediately inside `python/` form the maintained collection. The `todo/` directory contains 128 historical, unfinished conversions, including MATLAB-like code; these are **not runnable Python tours** and are not included in the catalogue or the passing execution count. Checkpoint copies are editor backups. The `nt_solutions/` directory remains only for compatibility with historical drafts; the maintained notebooks do not import or execute it.

See [the maintenance guide](../maintenance/README.md) for repeatable validation and website export.
