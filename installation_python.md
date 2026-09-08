---
layout: page
title: Run a Python tour
description: Start in your browser, or use a local Jupyter environment.
---

## Open in Colab

Choose a tour from the [Python collection]({{ '/python/' | relative_url }}) and select **Open in Colab**. Run the first setup cell, then run the remaining cells in order. The setup locates the companion toolbox and data, downloads the repository if needed, and installs missing dependencies.

Every worked example includes its complete implementation. You can change parameters, rerun the cells, and compare your results with the figures in the reading version.

## Run locally

Use Python 3.12 and a separate environment. Download the [complete repository](https://github.com/gpeyre/numerical-tours/archive/refs/heads/master.zip), or clone it:

```sh
git clone https://github.com/gpeyre/numerical-tours.git
cd numerical-tours
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r python/requirements.txt jupyterlab
jupyter lab python
```

On Windows, activate the environment with `.venv\Scripts\activate` instead. In Jupyter, choose the environment’s Python kernel and open a notebook from the `python` folder.

For the PyTorch tours, install `python/requirements-torch.txt`. For the JAX tour, install `python/requirements-jax.txt`. The setup cell also checks for these dependencies. The deep-texture tour downloads pretrained VGG-19 weights once; the other default examples use bundled data or simulated samples.

## Reproduce a result

Restart the kernel and run all cells from the beginning. The notebooks set random seeds and use CPU-friendly defaults. Changing a seed or an algorithm parameter may change the output; GPU implementations can also introduce small numerical differences.

The [repository validation guide](https://github.com/gpeyre/numerical-tours/tree/master/maintenance) explains how to execute the entire collection and inspect failures. A successful execution checks that all cells complete; numerical invariants additionally check key shared operations.
