"""Apply curated notebook introductions, reading lists, and portable setup cells."""

import json
from pathlib import Path
import re
from catalog_data import catalog, REFERENCES

ROOT = Path(__file__).resolve().parents[1]
SETUP = """# Locate the companion toolbox locally, or fetch it for a standalone/Colab copy.
from pathlib import Path
import importlib.util
import os
import subprocess
import sys

working = Path.cwd()
candidates = [working, working / "python", working.parent / "python"]
python_dir = next((p for p in candidates if (p / "nt_toolbox").is_dir()), None)
if python_dir is None:
    checkout = working / "numerical-tours-support"
    if not checkout.exists():
        subprocess.run(["git", "clone", "--depth", "1", "--branch", "master",
                        "https://github.com/gpeyre/numerical-tours.git", str(checkout)], check=True)
    python_dir = checkout / "python"
os.chdir(python_dir)
if str(python_dir) not in sys.path:
    sys.path.insert(0, str(python_dir))
{dependency_setup}
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(0)
plt.rcParams.update({{"figure.figsize": (8, 4), "figure.dpi": 100,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "font.size": 11, "image.cmap": "gray"}})
%matplotlib inline
"""


def cell(kind, text, id):
    value = dict(cell_type=kind, source=text, id=id, metadata={})
    if kind == "code":
        value.update(outputs=[], execution_count=None)
    return value


for slug, info in catalog().items():
    path = ROOT / "python" / (slug + ".ipynb")
    nb = json.loads(path.read_text())
    for original in nb["cells"]:
        original["source"] = "".join(original["source"])
    cells = [c for c in nb["cells"] if not c["id"].startswith("nt-")]
    # Existing citation keys are kept so that references in the mathematical text remain resolvable.
    existing = []
    for original in nb["cells"]:
        if original.get("id") == "nt-references":
            existing.extend(
                re.findall(
                    r"(?m)^[*\-] \[[^\]]+\].+(?:\n(?!\s*\n|[*\-] ).+)*",
                    "".join(original["source"]),
                )
            )
    cleaned = []
    for c in cells:
        source = c["source"]
        if c["cell_type"] == "markdown" and re.match(
            r"^#*\s*(Bibliography|References)\s*\n", source
        ):
            existing.extend(
                re.findall(r"(?m)^[*\-] .+(?:\n(?!\s*\n|[*\-] ).+)*", source)
            )
            continue
        if not cleaned and c["cell_type"] == "markdown":
            source = re.sub(r"^# .+\n?", "", source)
            source = re.sub(r"^[^\n]+\n[=\-]{3,}\n?", "", source)
        if c["cell_type"] == "markdown":
            source = source.replace("Exercise", "Worked example").replace(
                "exercice", "exemple commenté"
            )
            source = re.sub(r"(?i)(click|clic) [^\n]*(?:solution)[^\n]*", "", source)
            source = re.sub(
                r"(?m)^Write the code of (.+?)\.",
                r"The implementation below computes \1.",
                source,
            )
            source = re.sub(
                r"(?m)^Write the code (.+?)\.",
                r"The following cells implement the computation \1.",
                source,
            )
            source = source.replace("You need to", "We now").replace(
                "You should", "We can"
            )
        if c["cell_type"] == "code":
            source = (
                source.replace(
                    "import torch\n",
                    "import torch\ntorch.manual_seed(0)\ntorch.set_num_threads(1)\n",
                )
                if "torch.manual_seed" not in source
                else source
            )
        c["source"] = source
        if source.strip():
            cleaned.append(c)
    dependency_setup = """requirements = python_dir / "requirements.txt"
if any(importlib.util.find_spec(name) is None for name in ["numpy", "scipy", "matplotlib", "skimage", "sklearn", "pywt", "ipywidgets", "cvxpy", "skfmm", "autograd", "progressbar", "celer"]):
    subprocess.run([sys.executable, "-m", "pip", "install", "-r", str(requirements)], check=True)
"""
    code = "\n".join(c["source"] for c in cleaned if c["cell_type"] == "code")
    if "import torch" in code:
        dependency_setup += """if importlib.util.find_spec("torch") is None or importlib.util.find_spec("torchvision") is None:
    subprocess.run([sys.executable, "-m", "pip", "install", "-r", str(python_dir / "requirements-torch.txt")], check=True)
"""
    if "import jax" in code:
        dependency_setup += """if any(importlib.util.find_spec(name) is None for name in ["jax", "flax", "optax"]):
    subprocess.run([sys.executable, "-m", "pip", "install", "-r", str(python_dir / "requirements-jax.txt")], check=True)
"""
    setup_text = "## Run this tour\n\nRun the cells in order with a Python 3 kernel. The first cell locates the companion data and toolbox and installs missing dependencies when needed. All worked examples include their implementation directly in this notebook. Random seeds make comparisons reproducible; you can change them to explore other samples.\n"
    if slug == "ml_8_deep_texture_synthesis":
        setup_text += "\nThe VGG-19 experiment downloads pretrained ImageNet weights on its first run (about 550 MB). A CPU is sufficient for the default 96-pixel example; larger images benefit from a GPU.\n"
    elif slug.startswith("ml_12"):
        setup_text += "\nThe default experiment uses 2,000 particles and a small neural network so that the complete calculation remains practical on a CPU. Increase these sizes for a more accurate density estimate.\n"
    refs = existing[:8]
    for key in info["references"]:
        authors, title, publication, url, note = REFERENCES[key]
        if any(f"[{title}]".lower() in old.lower() or url in old for old in refs):
            continue
        refs.append(f"- {authors}. [{title}]({url}). {publication}. {note}")
        if len(refs) >= max(5, min(10, len(existing) + len(info["references"]))):
            break
    references = "## References and further reading\n\n" + "\n\n".join(refs) + "\n"
    title = "# " + info["title"] + "\n\n" + info["intro"] + "\n"
    if "canonical" in info:
        title += f"\nThis historical variant is retained for existing links. The main version is [{catalog()[info['canonical']]['title']}]({info['canonical']}.ipynb).\n"
    nb["cells"] = (
        [
            cell("markdown", title, "nt-introduction"),
            cell("markdown", setup_text, "nt-running"),
            cell("code", SETUP.format(dependency_setup=dependency_setup), "nt-setup"),
        ]
        + cleaned
        + [cell("markdown", references, "nt-references")]
    )
    nb["metadata"]["numerical_tours"] = {
        k: v for k, v in info.items() if k != "references"
    }
    path.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n")
print(f"Updated {len(catalog())} notebooks.")
