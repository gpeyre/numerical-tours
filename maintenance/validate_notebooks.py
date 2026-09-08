"""Execute each published Python tour in a fresh kernel; never skip failed cells."""

from __future__ import annotations
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def execute_one(path: Path, output: Path, timeout: int) -> dict:
    import nbformat
    from nbclient import NotebookClient

    start = time.monotonic()
    nb = nbformat.read(path, as_version=4)
    result = {
        "notebook": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    try:
        nbformat.validate(nb)
        NotebookClient(
            nb,
            timeout=timeout,
            kernel_name="python3",
            resources={"metadata": {"path": str(path.parent)}},
            allow_errors=False,
            record_timing=True,
        ).execute()
        result["status"] = "passed"
    except Exception as exc:
        result.update(status="failed", error=str(exc)[-6500:])
    result["seconds"] = round(time.monotonic() - start, 2)
    result["code_cells"] = sum(c.cell_type == "code" for c in nb.cells)
    result["executed_cells"] = sum(
        c.cell_type == "code" and c.execution_count is not None for c in nb.cells
    )
    output.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, output / path.name)
    (output / (path.stem + ".json")).write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "notebooks", nargs="*", help="Names or paths; default: all main Python tours"
    )
    parser.add_argument("--output", type=Path, default=ROOT / "reports" / "executed")
    parser.add_argument(
        "--timeout", type=int, default=180, help="Timeout per cell in seconds"
    )
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    paths = (
        [
            ROOT / "python" / (n if n.endswith(".ipynb") else n + ".ipynb")
            for n in args.notebooks
        ]
        if args.notebooks
        else sorted((ROOT / "python").glob("*.ipynb"))
    )
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/numerical-tours-matplotlib")
    os.environ.setdefault("IPYTHONDIR", "/private/tmp/numerical-tours-ipython")
    os.environ.setdefault("JUPYTER_RUNTIME_DIR", "/private/tmp/numerical-tours-jupyter")
    for var in [
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ]:
        os.environ[var] = "1"
    os.environ["PATH"] = (
        str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]
    )
    args.output = args.output.resolve()
    if args.worker:
        result = execute_one(paths[0], args.output, args.timeout)
        print(result["notebook"], result["status"], result["seconds"], flush=True)
        return

    def worker(path):
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            path.name,
            "--worker",
            "--output",
            str(args.output),
            "--timeout",
            str(args.timeout),
        ]
        try:
            process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=max(600, args.timeout * 5),
            )
            report_path = args.output / (path.stem + ".json")
            if process.returncode or not report_path.exists():
                raise RuntimeError(process.stderr[-2000:])
            result = json.loads(report_path.read_text())
        except Exception as exc:
            result = {
                "notebook": str(path.relative_to(ROOT)),
                "status": "failed",
                "error": str(exc),
            }
        print(path.stem, result["status"], result.get("seconds", ""), flush=True)
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(worker, paths))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output.parent / "execution.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(f"{sum(r['status'] == 'passed' for r in results)}/{len(results)} passed")
    sys.exit(any(r["status"] != "passed" for r in results))


if __name__ == "__main__":
    main()
