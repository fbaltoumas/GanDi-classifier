#!/usr/bin/env python3
"""Check whether this CPU supports the instruction baseline numpy's official
PyPI wheels require (SSE3/SSSE3/SSE4.1/SSE4.2/POPCNT, aka "x86-64-v2"), and if
not, rebuild numpy from source with a lower SIMD baseline so it can run here.

This CPU-feature gap is most commonly seen on virtual machines where the
hypervisor masks CPU features from the guest for live-migration compatibility
(e.g. Hyper-V's "Processor Compatibility Mode"), even though the physical CPU
fully supports them. If you can, fixing that VM/hypervisor setting is the
better long-term solution: it also restores the performance this workaround
gives up. See:
https://learn.microsoft.com/en-us/windows-server/virtualization/hyper-v/configure-processor-compatibility-mode

Usage:
    python3 scripts/check_cpu_baseline.py

Run this before `pip install .` / `pip install gandi` on a fresh setup, or
afterwards to repair an already-broken install (no need to reinstall the
`gandi` package itself afterwards).
"""
import platform
import subprocess
import sys


# Linux's /proc/cpuinfo reports SSE3 under its historical alias "pni"
# ("Prescott New Instructions"), never as a literal "sse3" flag.
REQUIRED_FLAGS = ["pni", "ssse3", "sse4_1", "sse4_2", "popcnt"]


def get_cpu_flags():
    if platform.system() != "Linux":
        return None
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("flags"):
                    return set(line.split(":", 1)[1].split())
    except OSError:
        return None
    return None


def rebuild_numpy():
    print("Rebuilding numpy from source with a lower SIMD baseline (this may take a minute)...")
    cmd = [
        sys.executable, "-m", "pip", "install", "--force-reinstall",
        "numpy", "--no-binary", "numpy", "-Csetup-args=-Dcpu-baseline=none",
    ]
    subprocess.run(cmd, check=True)
    print("Done. numpy has been rebuilt for this CPU.")


def main():
    flags = get_cpu_flags()
    if flags is None:
        print(
            "Could not read CPU flags (not on Linux, or /proc/cpuinfo unavailable) - skipping "
            "automatic detection. If you hit the numpy baseline RuntimeError, run this "
            "instead:\n"
            '  pip install --force-reinstall numpy --no-binary numpy -Csetup-args=-Dcpu-baseline="none"'
        )
        return

    missing = [flag for flag in REQUIRED_FLAGS if flag not in flags]
    if not missing:
        print(f"CPU supports the required baseline instructions ({', '.join(REQUIRED_FLAGS)}). Nothing to do.")
        return

    print(
        f"This CPU is missing: {', '.join(missing)} (required by numpy's official PyPI wheels).\n"
        "This is commonly caused by a hypervisor masking CPU features from this VM."
    )
    rebuild_numpy()


if __name__ == "__main__":
    main()
