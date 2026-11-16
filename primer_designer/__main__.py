"""Run primer-designer-lite as a module: ``python -m primer_designer ...``."""

from __future__ import annotations

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
