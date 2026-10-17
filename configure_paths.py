#!/usr/bin/env python3
"""Check the portable project layout; does not rewrite code or data."""
from pathlib import Path
import os
root = Path(os.environ.get('COWPEA_ROOT', '~/data/cowpealandrace_460')).expanduser()
print('Project root:', root)
analysis = root / 'manuscript_reanalysis'
if not analysis.is_dir():
    raise SystemExit('Extract both archives under ~/data/, or set COWPEA_ROOT to the extracted project directory.')
print('Analysis directory found. Run statistical scripts from their documented working directories.')
