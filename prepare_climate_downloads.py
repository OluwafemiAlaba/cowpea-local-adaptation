#!/usr/bin/env python3
"""Create portable WorldClim manifests; downloading is a separate explicit step."""
from pathlib import Path
import json
root=Path(__file__).resolve().parent/'manuscript_reanalysis/GF_validation_design_review'
for name in ['projection_climate_manifest.json','additional_climate_manifest.json','ssp245_climate_manifest.json']:
 entries=json.loads((root/name).read_text())
 for e in entries:
  url=e.get('url')
  if not url and e.get('source','').startswith('/vsicurl/'):url=e['source'][len('/vsicurl/'):]
  if not url or not url.startswith('https://geodata.ucdavis.edu/'):
   raise ValueError('Expected a documented WorldClim URL: '+str(e))
  e['source']='/vsicurl/'+url
 out=root/(Path(name).stem+'_portable.json')
 out.write_text(json.dumps(entries,indent=2)+'\n')
 print(out.name,len(entries),'tasks')
print('From the GF directory, run python3 crop_projection_climate.py TASK MANIFEST_portable.json for each zero-based task. Requires GDAL and network access. This script does not download data.')
