"""Crop selected bands on their native grid. Run one manifest row per task."""
import json, os, subprocess, sys
from pathlib import Path

task = int(sys.argv[1])
manifest = sys.argv[2] if len(sys.argv)>2 else 'projection_climate_manifest.json'
entry = json.loads(Path(manifest).read_text())[task]
out = Path(entry['output'])
out.parent.mkdir(exist_ok=True)
done = out.with_suffix('.complete.json')
if done.exists() and out.exists():
    print('Already complete:', out)
    sys.exit(0)
partial = out.with_suffix('.partial.tif')
env = dict(os.environ, GDAL_HTTP_MAX_RETRY='3', GDAL_HTTP_RETRY_DELAY='10',
           GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR', GDAL_CACHEMAX='256',
           CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif')
cmd = ['gdal_translate','-of','GTiff','-co','COMPRESS=DEFLATE',
       '-co','TILED=YES','-co','BIGTIFF=IF_SAFER','-projwin']
cmd += [str(x) for x in entry['extent']]
for band in entry['bands']:
    cmd += ['-b', str(band)]
cmd += [entry['source'], str(partial)]
print('Source:', entry['source'], 'bands:',entry['bands'], flush=True)
subprocess.run(cmd, check=True, env=env)
info = json.loads(subprocess.check_output(['gdalinfo','-json',str(partial)],env=env))
assert len(info['bands']) == len(entry['bands'])
assert abs(info['geoTransform'][1]-1/120)<1e-10
assert abs(info['geoTransform'][5]+1/120)<1e-10
partial.replace(out)
done.write_text(json.dumps(dict(manifest=entry,geometry=info),indent=2)+'\n')
print('COMPLETE',out,flush=True)
