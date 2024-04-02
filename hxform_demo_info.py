import json
from hxform import hxform as hx
from hxform import xprint

libs = hx.known_libs(info=False)
for lib in libs:
  xprint(f"{lib}")
  xprint(f"  coord systems: {hx.known_transforms(lib)}")

xprint("Library information:")
infos = hx.known_libs(info=True)
xprint(json.dumps(infos, indent=2))