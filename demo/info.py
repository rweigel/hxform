import json
import hxform

libs = hxform.libs()
for lib in libs:
  hxform.xprint(f"{lib}")
  hxform.xprint(f"  frames: {hxform.frames(lib)}")

hxform.xprint("Library information:")
infos = hxform.lib_info()
hxform.xprint(json.dumps(infos, indent=2))