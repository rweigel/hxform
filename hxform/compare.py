import sys

def compare(v, time, frame_in, frame_out, rep_in='car', rep_out='car', libs='all', libs_exclude=None, io=sys.stdout):

  def write(msg):
    results["log"] += msg + "\n"
    if io is not None:
      io.write(msg + "\n")

  import numpy as np
  import hxform

  libs_all = hxform.info.libs()
  if libs == 'all':
    libs = libs_all
  if libs_exclude is not None:
    libs = list(set(libs) - set(libs_exclude))

  libs.sort()

  results = {
    "time": time,
    "frame_in": frame_in,
    "frame_out": frame_out,
    "rep_in": rep_in,
    "rep_out": rep_out,
    "v_in": v,
    "v_out": {},
    "log": ""
  }

  write(f"Time: {time}")
  write(f"Transform: {frame_in} in {rep_in} => {frame_out} in {rep_out}")

  ml = np.max([len(lib) for lib in libs]) # Max length of lib name

  # TODO: Use PrettyTable or equivalent.
  spc = 9*" "
  pad = (ml+13)*" "
  templ = "\n{} x {} y {} z {} magnitude"
  write(templ.format(pad, spc, spc, 5*" "))
  write(83*"-")

  pad = (ml-5)*" "
  templ = "Input ({}): {} {:11.8f} {:11.8f} {:11.8f} {:11.8f}\n"
  write(templ.format(frame_in, pad, *v, np.linalg.norm(v)))
  write("Output ({}): {} ∠° wrt Input".format(frame_out, 55*" "))
  write(83*"-")
  vp = np.full((len(libs), 3), fill_value=np.nan)
  outputs = np.full((len(libs), 5), fill_value=np.nan)

  for i, lib in enumerate(libs):
    vp[i,:] = hxform.transform(np.array(v), time, frame_in, frame_out, lib=lib)
    mag = np.linalg.norm(vp[i,:])
    results["v_out"][lib] = vp[i,:]
    pad = (ml-len(lib)+5)*" "
    angle = (180/np.pi)*np.arccos(np.dot(v, vp[i,:])/(np.linalg.norm(v)*mag))
    outputs[i,0:3] = vp[i,:]
    outputs[i,3] = mag
    outputs[i,4] = angle
    templ = "{} {}   {:11.8f} {:11.8f} {:11.8f} {:11.8f} {:11.8f}"
    write(templ.format(lib, pad, *outputs[i,:]))

  max_diff = np.max(outputs, axis=0)
  min_diff = np.min(outputs, axis=0)
  write("\n")
  write("max-min:              {:11.8f} {:11.8f} {:11.8f} {:11.8f} {:11.8f}".format(*(max_diff-min_diff)))
  templ = "100*|max-min|/|max|:   {:10.4f}% {:10.4f}% {:10.4f}% {:10.4f}% {:10.4f}%"
  write(templ.format(*(100*np.abs(max_diff-min_diff)/np.abs(max_diff))))

  return results