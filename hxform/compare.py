def compare(v, time, frame_in, frame_out, rep_in='car', rep_out='car', libs='all', libs_exclude=None):

  import numpy as np
  import hxform

  libs_all = hxform.info.known_libs(info=False)
  if libs == 'all':
    libs = libs_all
  if libs_exclude is not None:
    libs = list(set(libs) - set(libs_exclude))

  hxform.xprint(f"Time: {time}")
  hxform.xprint(f"Transform: {frame_in} in {rep_in} => {frame_out} in {rep_out}")

  ml = np.max([len(lib) for lib in libs]) # Max length of lib name

  # TODO: Use PrettyTable or equivalent.
  spc = 9*" "
  pad = (ml+13)*" "
  templ = "\n{} x {} y {} z {} mag {} angle [deg]"
  hxform.xprint(templ.format(pad, spc, spc, spc, 2*" "))
  hxform.xprint(83*"-")

  pad = (ml+1)*" "
  hxform.xprint("Input: {} {:11.8f} {:11.8f} {:11.8f} {:11.8f}\n".format(pad, *v, np.linalg.norm(v)))
  hxform.xprint("Output:\n" + 83*"-")
  vp = np.full((len(libs), 3), fill_value=np.nan)
  outputs = np.full((len(libs), 5), fill_value=np.nan)

  out_dict = {}
  for i, lib in enumerate(libs):
    vp[i,:] = hxform.transform(np.array(v), time, frame_in, frame_out, lib=lib)
    mag = np.linalg.norm(vp[i,:])
    out_dict[lib] = vp[i,:]
    pad = (ml-len(lib)+5)*" "
    angle = (180/np.pi)*np.arccos(np.dot(v, vp[i,:])/(np.linalg.norm(v)*mag))
    outputs[i,0:3] = vp[i,:]
    outputs[i,3] = mag
    outputs[i,4] = angle
    templ = "{} {}   {:11.8f} {:11.8f} {:11.8f} {:11.8f} {:11.8f}"
    hxform.xprint(templ.format(lib, pad, *outputs[i,:]))

  max_diff = np.max(outputs, axis=0)
  min_diff = np.min(outputs, axis=0)
  hxform.xprint("\n")
  hxform.xprint("max-min:              {:11.8f} {:11.8f} {:11.8f} {:11.8f} {:11.8f}".format(*(max_diff-min_diff)))
  hxform.xprint("100*|max-min|/|max|:   {:10.4f}% {:10.4f}% {:10.4f}% {:10.4f}% {:10.4f}%".format(*(100*np.abs(max_diff-min_diff)/np.abs(max_diff))))

  return out_dict