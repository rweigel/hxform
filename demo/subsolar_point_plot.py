import utilrsw

def plot(times, thetas, frame='GEO', libs=None, lib_ref=None):

  import matplotlib.pyplot as plt
  from matplotlib.dates import DateFormatter

  def label_escaped(lib, lib_ref=None):
    lib_escaped = lib.replace('_', r'\_')
    lib_escaped = f'\\text{{{lib_escaped}}}'
    label = fr'$λ_{{{lib_escaped}}}$'
    if lib_ref is not None:
      lib_ref_escaped = f'\\text{{{lib_ref}}}'
      label = fr'{label} - $λ_{{{lib_ref_escaped}}}$'
    return label

  for lib in libs:
    if lib == lib_ref:
      continue

    gs = plt.gcf().add_gridspec(2, hspace=0.07)
    axes = gs.subplots(sharex=True)

    axes[0].plot(times, thetas[lib], label=label_escaped(lib))
    axes[0].set_ylabel('(deg)')
    axes[0].set_title(f'Subsolar colatitude in {frame} frame')
    axes[0].grid(True)
    axes[0].legend()
    axes[0].get_yaxis().get_major_formatter().set_useOffset(False)

    axes[1].plot(times, thetas[lib]-thetas[lib_ref], label=label_escaped(lib, lib_ref))
    axes[1].grid(True)
    axes[1].set_ylabel('(deg)')
    axes[1].legend()
    axes[1].get_yaxis().get_major_formatter().set_useOffset(False)
    axes[1].xaxis.set_major_formatter(DateFormatter('%Y'))
    axes[1].set_xlabel('Year')

    utilrsw.mpl.savefig(f'subsolar_point/{frame}_{lib}', formats=['png'])


with open('subsolar_point/subsolar_point.pkl', 'rb') as f:
  import pickle
  print("Reading subsolar_point/subsolar_point.pkl")
  data = pickle.load(f)

times = data['thetas'].index.to_pydatetime()
opts = data['opts']

plot(times, data['thetas'],
     frame=opts['frame'], libs=opts['libs'], lib_ref=opts['lib_ref'])
