from astropy.time import Time
# https://docs.astropy.org/en/stable/time/index.html

from hxform import xprint # print to console and astropy_time_demo.log

time_string = '2009-06-17T12:00Z'
time_object = Time(time_string, precision=9)

if False:
  # Basic usage
  t = t.tt   # https://docs.astropy.org/en/stable/time/index.html#id6
  val = t.jd # https://docs.astropy.org/en/stable/time/#time-format
  print(val)

# Print all formats and scales for t
xprint(f"\nTime: {time_object}\n")
maxlen = max([len(f) for f in list(Time.FORMATS.keys())])
xprint(f"Format        Scale Value")
xprint(f"------------------------------------------")
for FORMAT in [*list(Time.FORMATS.keys()),'jd1','jd2']:
  for SCALE in Time.SCALES:
    if SCALE == 'local': continue
    # Set scale. time_object = getattr(time_object, SCALE) is equivalent to doing
    #    time_object = time_object.SCALE
    # where SCALE is one of scales listed at
    # https://docs.astropy.org/en/stable/time/index.html#id6
    time_object = getattr(time_object, SCALE)

    # Set format. getattr(t, FORMAT) is equivalent to doing
    #    t = Time(time_string)
    #    val = t.format
    # where FORMAT is one of the formats listed at
    # https://docs.astropy.org/en/stable/time/#time-format
    val = getattr(time_object, FORMAT)

    xprint(f"{FORMAT:13} {SCALE:5} {val}") # https://docs.astropy.org/en/stable/time/#id6
