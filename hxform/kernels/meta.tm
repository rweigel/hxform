KPL/MK

Simple meta-kernel to load kernels necessary to demonstrate basic
dynamic frame transformations.

Start by defining basic variables to reference paths to kernels.
This is pretty basic, just point to the root directory containing
the kernels of interest.  This allows a single change in this
meta-kernel to "re-root" its contents--a common best-practice.

\begindata
   PATH_SYMBOLS += 'ROOT'
   PATH_VALUES  += './spiceypy/kernels'
\begintext

Now list the kernels to load in reverse priority order.  Last
always takes precedence over what is loaded first in the event
of a conflict.  This set of kernels is pretty basic, and there
are no conflicts.

Load the leapseconds kernel.  This contains various constants
and a tabulation of leapseconds that are necessary for
conversion between various time systems utilized within SPICE.

\begindata
   KERNELS_TO_LOAD += '$ROOT/naif0012.tls'
\begintext

Load the frame kernel of interest.  In this case we are using
a Van Allen Probes kernel that defines GSE and GSM.  The fit
to the IGRF magnetic dipole is stale and would need to be
updated.

\begindata
   KERNELS_TO_LOAD += '$ROOT/rbsp_general011.tf'
\begintext

The GSM/GSE frames need ephemeris data from the Earth and
the Sun to evaluate.  Provide it by loading a standard NAIF
distributed planetary SPK:

\begindata
   KERNELS_TO_LOAD += '$ROOT/de440s.bsp'
\begintext

The GEO frame requires the coefficients of the IAU orientation
model for the Earth.  The text PCK supplies these, so load that.

\begindata
   KERNELS_TO_LOAD += '$ROOT/pck00011.tpc'
\begintext

