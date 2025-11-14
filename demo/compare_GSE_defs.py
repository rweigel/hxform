import numpy
import hxform

"""Two definitions of GSE compared.

[1] https://docs.google.com/presentation/d/1DpX0xfcwDGiyf148AF7Lo0cEHVXVJ3v5
[2] Meeus 1991 https://dn710207.ca.archive.org/0/items/astronomicalalgorithmsjeanmeeus1991/Astronomical%20Algorithms-%20Jean%20Meeus%20%281991%29.pdf

GSE_A:
  u unit vector from center of Earth to center of Sun.
  v unit vector perpendicular to ecliptic pole, positive northward.
  The angle between u and v is not 90 degrees. It can differ by up to
  1.2 arcseconds (3.3 x 10^-4 degrees) (Meeus 1991, [2], p152).

  X = u ("X primary")

  Y = v ⨯ X/|v ⨯ X| (Divide by norm because v and u are not exactly
                     perpendicular.)

  Z = X ⨯ Y (Z is not exactly aligned with ecliptic pole, v; X and Y are 
             perpendicular and unit length by construction, so division by
             norm is not needed.)

GSE_B:
  Z = v ("Z primary"; aligned with ecliptic pole, v)

  Y = Z ⨯ u/|Z ⨯ u| (So Y is perpendicular to Z and u.)

  X = Y ⨯ Z (X is not exactly aligned with unit vector from center of Earth
             to center of Sun, u)

In the following, we transform the unit vector along the Z axis in GSE_A to
its transform in GSE_B, and compute the angle between the two vectors. We find
0.69 arcseconds, consistent with the claim that the angle between u and v
differs from 90 degrees by up to 1.2 arcseconds and near the value of 0.72
in the example on page 153 of Meeus 1991 [2].

To compute the angle between u and v, we transform Z in GSE_A to GSE_B. The angle
between Z in A and its transform in B is the same as the angle between u and v.
"""

lib     = 'spiceypy2' # Uses kernel 
time    = '1992-10-13T12:00:00Z' # (Meeus 1991, [2], p153).
input   = [0, 0, 1]
initial = 'GSE'             # GSE_A
final   = 'GSE_Z_PRIMARY'   # GSE_B

time = hxform.timelib.iso2ints(time)

output = hxform.transform(input, time, initial, final, lib=lib)

hxform.xprint(f"Time:   {time}")
hxform.xprint(f"Input:  {input}")
hxform.xprint("")
hxform.xprint("Transform: Z in GSE_A => GSE_B")
hxform.xprint(f"  Output:     {output}")
hxform.xprint(f"  1-|Output|: {1 - numpy.linalg.norm(output):.2e}")

mag1 = numpy.linalg.norm(input)
mag2 = numpy.linalg.norm(output)
angle = (180/numpy.pi)*numpy.arccos(numpy.dot(input, output)/(mag1*mag2))
hxform.xprint("")
hxform.xprint("Angle between Z in GSE_A and its transform to GSE_B:")
hxform.xprint(f"  {(angle):.4e} [degrees]")
hxform.xprint(f"  {(angle)*3600:.4e} [arcseconds]")


"""
Time:   [1992, 10, 13, 12, 0, 0]
Input:  [0, 0, 1]

Transform: Z in GSE_A => GSE_B
  Output:     [-3.341712807303221e-06, -5.551115123125783e-17, 0.9999999999944162]
  1-|Output|: 2.22e-16

Angle between Z in GSE_A and its transform to GSE_B:
  1.9147e-04 [degrees]
  6.8928e-01 [arcseconds]
"""