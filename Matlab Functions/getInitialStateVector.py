#==============================================================================#
# PRORGAM: getInitialStateVector.py
# AUTHOR: Justin Mansell (jmansell@purdue.edu), 2021
# DESCRIPTION: Turns a TLE into an inertial state vector
#==============================================================================#
from sgp4.api import Satrec

#---------------------------------INPUT TLE------------------------------------#
tleLine1 = '1 44420U 19036AC  21315.07330041  .00012333  00000-0  23653-2 0  9996'
tleLine2 = '2 44420  24.0061  99.8774 0017422   1.1246 358.9285 14.61908423125076'
#------------------------------------------------------------------------------#

# Create new satellite object
tle = Satrec.twoline2rv(tleLine1, tleLine2)

# Check initial epoch
jd = tle.jdsatepoch # Whole Julian day
jdf = tle.jdsatepochF # Fraction of Julian day

# Get cartesian state vector at TLE epoch
e, r, v = tle.sgp4(jd, jdf)

# Print
print("Epoch is: ", jd+jdf)
print("Position is:", r)
print("Velocity is:", v)

# Done
