################################################################################
#                                                                              #
# This file is distributed as part of the FIEStool data reduction package.     #
#                                                                              #
# Copyright (c) 2005 NOTSA; All Rights Reserved;                               #
# Original author: Eric Stempels; See COPYRIGHT file for details               #
#                                                                              #
################################################################################

################################################################################
#
# FIEStool
#
# checkImport.py
#
# Last update 18/12/2007
#
################################################################################

import sys

deferror = "FIEStool failed to start!"

def numericVersion(version):
    """
    Returns the first part of a version number that comprises only numbers and
    (optionally) dots.
    """
    for candidate in version.strip().split():
        try:
            lst = [int(x) for x in candidate.split('.')]
            break
        except ValueError:
            pass
    else:
        raise ValueError('Not a valid version "%s"' % version)

    return lst

def checkImport(package, minversion):

  """
  Checks that the minimum version of required software packages is installed
  """

  try:
    thispackage = __import__(package)
  except ImportError:
    print deferror
    print "Required external package '%s' not installed" % package
    raise

  complies = True

  try:
    required = numericVersion(minversion)
    current =  numericVersion(thispackage.__version__)
    maxl = max(len(required), len(current))
    dr, dc = maxl - len(required), maxl - len(current)
    if dr > 0: required = required + ([0] * dr)
    if dc > 0: current = current + ([0] * dc)

    for x, y in zip(required, current):
        if y > x:
            break
        if x > y:
            complies = False
            break
  except ValueError:

    print "ValErr"
    # Fallback: Compare strings, as the integer didn't work.
    #           Not reliable, but is better than nothing
    complies = not (thispackage.__version__ < minversion)

  if not complies:
    print deferror
    print "Version of external package '%s' lower than minimum required version of %s" % (package, minversion)
    print "Your version of '%s' is %s" % (package, thispackage.__version__)
    raise ImportError

  pass





