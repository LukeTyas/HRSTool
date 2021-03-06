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
# fileUtils.py
#
# Last update 26/09/2005
#
################################################################################

"""
   File-related utility routine(s) that make life easier.
"""

# Import external modules

import os
import shutil
import fnmatch
import messageLog


################################################################################


def removefiles(*patterns):

  """
     Remove files with UNIX-like filename pattern matching.
     No questions asked, so be careful!
  """

  # 'patterns' is an array or tuple of filename patterns. Can, of course,
  # also be inidividual file names.

  # Loop over patterns
  for pattern in patterns:

    # Split off the directory path
    (dirname, filepattern) = os.path.split(pattern)
    
    # Get a list of all files in the directory
    try:
      filelist = os.listdir(dirname)
    except OSError, errstr:
      # Succeed even if the directory was not there (and put warning in log)
      messageLog.put_warning(errstr, 9)
      return

    # Check each file in the directory list
    for file in filelist:

      # And see if it's name matches the pattern
      if fnmatch.fnmatch(file, filepattern):

        # If yes, (try to) remove it from the system
	messageLog.put('Removing file : "%s"' % file, 7)
        try:
          os.remove(os.path.join(dirname, file))
	except OSError, errstr:
          # Succeed even if there were no files (and put warning in log)
	  messageLog.put_warning(errstr, 9)


def movefile(oldname, newname):

  """
     Move file to new name, overwriting and no error reporting
  """
  
  try:
    shutil.move(oldname, newname)
  except:
    pass

