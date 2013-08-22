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
# tasks.py
#
# Last update 08/01/2009
#
################################################################################

"""
   Main module containing all the steps that are involved to properly reduce
   an observed frame of FIES. The routines are valid for the general case of
   cross-dispersed echelle spectrographs, and should work for other instruments.
   
   These routines are written in modular form (single pass), to increase
   legibility and make it easier to modify the software in the future.
   
   Additional tasks can be defined by deriving a new class from the 'task'
   superclass, and filling in the procedure.
   
   Docstrings (similar to this one) in the tasks are transferred to the on-line
   help popup windows, so these should be good and accurate.
"""

# Import necessary modules

# System modules
import os
import config

# Interact with FITS files
import pyfits

# Math module for efficient array processing
import numpy
#import numarray.nd_image

# Other modules with which communication is needed (the logging facility
# and the helper to autoload configuration files)
import messageLog
import autoLoader

# Utilities
import frameUtils
import waveUtils
import plotUtils
import fileUtils

# And, most important, the module that provides PyRAF-based wrapping routines
# to interact with IRAF.
import irafWrappers

# Extract the version number from FIEStool
from FIEStool import _version as FIEStool_version


################################################################################


class taskAbort(Exception): 
  """Module-specific error, to be raised when an the reduction is to be aborted
     in a friendly way. This error should be caught by the task manager, who
     then stops the processing.
  """
  def __init__(self, *args):
    self.args = args


class taskError(Exception): 
  """Module-specific error, to be raised when an 'expected' exception occurs
     during the processing of a task. This error should then be caught by the
     task manager, avoiding nasty error output.
  """
  def __init__(self, *args):
    self.args = args


################################################################################
# Top-level class 'task' provides skeleton for tasks and basic functionality
################################################################################

class task:
  """Description not given"""

  # The name with which the task will be referred to, normally identical to
  # the class name.
  name 	       = 'undefined'

  # The text that should appear on the taskbutton
  buttonText   = 'Not defined'

  # This tag is used by the pipeline processing to append a suffix
  # to the name of the file that is being reduced.
  suffix       = None

  # By default, tasks can be executed inside threads. Currently some tasks
  # MUST be executed by the main thread, because they fail otherwise. This is
  # a known (reported) bug in PyRAF, and is supposed to be fixed in the
  # autumn of 2005.
  inthread     = 1

  # Does this task produce an output file, or is it just a sidestep in the
  # pipeline processing (for example, a plotting routine does not produce a
  # new FITS file, only the on-screen plot)
  output       = 1

  # A list of tasks that have to be completed in the pipeline before the current
  # task can be invoked.
  prereq       = []


  # Hook. Here the routines for executing the tasks will be attached.
  def __init__(self):
    pass



################################################################################
# Dummy/Testing/Example tasks for development purposes
################################################################################

class testtask(task):

  """
     This task does nothing (except for writing to the log file). The help text
     that should appear in the interface should go here
  """

  # See the task superclass for info on which properties may be defined. The
  # following to are mandatory
  name         = 'skeleton'
  buttonText   = 'Does really nothing!'

 
  # Task that run in the pipeline (i.e. not as a 'calibration task') must
  # have an inframe and an outframe parameter. Calibration tasks (see further
  # below) only have an inframe defined. Their output files are defined through
  # the list of configuration options.
 
  def run(self, inframe, outframe):

    # Clear, isn't it?
    messageLog.put('This task does nothing (except for this echo)...')


class testerror(task):

  "This task only causes an error..."

  name         = 'testtask2'
  buttonText   = 'Cause error!'

  def run(self, inframe, outframe):
    
    # Slightly more difficult, but still quite simple
    messageLog.put('This task only causes an error...')
    raise taskError, 'Testing error handling!'

    

################################################################################
# plotspec is the main plotting routine, using either Biggles or IRAF
################################################################################

class plotspec(task):
  """Plot (part of) an extracted spectrum.
  """

  name         = 'plotspec'
  buttonText   = 'Plot reduced spectrum'
  inthread     = 0

  def run(self, inframe, outframe=None, usebiggles=1, *dummy):

    data = pyfits.open(inframe)
    npixels = data[0].header['NAXIS1']
    
    wvlrange = [ config.options['plot_startwave'].value,
                 config.options['plot_endwave'].value ]

    order = config.options['plot_defaultorder'].value

    # Is this a two-dimensional spectrum?
    if data[0].header.has_key('NAXIS2'):

      # Yes, we're dealing with a two-dimensional spectrum
      norders = data[0].header['NAXIS2']
      
      # Read the wavelength definition
      wavedef = waveUtils.getwavedef2D(data)

      # Did we successfully read the wavelength definition?
      if (wavedef is not None):

	# Disregards from the order numbering. Here we're only interested in the wavelengths
        wavedef = wavedef[:, 2:]

	# Yes - so make an array of the starting and ending wavelengths
	# of each order
	wavedefrange = wavedef.copy()
	wavedefrange[:, 1] = wavedefrange[:, 0] + npixels * wavedef[:, 1]

	# If either wavelength range is zero, make it min or max
	if (wvlrange[0] == 0):
	  wvlrange[0] = wavedefrange.min()
	if (wvlrange[1] == 0):
	  wvlrange[1] = wavedefrange.max()

	if (wvlrange[1] - wvlrange[0] <= 0):
	  raise taskError, 'Cannot plot: Invalid wavelength range (%.1f - %.1f)' % tuple(wvlrange)

	# Determine the lowest and highest order number to plot for this range.
	# If the input cannot be read (for example due to wrong format), take
	# lowest and/or highest available order number
	try:
	  minorder = (wavedefrange[:, 1] > wvlrange[0]).nonzero()[0].min()
	except ValueError:
          minorder = 0

	try:
	  maxorder = (wavedefrange[:, 0] < wvlrange[1]).nonzero()[0].max()
	except ValueError:
          maxorder = norders - 1

	# Now create the vector of wavelengths for these orders
	wave = numpy.zeros((maxorder-minorder+1, npixels), dtype=numpy.float64)
	for i in range(minorder, maxorder+1):
	  wave[i-minorder, :] = wavedef[i, 0] + wavedef[i, 1] * numpy.arange(npixels)


      else:
        # No, we could not read the wavelength definition
	# Use order number instead

        wave = None
        wvlrange = None

        # Does the requested order number exist at all? If not, use default (1)
	if order-1 not in range(norders):
	  messageLog.put_warning('Order number out of range. Using default : 1', 5)
	  order = 1

        # This time, plot only one single order
	minorder = maxorder = order - 1


      # Extract the orders of the spectrum that will be plotted from the dataset
      if data[0].header.has_key('NAXIS3'):
        spec = data[0].data[0, minorder:maxorder+1, :]
      else:
        spec = data[0].data[minorder:maxorder+1, :]


    else:
      # Apparently this is not a two-dimensional, but a one-dimensional (merged) spectrum
      # This case is easier!
    
      norders = 1
      
      # The value below is required when calling IRAFs plotting task (splot) on 1D spectra
      minorder = 0

      # Get the 1D wavelength definition
      wavedef = waveUtils.getwavedef1D(data)

      # Did we successfully read the wavelength definition?
      if wavedef is not None:

	# Disregards from the order numbering. Here we're only interested in the wavelengths
	wavedef = wavedef[:, 2:]

        # Construct the wavelength vector
	allwave = numpy.zeros((1, npixels), dtype=numpy.float64)
	allwave[0, :] = wavedef[0, 0] + wavedef[0, 1] * numpy.arange(npixels)

	# If either wavelength range is zero, make it min or max
	if (wvlrange[0] == 0):
	  wvlrange[0] = allwave.min()
	if (wvlrange[1] == 0):
	  wvlrange[1] = allwave.max()

	if (wvlrange[1] - wvlrange[0] <= 0):
	  raise taskError, 'Cannot plot: Invalid wavelength range (%.1f - %.1f)' % tuple(wvlrange)

	# Determine the lowest and highest pixel for this wavelength range
	try:
	  minindex = (allwave[0, :] < wvlrange[0]).nonzero()[0].max()
	except ValueError, IndexError:
          minindex = 0

	try:
	  maxindex = (allwave[0, :] > wvlrange[1]).nonzero()[0].min()
	except ValueError, IndexError:
          maxindex = npixels

	# And create the vector of wavelengths
        wave = allwave[:, minindex:maxindex]

      else:
	# No, we could not read the wavelength definition
	# Plot entire spectrum instead

	wave = None
	wvlrange = None

	minindex = 0
	maxindex = npixels

      # Create a vector containing the spectrum
      allspec = numpy.zeros((1, npixels), dtype=numpy.float64)
      allspec[0, :] = data[0].data

      # Extract the orders of the spectrum that will be plotted from the dataset
      spec = allspec[:, minindex:maxindex]
      

    # Now (finally) call the plotting routine

    messageLog.put('Plotting spectrum', 5)

    if usebiggles:
      # Biggles :

      title = "%s - %s" % (os.path.basename(inframe), data[0].header['OBJECT'])

#     plotUtils.plot(spec, wave, title=title, xrange=wvlrange, psfile=inframe+'.png')
      plotUtils.plot(spec, wave, title=title, xrange=wvlrange)
    else:
      # Not Biggles :
      irafWrappers.specplot(inframe, order=minorder+1)

    data.close()


################################################################################
# The coming group of tasks are routines that will be used in the pipeline
################################################################################

class autoload(task):
  """Automatically load a pipeline configuration depending on
     the values of FITS headers. This is useful when processing a set
     of frames with different settings for binning or windowing
  """

  name         = 'autoload'
  buttonText   = 'AutoLoad configuration'

  def run(self, inframe, outframe):

    # Restore the active ruleset from the configuration    
    currentrules = config.options['autoloadrules'].value

    try:
      messageLog.put('Checking rules for autoloading')
      # Put the active ruleset into the autoLoader
      ruleset = autoLoader.ruleset(currentrules)
      # And execute the ruleset. If succesful, it will return the name of a
      # config file that should be loaded.
      newconfig = ruleset.run(inframe)
      if newconfig == 'abort':
	messageLog.put('Autoloading rules stops further processing of this frame')
        raise taskAbort
      elif newconfig:
        # Load data from config, but do not overwrite in- and outpath parameters
	messageLog.put('Autoloading configuration from %s' % newconfig)
	config.load(newconfig, ignore_options=['inpath', 'outpath',
			'filename_filter', 'plot_startwave', 'plot_endwave',
			'plot_defaultorder', 'currentmode'])
      else:
        pass
    except autoLoader.RuleError, errString:
      # AutoLoading didn't work out
      messageLog.put_warning('Could not perform autoloading :')
      messageLog.put_warning('Reason : %s' % errString)
      messageLog.put_warning('Trying to continue')

################################################################################

class preproc(task):
  """Do the initial processing of the obtained frames. Includes FITS file
     verification, clipping of pre- and overscan, padding of windowed frames
     to full size and adjustment with bad pixel mask
  """

  name	 	= 'preproc'
  buttonText	= 'Preprocess frame'
  suffix	= 'prep'

  def run(self, inframe, outframe):


    messageLog.put('Preprocessing object frame')


    # Convert from MEF to non-MEF if necessary
    messageLog.put('Extracting data from Multi-Extension FITS (block %i)' % config.options['mef_dataext'].value, 5)
    indata  = pyfits.open(inframe)
    outdata = frameUtils.extractMEF(indata, extension=config.options['mef_dataext'].value)

    if (outdata is None):
      raise taskError('FITS Extension %i not found in this frame' % config.options['mef_dataext'].value)

    # Test the file for FITS integrity, and fix headers as far as possible
    # This test is essential, because FIES files are currently (May 2005)
    # NOT up to FITS standard
    messageLog.put('Testing file for FITS integrity', 5)
    try:
      outdata.verify('exception')
    except: 
      messageLog.put_warning('FITS verification failed! Trying to fix...')
      try:
        outdata.verify('fix')
      except:
        # Could not fix... giving up!
	outdata.close()
	raise taskError('Could not fix FITS headers! Is this really a FITS file?')
      messageLog.put_warning('Fixing headers succeeded!')


    # Add history header card giving FIEStool version
    outdata[0].header.add_history('This is FIEStool vesion %s' % FIEStool_version)

    # Add history header card stating name of original frame
    outdata[0].header.add_history('Original file : %s' % os.path.basename(inframe))


    # Clip the overscan and other non-data regions from the image. The clipping
    # range is specified in the config settings.
    if (config.options['frameorientation'].value != 0) :
      messageLog.put('Rotating/flipping frame', 5)
      try:
        frameUtils.flipframe(outdata, config.options['frameorientation'].value)
      except:
        raise taskError('Could not flip frame - something wrong with input data')
      outdata[0].header.add_history('Frame rotated/flipped (orientation %i)' % config.options['frameorientation'].value)

    messageLog.put('Clipping frame', 5)
    frameUtils.clipframe(outdata)
    outdata[0].header.add_history('Frame clipped')


    # Check for possible windowing, and pad (extend) the frame if necessary...
    padded = frameUtils.padframe(outdata)

    # If the frame was padded to full size, set a header flag that this was
    # done. Padded frames will not have scattered light subtraction, because
    # it will give very bad fits.
    if padded:
      messageLog.put('Padding frame', 5)
      outdata[0].header.add_history('Frame padded with zeros (only if windowed)')
      outdata[0].header.update('PADDED', 1)
      messageLog.put_warning('Windowed frames may give unreliable background subtraction')

    # Use a pixel mask file to remove bad pixels from the image. Which way
    # bad pixels are treated is still pending (currently: give zero value)
    try:
      mask = pyfits.open(config.options['pixelmask'].value)
      mask = frameUtils.extractMEF(mask, extension=config.options['mef_dataext'].value)
      messageLog.put('Rotating/flipping bad pixel mask', 5)
      try:
        frameUtils.flipframe(mask, config.options['frameorientation'].value)
      except:
        raise taskError('Could not flip bad pixel mask - something wrong with input data')
      messageLog.put('Clipping bad pixel mask', 5)
      frameUtils.clipframe(mask)
      messageLog.put('Multiplying with bad pixel mask (putting bad pixels to zero)')
      try:
        outdata[0].data = outdata[0].data * mask[0].data
      except ValueError:
        messageLog.put_warning('Data size (%i, %i)' % outdata[0].data.shape)
        messageLog.put_warning('Bad pixel mask size (%i, %i)' % mask[0].data.shape)
        raise taskError('Incompatible sizes of object and bad pixel mask. Did you start binning?')
      outdata[0].header.add_history('Multiplied with bad pixel mask : %s' % os.path.basename(config.options['pixelmask'].value))
      mask.close()
    except IOError: 
      messageLog.put_warning('Could not open bad pixel mask. Not critical, continuing.')


    # Calculate the heliocentric velocity correction from the FITS headers
    # and store the result in a header card (doesn't hurt to do this, does it?) 
    messageLog.put('Calculating heliocentric velocity correction', 5)
    waveUtils.helcorr(outdata)


    # Remove outframes if it exists
    fileUtils.removefiles(outframe)


    # Write output to outframe (outdata object actually still points to input data)
    try:
      outdata.writeto(outframe)
    except IOError:
      raise taskError('Cannot write output to %s' % outframe)

    # Clean up
    outdata.close()


################################################################################

class updatethar(task):

  """Update the current wavelength solution if the frame being processed
     is a ThAr frame.
  """

  # See the task superclass for info on which properties may be defined. The
  # following to are mandatory
  name         = 'updatethar'
  buttonText   = 'Check for new ThAr solution'
  prereq       = ['preproc']
 
  def run(self, inframe, outframe):

    # Open the input frame
    data  = pyfits.open(inframe)

    # If this is a ThAr frame, try to refind the wavelength solution    
    if data[0].header['FILMP7'] == 1:

      # No need to keep this in memory
      data.close()

      # Do some sanity checks
      if not os.path.exists(config.options['masterwaveref'].value) :
        raise taskError, 'No master wavelength reference defined - cannot redetermine solution'
      if (config.options['masterwaveref'].value == config.options['waveref'].value) :
        messageLog.put_warning('Action would destroy master wavelength reference data!')
        raise taskError, 'Filenames of "master wavelength reference" and "wavelength reference" must be different'

      # OK, create a temporary frame to hold intermediate output
      tempframe = os.path.join(os.path.dirname(inframe), 'tempframe.fits')
      fileUtils.removefiles(tempframe)

      messageLog.put('Extracting orders from %s' % inframe)
      irafWrappers.extractwave(inframe,
                               tempframe,
			       config.options['masterorderdef'].value)

      messageLog.put('Determining the wavelength solution from %s' % inframe)
      (pixshift, wvlshift, ratio) = irafWrappers.findwaveshift(tempframe, config.options['masterwaveref'].value)

      # Do not do any correction if less than 10% of the lines were identified
      if ratio < 0.1:
	messageLog.put_warning('Could only reidentify %s %% of the lines - ignoring shift' % (ratio*100))
      else:
        # Since ratio is OK, not do it for real!
	messageLog.put('Reidentified %s %% of the lines' % (ratio*100))
	messageLog.put('Saving new wavelength solution')
	fileUtils.removefiles(config.options['waveref'].value)
	os.rename(tempframe, config.options['waveref'].value)
        irafWrappers.refindwavesol(config.options['waveref'].value, config.options['masterwaveref'].value)

      # Finish further processing (since this is a ThAr frame)
      messageLog.put('This ThAr frame needs no further processing... stopping here')
      raise taskAbort

    else:

      messageLog.put('This is not a ThAr frame... skipping task')



################################################################################

class subtractbias(task):
  """Subtract combined BIAS from image.
  """

  name         = 'subtractbias'
  buttonText   = 'Subtract BIAS'
  suffix       = 'bias'
  prereq       = ['preproc']

  def run(self, inframe, outframe):

    # Open the input frame
    data  = pyfits.open(inframe)

    # Read combined BIAS from disk
    try:
      bias = pyfits.open(config.options['masterbias'].value)
    except IOError:
      raise taskError('Cannot open combined BIAS : "%s"' % config.options['masterbias'].value)

    if (data[0].data.shape != bias[0].data.shape):
      messageLog.put_warning('Data size (%i, %i)' % data[0].data.shape)
      messageLog.put_warning('Bias size (%i, %i)' % bias[0].data.shape)
      raise taskError('Incompatible sizes of object and combined bias frame. Did you start binning?')

    # Replace all pixels values that are currently 0 (i.e bad pixels and padded
    # pixels) with the master BIAS value of that pixel.
    messageLog.put('Replace zero-valued pixels with combined bias values', 5)
    data[0].data = numpy.choose(data[0].data == 0, (data[0].data, bias[0].data))
    data[0].header.add_history('Replaced zero-valued pixels with combined bias values')

    # Subtract the combined BIAS
    messageLog.put('Subtracting combined bias')
    data[0].data = data[0].data - bias[0].data
    data[0].header.add_history('Subtracted combined bias : %s' % os.path.basename(config.options['masterbias'].value))

    # Write output to outframe (data object actually still points to input data)
    try:
      data.writeto(outframe)
    except IOError:
      raise taskError('Cannot write output to %s' % outframe)

    # Clean up
    bias.close()
    data.close()

################################################################################

class scattering(task):
  """Model and subtract the scattered light in an image. The level of scattered
     light is estimated from the level between the individual orders.
  """

  name		= 'scattering'
  buttonText	= 'Subtract scattering'
  suffix	= 'scat'
  prereq	= ['preproc']
  
  newmethod = True

  def run(self, inframe, outframe):

    # Remove outframe if it exists
    fileUtils.removefiles(outframe)

    # Has this frame been padded?
    indata = pyfits.open(inframe)
    padded = indata[0].header.get('PADDED', None)
    #indata.close()

    if padded:
      # Skip if it has been padded
      messageLog.put_warning('Frame is windowed! Skipping scattered light subtraction', 3)
      return

    # Use PyRAF wrapper routine to fit and subtract scattered light
    messageLog.put('Modelling and subtracting scattered light (be patient)')
    if not os.path.exists(config.options['masterorderdef'].value):
      raise taskError, 'Cannot find order definition reference frame : "%s"' % config.options['masterorderdef'].value

    if self.newmethod:


      messageLog.put('Using new subtraction method')
      orderdef = frameUtils.getorderdef(config.options['masterorderdef'].value)
      frameUtils.remove_bg(indata, orderdef=orderdef)

      # Write output to outframe (data object actually still points to input data)
      try:
	indata.writeto(outframe)
      except IOError:
	raise taskError('Cannot write output to %s' % outframe)

      indata.close()

    else:
    
      irafWrappers.scattering(inframe, outframe, config.options['masterorderdef'].value)

      # Update history header in the output frame
      outdata = pyfits.open(outframe)
      outdata[0].header.add_history('Modelled and subtracted scattered light')

      outdata.close()

################################################################################

class flatfield(task):
  """Divide by 2-dimensional normalized FLAT.
  """

  name         = 'flatfield'
  buttonText   = 'Divide by 2-D flat'
  suffix       = 'flat'
  prereq       = ['preproc']

  def run(self, inframe, outframe):

    # Open the 2D normalized combined FLAT
    try:
      flat = pyfits.open(config.options['masternormflat'].value)
    except IOError:
      raise taskError('Cannot open 2D-normalized combined FLAT : "%s"' % config.options['masternormflat'].value)

    # Open the input frame
    data = pyfits.open(inframe)

    if (data[0].data.shape != flat[0].data.shape):
      messageLog.put_warning('Data size (%i, %i)' % data[0].data.shape)
      messageLog.put_warning('Flat size (%i, %i)' % flat[0].data.shape)
      raise taskError('Incompatible sizes of object and normalized combined flat frame. Did you start binning?')

    # Perform the (2D) division
    messageLog.put('Divide by 2D-normalized flat')
    data[0].data = data[0].data / flat[0].data
    
    # Update history header
    data[0].header.add_history('Divided by normalized combined flat : %s' % os.path.basename(config.options['masternormflat'].value))

    # Save result in outframe
    try:
      data.writeto(outframe)
    except IOError:
      raise taskError('Cannot write output to %s' % outframe)

    # Clean up
    flat.close()
    data.close()


################################################################################

class plotcross(task):
  """Plot the cross-order profile. The profile is an average of the central 5
     columns of the frame.
  """

  name         = 'plotcross'
  buttonText   = 'Plot cross-order profile'
  prereq       = ['preproc']
  
  # This task produces no output file
  output       = 0


  def run(self, inframe, outframe):

    messageLog.put('Plotting cross-order profile')

    # Open input frame and determine data size
    data = pyfits.open(inframe)
    xsize = len(data[0].data[0, :])
    ysize = len(data[0].data[: ,0])
    
    # Select a slice in the center of the detector (10 pixels wide)
    centerslice = data[0].data[ysize/2-5:ysize/2+5, :]

    # Do the same for the other two slices
    lowslice  = data[0].data[ysize/4-5:ysize/4+5, :]
    highslice = data[0].data[3*ysize/4-5:3*ysize/4+5, :]

    # Collapse into 1D by taking row averate
    averrow  = numpy.mean(centerslice, axis=0)
    averlow  = numpy.mean(lowslice, axis=0)
    averhigh = numpy.mean(highslice, axis=0)

    # Plot the result
    title = "%s - Cross-order profile" % (os.path.basename(inframe,))
    plotUtils.plot(averrow, title=title)
    
    # Uncomment these lines to show more info about the scattered light
    # subtraction in other regions of the frame
    #
    # title = "%s - Cross-order profile (low slice) " % (os.path.basename(inframe,))
    # plotUtils.plot(averlow, title=title)
    # title = "%s - Cross-order profile (high slice)" % (os.path.basename(inframe,))
    # plotUtils.plot(averhigh, title=title)

################################################################################

class getspecshift(task):
  """Determine the spectrum shift by extracting the interlaced spectrum
     and comparing this to the interlaced wavelength reference frame.
  """

  # This task has to be performed before the (object) spectrum is extracted.
  # The result is stored in the FITS header, which is later accessed by the
  # task 'adjustwave'

  name         = 'getspecshift'
  buttonText   = 'Determine spectrum shift'
  suffix       = 'specshift'
  prereq       = ['preproc']

  def run(self, inframe, outframe):

    messageLog.put('Determining the wavelength shift of interlaced spectrum')

    # Remove outframe if it exists
    fileUtils.removefiles(outframe)

    # Create a temporary frame to hold intermediate output
    tempframe = os.path.join(os.path.dirname(inframe), 'tempframe.fits')
    #tempframe = inframe+'_interl.fits'
    fileUtils.removefiles(tempframe)

    # Check that an interlaced order definition exists
    if not os.path.exists(config.options['masterinterlacedorderdef'].value):
      raise taskError, 'Cannot open interlaced order definition reference: "%s"' % config.options['masterinterlacedorderdef'].value

    # Fix GAIN and RDNOISE headers if necessary
    frameUtils.fixheaders(inframe)

    # Extract the interlaced spectrum
    messageLog.put('Extracting the interlaced orders')
    irafWrappers.extractspec(inframe, tempframe, config.options['masterinterlacedorderdef'].value, clean='no')

    # Determine the (zeroth-order) vertical wavelength shift of the extracted
    # spectrum with respect to the default interlaced wavelength definition
    messageLog.put('Determining the wavelength solution in %s' % outframe)
    (pixshift, wvlshift, ratio) = irafWrappers.findwaveshift(tempframe, config.options['interlacedwaveref'].value)

    # Do not do any correction if less than 10% of the lines were identified
    if ratio < 0.1:
      messageLog.put_warning('Could only reidentify %s %% of the lines - ignoring shift' % (ratio*100))
      pixshift = 0.0
      wvlshift = 0.0
    else:
      messageLog.put('Reidentified %s %% of the lines' % (ratio*100))

    # Report to user
    messageLog.put('Found a shift in interlaced spectrum of %4.3f pixels' % pixshift)

    # Store the result in the frame header
    outdata = pyfits.open(inframe)
    outdata[0].header.update('WVLSHIFT', wvlshift, 'Interlaced wvl shift (fundamental order)')
    outdata[0].header.update('PIXSHIFT', pixshift, 'Interlaced pixel shift')
    outdata.writeto(outframe)
    outdata.close()

    # Clean up
    #irafWrappers.dispcor(tempframe, inframe+'_interl.fits', config.options['waveref'].value)
    fileUtils.removefiles(tempframe)

################################################################################

class extspec(task):
  """Extract individual orders from a two-dimensional frame using IRAF. The
     method used is the so-called 'optimal extraction'.
  """

  name         = 'extspec'
  buttonText   = 'Extract spectrum'
  suffix       = 'ext'
  prereq      = ['preproc']

  def run(self, inframe, outframe):

    # Remove outframe if it exists
    fileUtils.removefiles(outframe)

    # Make sure there is an order definition that can be used
    if not os.path.exists(config.options['masterorderdef'].value):
      raise taskError, 'Cannot open order definition reference : "%s"' % config.options['masterorderdef'].value

    # Fix GAIN and RDNOISE headers if necessary
    frameUtils.fixheaders(inframe)

    # Extract the orders using the PyRAF wrapper code
    messageLog.put('Extracting the orders')
#   irafWrappers.extractspec(inframe, outframe, config.options['masterorderdef'].value)
    irafWrappers.extractspec(inframe,
                             outframe,
			     config.options['masterorderdef'].value)


    # Update history card information
    outdata = pyfits.open(outframe, 'update')
    outdata[0].header.add_history('Extracted orders according to definition in %s' % os.path.basename(config.options['masterorderdef'].value))

    # Clean up
    outdata.close()


################################################################################

class blazecorr(task):
  """Divide the extracted orders by the blaze shape, producing flattened spectra.
     There is no additional continuum fitting performed. The blaze shape is taken from
     a pre-determined fit to the combined flat frame.
  """

  name         = 'blazecorr'
  buttonText   = 'Correct for blazeshape'
  suffix       = 'blzcorr'
  prereq      = ['extspec']

  def run(self, inframe, outframe):

    messageLog.put('Correcting for blaze shape')

    # Make sure there exists a blaze shape to correct with
    if not os.path.exists(config.options['blazeshape'].value):
      raise taskError, 'Cannot open fitted blaze shape : "%s"' % config.options['blazeshape'].value

    # Open the input frame and the blaze shape
    outdata = pyfits.open(inframe)
    blzshape = pyfits.open(config.options['blazeshape'].value)

    # And divide the two
    outdata[0].data = outdata[0].data / blzshape[0].data
    
    # Update history card information
    outdata[0].header.add_history('Divided out blaze shape using %s' % os.path.basename(config.options['blazeshape'].value))

    # Remove outframe if it exists and write to disk
    fileUtils.removefiles(outframe)
    outdata.writeto(outframe)

    # Clean up
    outdata.close()

################################################################################

class addwave(task):
  """Append wavelength information to the extracted orders. The wavelength definition
     is taken from the pre-determined solution of the wavelength reference frame.
  """

  name         = 'addwave'
  buttonText   = 'Add wavelengths'
  suffix       = 'wave'
  prereq      = ['extspec']

  def run(self, inframe, outframe):

    messageLog.put('Adding wavelength definition to spectrum')

    # Make sure there is a wavelength definition to apply in the first place
    if not os.path.exists(config.options['waveref'].value):
      raise taskError, 'Cannot find wavelength reference frame : "%s"' % config.options['waveref'].value

    # Remove outframe if it exists
    fileUtils.removefiles(outframe)

    # Use PyRAF wrapper code to transfer the wavelength definition to the
    # output frame
    irafWrappers.dispcor(inframe, outframe, config.options['waveref'].value)

    # Update history card information
    if os.path.exists(outframe):
      outdata = pyfits.open(outframe, 'update')
      outdata[0].header.add_history('Added wavelength definition using %s'
                          % os.path.basename(config.options['waveref'].value))
      outdata.close()
    else:
      raise taskError, 'Could not apply wavelength definition - no definition found'


################################################################################

class adjustwave(task):
  """Adjust the wavelength definition by correcting for possible shifts in the
     simultaneously observed interlaced calibration lamp spectrum.
  """

  name         = 'adjustwave'
  buttonText   = 'Adjust wavelengths'
  suffix       = 'adjwave'
  prereq      = ['getspecshift', 'addwave']

  def run(self, inframe, outframe):

    messageLog.put('Adjusting for shifts in wavelength definition')

    # Read input file and its wavelength definition
    indata  = pyfits.open(inframe)
    wavedef = waveUtils.getwavedef2D(indata)

    if wavedef is None:
#     print indata[0].header
      raise taskError, 'Could not read wavelength definition!'

    # Determine number of orders in frame
    norders = len(wavedef[:, 0])

    # Get determined shift from header card
    wvlshift = indata[0].header['WVLSHIFT']

    # Write the input data to the output file (YES, indeed before the
    # shifting routine is called. Essentially this is just copying the input
    # to the output file. Then the next operation takes care of operating on the
    # 'output' file
    indata.writeto(outframe)    
    indata.close()

    for order in range(norders):

      # Shift the order offset, one order at a time, where the shift in
      # wavelength is derived from the earlier determined fundamental
      # (zeroth) order wavelength shift.
      onumber = wavedef[order, 0]
      beamno = wavedef[order, 1]
      offset  = wavedef[order, 2]
      newoffset = offset - (wvlshift / beamno)
      messageLog.put("Shifting order %i from %f to %f" % (onumber, offset, newoffset), 7)
      irafWrappers.shiftorder(outframe, order=onumber, newoffset=newoffset)


################################################################################

class ordermerge(task):
  """Merge extracted orders into one 1-dimensional spectrum. Because the wavelength
     grid of the individual orders may be different, merging implies rebinning
     of the spectrum.
  """

  name        = 'ordermerge'
  buttonText  = 'Merge orders'
  suffix      = 'merge'
  prereq      = ['blazecorr', 'addwave']

  def run(self, inframe, outframe):

    # Use PyRAF wrapper to merge the orders. This is a time-consuming operation
    messageLog.put('Merging orders into one single spectrum (takes time - be patient)')

    # Output a general warning that this operation actually destroyed information!
    messageLog.put_warning('Flux not conserved after rebinning on new wavelength grid', 5)

    irafWrappers.mergeorders(inframe, outframe, config.options['waveref'].value)

    # Update history card information
    if os.path.exists(outframe):
      outdata = pyfits.open(outframe, 'update')
      outdata[0].header.add_history('Individual orders merged from %s' % outframe)
      outdata.close()


################################################################################
# The next group of tasks is used for calculating calibration (reference) frames
# These tasks are typically run on a one-off basis.
################################################################################

class sumbias(task):
  """Combine a set of BIAS frames into the combined BIAS frame. The currently
     implemented method is median filtering of the images.
  """

  name         = 'sumbias'
  buttonText   = 'Create combined BIAS'

  def run(self):
    
    # Get the user-defined list of bias frames
    framelist = config.options['biaslist'].value

    # Determine the number of bias frames to combine
    try:
      nframes = len(framelist[0])
    except IndexError:
      raise taskError('No BIAS frames defined')

    if not os.path.exists(os.path.dirname(config.options['masterbias'].value)):
      raise taskError, 'Directory of combined BIAS frame does not exist'
    if not config.options['masterbias'].value :
      raise taskError('Combined BIAS frame not defined')

    # Create a 3D dataset of the bias frames
    datacube = frameUtils.createdatacube(framelist, extension=config.options['mef_dataext'].value)
    if datacube is None: raise taskError('Could not process list of BIAS frames - Stopped')
    nframes = datacube.shape[0]

    nsetone = (nframes / 2)	# Important to divide integer by integer
    nsettwo = nframes - nsetone	# so we get (floored-)integer results

    firstmean  = numpy.sum(datacube[0:nsetone,  :, :], axis=0, dtype='float64') / float(nsetone)
    secondmean = numpy.sum(datacube[nsetone:,   :, :], axis=0, dtype='float64') / float(nsettwo)

    # Determine the standard deviation (unfiltered RON) of the difference between the frames
    messageLog.put('Calculating standard deviation of BIAS frames')
    unfilteredRON = numpy.std(firstmean - secondmean, dtype='float64') * numpy.sqrt((nsetone*nsettwo)/float(nframes))
    messageLog.put('Standard deviation is %f ADU' % unfilteredRON)
    del(firstmean, secondmean)

    # Determine the mean bias value
    messageLog.put('Calculating mean bias level')
    biaslevel = numpy.mean(datacube, dtype='float64')
    messageLog.put('Mean bias level is %f ADU' % biaslevel)

    # Set a highest and lowest value for plotting bias histogram
    blow  = int(biaslevel-5.0*unfilteredRON)
    bhigh = int(biaslevel+5.0*unfilteredRON)
    nbins = bhigh - blow

    # In principle, one should do pixel rejection on a pixel by pixel basis,
    # and not use a global average. Still, because one can expect the BIAS
    # to be reasonably flat, one can use the global average, and only reject
    # a handful of pixels too many.
    # If there would be a global structure (for example, a ramp) in the image,
    # one cannot use the global average, but should use 2-d (pixel-by-pixel)
    # averages.

    # Set highest and lowest acceptable values for bias pixels
    messageLog.put('Fixing pixels with values above and below 5 times standard deviation')
    rlow  = biaslevel-5.0*unfilteredRON
    rhigh = biaslevel+5.0*unfilteredRON

    # Loop to conserve memory (the following requires a lot of memory)
    for i in range(nframes):
    
      datacube[i,:,:] = numpy.select([datacube[i,:,:] < rlow, datacube[i,:,:] <= rhigh,  datacube[i,:,:] > rhigh], [biaslevel, datacube[i,:,:], biaslevel])

#     THIS IS THE PROPER WAY TO DO IT, BUT IS EXTREMELY RESOURCE INTENSIVE 
#     USING JUST THE OLD BIAS LEVEL WILL NOT INTRODUCE LARGE ERRORS	   
#     (ABOUT 0.01 ADU IN THE REPLACED PIXELS)				   
#     # (Re)calculate mean for 'good' pixels				   
#     messageLog.put('Calculating mean value for non-deviant pixels', 7)   
#     biaslevel = numpy.mean(datacube[noreplaceindex])  

      messageLog.put('Fixed deviant pixels in frame %i of %i' % (i+1, nframes))


    # Is this asking for too much?
    #RON = numpy.std(datacube)
    #messageLog.put('Read-out noise level is %f ADU' % RON)

    # Diagnostic plot...
    # Embraced with try/except for backward compatibilty with Python 2.2
    # and Numarray 0.9 (that doesn't have the histogram method)
    
    #try:
    #  # Determine and plot the bias pixel histogram
    #  histogram = numpy.histogram(datacube, nbins, range=(blow, bhigh))
    #  histvect  = numpy.arange(nbins) + blow
    #  plotUtils.plot(histogram, histvect, title='Histogram of (%i) adjusted BIAS frames - RON = %f' % (nframes, RON))
    #  # Clear memory
    #  del(histogram)
    #  del(histvec)
    #except: pass


    # Create an output frame
    messageLog.put('Creating combined BIAS frame')
    outdata = frameUtils.createframe(object='Combined BIAS')

    # Fill the output frame with data and calculate average
    messageLog.put('Averaging BIAS frames', 7)
    outdata[0].data = numpy.mean(datacube, axis=0, dtype='float64')
    outdata[0].header.add_history('Averaged %i frames to obtain combined BIAS' % nframes)

    # Clear memory again
    del(datacube)

    # Save combined bias to disk
    fileUtils.removefiles(config.options['masterbias'].value)
    messageLog.put('Saving combined BIAS to %s' % config.options['masterbias'].value, 5)
    outdata.writeto(config.options['masterbias'].value)
    
    # Clean up
    outdata.close()

################################################################################

class sumflat(task):
  """Combine a set of FLAT frames into the combined FLAT frame. The
     currently implemented method is averaging of the images.
  """

  name        = 'sumflat'
  buttonText  = 'Create combined FLAT'

  def run(self):
    
    framelist = config.options['flatlist'].value

    try:
      nframes = len(framelist)
    except IndexError:
      raise taskError('No FLAT frames defined')


    if not os.path.exists(os.path.dirname(config.options['masterflat'].value)):
      raise taskError, 'Directory of combined FLAT frame does not exist'
    if not config.options['masterflat'].value :
      raise taskError('Combined FLAT frame not defined')

    bias = pyfits.open(config.options['masterbias'].value)
    biasdata = bias[0].data

    messageLog.put('Creating combined FLAT frame')
    outdata = frameUtils.createframe(object='Combined FLAT')
    outdata[0].data = biasdata * 0.0


    # Loop over frames to conserve memory
    for i in range(nframes):

      messageLog.put('Processing frame %i of %i (%s)' % (i+1, nframes, os.path.basename(framelist[i])))

      try:
	inframe = pyfits.open(framelist[i])
      except IOError: raise taskError('Cannot open file: %s' % framelist[i])

      inframe = frameUtils.extractMEF(inframe, extension=config.options['mef_dataext'].value)

      try:
	frameUtils.flipframe(inframe, config.options['frameorientation'].value)
        frameUtils.clipframe(inframe)
      except:
	raise taskError('Could not flip or clip frame - something wrong with input data')

      indata = inframe[0].data
      
      if biasdata.shape != indata.shape:
	raise taskError('FLAT frame is not same size as combined BIAS frame')

      outdata[0].data = outdata[0].data + (indata - biasdata)


    bias.close()

    # Copy relevant FITS headers to combined FLAT
    flatframe = pyfits.open(framelist[0])
    flatframe = frameUtils.extractMEF(flatframe, extension=config.options['mef_dataext'].value)
    try:
      outdata[0].header.update('GAIN',    flatframe[0].header['GAIN'])
      outdata[0].header.update('RDNOISE', flatframe[0].header['RDNOISE'])
    except KeyError:
      outdata[0].header.update('GAIN',    1.0)
      outdata[0].header.update('RDNOISE', 0.0)
    flatframe.close()

    messageLog.put('Computing average')
    outdata[0].data = outdata[0].data / float(nframes)
    outdata[0].header.add_history('Combined images by averaging (%i files) ' % nframes)
    outdata[0].header.update('NSUMMED', nframes)

    fileUtils.removefiles(config.options['masterflat'].value)
    messageLog.put('Saving combined FLAT to %s' % config.options['masterflat'].value, 5)
    outdata.writeto(config.options['masterflat'].value)
    outdata.close()

################################################################################

class findord(task):
  """Interactively define and trace the location of spectral orders using a well-exposed
     order-definition frame. 
  """

  name        = 'findord'
  buttonText  = 'Find order locations'
  inthread    = 0
  
  interactive = 'yes'

  def run(self):

    messageLog.put('Reading order definition frame %s' % config.options['orderdef'].value, 5)
    try:
      indata = pyfits.open(config.options['orderdef'].value)
    except IOError: raise taskError('Cannot read order definition frame')

    if not os.path.exists(os.path.dirname(config.options['masterorderdef'].value)):
      raise taskError, 'Directory of order definition reference frame does not exist'
    if not config.options['masterorderdef'].value :
      raise taskError('Order definition reference frame not defined')

    indata = frameUtils.extractMEF(indata, extension=config.options['mef_dataext'].value)
    try:
      frameUtils.flipframe(indata, config.options['frameorientation'].value)
    except:
      raise taskError('Could not flip frame - something wrong with input data')

    frameUtils.clipframe(indata)

    messageLog.put('Performing a rough background subtraction')
    frameUtils.remove_bg(indata, nsteps=10, poly_order_y=5)

    fileUtils.removefiles(config.options['masterorderdef'].value)
    messageLog.put('Saving order definition reference to %s' % config.options['masterorderdef'].value, 5)
    indata.writeto(config.options['masterorderdef'].value)

    indata.close()
    
    messageLog.put('Finding the order locations')
    messageLog.put('Source file is %s' % config.options['masterorderdef'].value, 7)

    irafWrappers.findtraceord(config.options['masterorderdef'].value, interactive=self.interactive)

################################################################################

class findinterlacedord(task):
  """Interactively define and trace the location of interlaced spectral orders
     using a well-exposed order-definition frame. 
  """

  name        = 'findinterlacedord'
  buttonText  = 'Find interlaced order locs.'
  inthread    = 0

  interactive = 'yes'

  def run(self):

    messageLog.put('Reading order definition frame %s' % config.options['interlacedorderdef'].value, 5)
    try:
      indata = pyfits.open(config.options['interlacedorderdef'].value)
    except IOError: raise taskError('Cannot read interlaced order definition frame')

    if not os.path.exists(os.path.dirname(config.options['masterinterlacedorderdef'].value)):
      raise taskError, 'Directory of interlaced order definition reference (output) frame does not exist'
    if not config.options['masterinterlacedorderdef'].value :
      raise taskError('Interlaced order definition reference (output) frame not defined')

    indata = frameUtils.extractMEF(indata, extension=config.options['mef_dataext'].value)
    try:
      frameUtils.flipframe(indata, config.options['frameorientation'].value)
    except:
      raise taskError('Could not flip frame - something wrong with input data')
    frameUtils.clipframe(indata)

    messageLog.put('Performing a rough background subtraction')
    frameUtils.remove_bg(indata, nsteps=10)

    messageLog.put('Creating interlaced order definition reference frame')
    outdata = frameUtils.createframe(object="Interlaced order definition reference")
    outdata[0].data = indata[0].data
    outdata[0].header.add_history('Data originates from %s' % config.options['interlacedorderdef'].value)

    fileUtils.removefiles(config.options['masterinterlacedorderdef'].value)
    messageLog.put('Saving interlaced order definition reference to %s' % config.options['masterinterlacedorderdef'].value, 5)
    outdata.writeto(config.options['masterinterlacedorderdef'].value)

    indata.close()
    outdata.close()

    messageLog.put('Now shift the order definition in place')
    irafWrappers.findtraceord(config.options['masterinterlacedorderdef'].value,
                         reference=os.path.basename(config.options['masterorderdef'].value),
			 trace='no', interactive=self.interactive)

################################################################################

class plotorders(task):

  """Show (schematically) the recovered location of the spectral orders.
  """

  name        = 'plotorders'
  buttonText  = 'Plot order locations'

  def run(self):

    messageLog.put('Plotting order definition from %s' % config.options['masterorderdef'].value, 5)
    if not os.path.exists(config.options['masterorderdef'].value) :
      raise taskError, 'Order definition reference frame does not exist (%s)' % config.options['masterorderdef'].value

    data = pyfits.open(config.options['masterorderdef'].value)
    xsize = len(data[0].data[0, :])
    ysize = len(data[0].data[: ,0])
    data.close()

    try:
      orders = frameUtils.getorderdef(config.options['masterorderdef'].value)
    except TypeError: raise taskError, 'Order definition is not based on Chebyshev polynomials'
#    except: raise taskError, 'Could not read order definition file'

    if not os.path.exists(config.options['masterinterlacedorderdef'].value) :
      messageLog.put('Interlaced order definition reference frame does not exist (%s) - skipped' % config.options['masterorderdef'].value, 5)
      interlacedorders = None
    else:
      try:
	interlacedorders = frameUtils.getorderdef(config.options['masterinterlacedorderdef'].value)
      except TypeError: raise taskError, 'Interlaced order definition is not based on Chebyshev polynomials'
      except IOError:
	messageLog.put_warning('Could not read interlaced order definition file')
	interlacedorders = None

    plotUtils.plotorders(orders, interlacedorders=interlacedorders, xsize=xsize, ysize=ysize, title='Order locations')


################################################################################

class normflat(task):
  """Calculate a 2D-normalized combined FLAT by fitting the 2-dimensional shape of the
     orders. The normalized combined FLAT can be used to correct for local changes in 
     responsivity of the CCD, for example fringes. In addition to the normalized FLAT,
     a frame containing the fitted blaze shape is calculated. 
  """

  name        = 'normflat'
  buttonText  = 'Create normal. comb. FLAT'
  prereq      = ['sumflat', 'findord']

  newmethod = True

  def run(self):

    if not os.path.exists(config.options['masterflat'].value):
      raise taskError, 'Cannot find combined FLAT frame "%s"' % config.options['masterflat'].value
    if not os.path.exists(config.options['masterorderdef'].value):
      raise taskError, 'Cannot find order definition reference frame : "s"' % config.options['masterorderdef'].value
    if not os.path.exists(os.path.dirname(config.options['masternormflat'].value)):
      raise taskError, 'Directory of 2D-normalized combined flat does not exist'
    if not config.options['masternormflat'].value :
      raise taskError('2D-normalized combined flat frame not defined')
    if not os.path.exists(os.path.dirname(config.options['blazeshape'].value)):
      raise taskError, 'Directory of Fitted blaze shape frame does not exist'
    if not config.options['blazeshape'].value :
      raise taskError('Fitted blaze shape frame not defined')

    # Prepare all temporary files
    messageLog.put('Prepare temporary files', 7)
    tempframe1 = os.path.join(os.path.dirname(config.options['masternormflat'].value), 'tempframe1.fits')
    tempframe2 = os.path.join(os.path.dirname(config.options['masternormflat'].value), 'tempframe2.fits')
    tempframe3 = os.path.join(os.path.dirname(config.options['masternormflat'].value), 'tempframe3.fits')
    tempframe4 = os.path.join(os.path.dirname(config.options['masternormflat'].value), 'tempframe4.fits')
    tempframe5 = os.path.join(os.path.dirname(config.options['masternormflat'].value), 'tempframe5.fits')
    fileUtils.removefiles(tempframe1, tempframe2, tempframe3, tempframe4, tempframe5)

    # Obtained number of summed FLAT frames from header of masterflat
    masterflat = pyfits.open(config.options['masterflat'].value)
    n_summed_flats = masterflat[0].header['NSUMMED']
    del(masterflat)

    # Remove output file if it already exists
    fileUtils.removefiles(config.options['masternormflat'].value)

    # Create output object
    masternormflat = frameUtils.createframe('Normalized combined FLAT')


    if (self.newmethod):

      # Subtract scattered light using new method
      messageLog.put('Modelling and subtracting scattered light - new method')
      masterflat = pyfits.open(config.options['masterflat'].value)
      orderdef = frameUtils.getorderdef(config.options['masterorderdef'].value)
      frameUtils.remove_bg(masterflat, orderdef=orderdef)
      masterflat.writeto(tempframe1)
      masterflat.close()
    
    else:
    
      # Subtract scattered light and store in tempframe1 (non-MEF!)
      messageLog.put('Modelling and subtracting scattered light from combined flat')
      irafWrappers.scattering(config.options['masterflat'].value,
                              tempframe1,
			      config.options['masterorderdef'].value)

    # Fix GAIN and RDNOISE headers if necessary
    frameUtils.fixheaders(tempframe1)

    # The number 10000 in the statement below indicates that the combined FF
    # should have a S/N radio of at least 100. Pixels above this threshold will
    # be corrected for pixel-to-pixel varations, while those below are not.
    threshold = (10000 / n_summed_flats)

    # Calculate the 2D-normalized flat and store tempframe2 (non-MEF!)
    messageLog.put('Calculating 2D-normalized combined FLAT (may take time, be patient)')
    irafWrappers.normalize(tempframe1,
                           tempframe2,
			   config.options['masterorderdef'].value,
			   threshold)

    # Bail out if no outframe generated (IRAF exists without complaining...)
    if not os.path.exists(tempframe2):
      messageLog.put_warning('The order extraction failed for at least one order')
      messageLog.put_warning('Try to review you order definition, ')
      messageLog.put_warning('and maybe remove orders close to the CCD edges.')
      raise taskError('Order extraction failed')

    # Copy the contents of tempframe2 to the combined normalized flat
    tf2 = pyfits.open(tempframe2)
    masternormflat[0].data = tf2[0].data
    # And save the master normalizrd flat
    masternormflat.writeto(config.options['masternormflat'].value)
    tf2.close()
    
    # Remove tempframe2, because we'll need to create a new one
    fileUtils.removefiles(tempframe2)

    # Determine the remaining blaze variations by dividing flat by normalized flat
    messageLog.put('Calculate remaining blaze function', 7)
    tf1 = pyfits.open(tempframe1)
    tf1[0].data = tf1[0].data / masternormflat[0].data
    tf1.writeto(tempframe2)
    tf1.close()

    # Remove output file if i already exists
    fileUtils.removefiles(config.options['blazeshape'].value)

    # Fix GAIN and RDNOISE headers if necessary
    frameUtils.fixheaders(tempframe2)

    # And extract the 1-dimensional blaze functions (result in tempframe2, non-MEF)
    messageLog.put('Extracting blaze shape')
    irafWrappers.extractspec(tempframe2, tempframe3, config.options['masterorderdef'].value)

    # Bail out if no outframe generated (IRAF exists without complaining...)
    if not os.path.exists(tempframe3):
      messageLog.put_warning('The order extraction failed for at least one order')
      messageLog.put_warning('Try to review you order definition, ')
      messageLog.put_warning('and maybe remove orders close to the CCD edges.')
      raise taskError('Order extraction failed')

    # Fit a smooth curve to the extracted orders
    irafWrappers.fit1d(tempframe3, tempframe4)

    # 'Create' the frame containing the fitted blaze shape
    blzim = pyfits.open(tempframe3)

    # Perform light Gaussian filtering of the blaze shape to exclude the FF-noise
    lightfilter = blzim[0].data.copy()

    kernel = numpy.arange(31)/3. - 5
    kernel = 1 / (3 * numpy.sqrt(numpy.pi)) * numpy.exp(-kernel*kernel)
    
    npix = blzim[0].data.shape[1]
    for i in numpy.arange(blzim[0].data.shape[0]):
      lightfilter[i, :]  = numpy.convolve(blzim[0].data[i, :], kernel, mode='same')
      lightfilter[i, 0:5] = blzim[0].data[i, 0:5]
      lightfilter[i, npix-5:npix] = blzim[0].data[i, npix-5:npix]

    # Fit a smooth curve to the blaze function
    fitim = pyfits.open(tempframe4)
    strongfilter = fitim[0].data

    # Use the smooth curve for low S/N data, and the filtered blaze for high S/N data
    # (same threshold as used in FF-normalization)
    blzim[0].data = numpy.choose(blzim[0].data < threshold, (lightfilter, strongfilter))

    # Put low values of the FF to 1, to avoid low (=bad) signal attenuation
    blzim[0].data = numpy.choose(blzim[0].data < 1, (blzim[0].data, 1))

    # Write the blaze image to disk
    messageLog.put('Saving extracted blaze shape to %s' % config.options['blazeshape'].value)
    blzim.writeto(config.options['blazeshape'].value)
    blzim.close()
    fitim.close()

    # (Master-normflat is still open)
    messageLog.put('Saving 2D-normalized combined FLAT to %s' % config.options['masternormflat'].value)
    masternormflat.close()

    # Cleanup
    fileUtils.removefiles(tempframe1, tempframe2, tempframe3, tempframe4)


################################################################################

class wavecal(task):
  """Interactively determine a 2-dimensional wavelength definition from a
     well-exposed wavelength definition frame, normally an ThAr frame.
  """

  name        = 'wavecal'
  buttonText  = 'Find wavelength solution'
  prereq      = ['findord']
  inthread    = 0

  interactive = 'yes'

  def run(self):

    if not os.path.exists(config.options['wavedef'].value):
      raise taskError, 'Cannot find wavelength definition frame : "%s"' % config.options['wavedef'].value
    if not os.path.exists(os.path.dirname(config.options['waveref'].value)):
      raise taskError, 'Directory of wavelength reference frame does not exist'
    if not os.path.exists(config.options['masterorderdef'].value):
      raise taskError, 'Cannot find order definition reference frame : "%s"' % config.options['masterorderdef'].value

    messageLog.put('Prepare temporary file', 7)
    tempframe1 = os.path.join(os.path.dirname(config.options['waveref'].value), 'tempframe1.fits')
    fileUtils.removefiles(tempframe1)

    indata = pyfits.open(config.options['wavedef'].value)
    indata = frameUtils.extractMEF(indata, extension=config.options['mef_dataext'].value)

    try:
      frameUtils.flipframe(indata, config.options['frameorientation'].value)
    except:
      raise taskError('Could not flip frame - something wrong with input data')
    frameUtils.clipframe(indata)

    messageLog.put('Creating wavelength reference frame')
    outdata = frameUtils.createframe(object="Wavelength reference")
    outdata[0].data = indata[0].data
    outdata[0].header.add_history('Data originates from %s' % config.options['wavedef'].value)

    try:
      bias = pyfits.open(config.options['masterbias'].value)
      outdata[0].data = outdata[0].data - bias[0].data
      bias.close()
    except IOError:
      messageLog.put_warning('Could not subtract combined BIAS - trying to continue')

    fileUtils.removefiles(config.options['waveref'].value)
    outdata.writeto(config.options['waveref'].value)

    indata.close()
    outdata.close()

    messageLog.put('Extracting orders from wavelength reference frame')
    irafWrappers.extractwave(config.options['waveref'].value,
                             tempframe1,
			     config.options['masterorderdef'].value)


    os.rename(tempframe1, config.options['waveref'].value)      

    messageLog.put('Determining the wavelength solution')
    if os.path.exists(config.options['masterwaveref'].value) :
      if (config.options['masterwaveref'].value != config.options['waveref'].value) :
        messageLog.put('Re-using wavelength solution from %s' % config.options['masterwaveref'].value)
        irafWrappers.refindwavesol(config.options['waveref'].value, config.options['masterwaveref'].value)
        if self.interactive == 'yes':
	  irafWrappers.findwavesol(config.options['waveref'].value, keep_old_solution=1)
      else :
        messageLog.put_warning('Action would destroy master wavelength reference data!')
        raise taskError, 'Filenames of "master wavelength reference" and "wavelength reference" are identical'
    else :
      irafWrappers.findwavesol(config.options['waveref'].value)
#   irafWrappers.refindwavesol(config.options['waveref'].value, config.options['waveref'].value)

    # Assign the wavelength reference to the output reference frame
#   HCS 061231: removed due to "error: overwrite existing image ?"
#   irafWrappers.dispcor(config.options['waveref'].value, config.options['waveref'].value, config.options['waveref'].value)

    fileUtils.removefiles(tempframe1)


################################################################################

class interlacedwavecal(task):
  """Interactively determine a 2-dimensional wavelength definition from a
     well-exposed interlaced wavelength definition frame, normally an ThAr frame.
  """

  name        = 'interlacedwavecal'
  buttonText  = 'Find interlaced wavel. sol.'
  prereq      = ['findinterlacedord']
  inthread    = 0

  def run(self):

    if not os.path.exists(config.options['interlacedwavedef'].value):
      raise taskError, 'Cannot find interlaced wavelength definition frame : "%s"' % config.options['interlacedwavedef'].value
    if not os.path.exists(os.path.dirname(config.options['interlacedwaveref'].value)):
      raise taskError, 'Directory of interlaced wavelength reference frame does not exist'
    if not os.path.exists(config.options['masterorderdef'].value):
      raise taskError, 'Cannot find interlaced order definition reference frame : "%s"' % config.options['masterinterlacedorderdef'].value

    messageLog.put('Prepare temporary file', 7)
    tempframe1 = os.path.join(os.path.dirname(config.options['interlacedwaveref'].value), 'tempframe1.fits')
    fileUtils.removefiles(tempframe1)

    indata = pyfits.open(config.options['interlacedwavedef'].value)
    indata = frameUtils.extractMEF(indata, extension=config.options['mef_dataext'].value)

    try:
      frameUtils.flipframe(indata, config.options['frameorientation'].value)
    except:
      raise taskError('Could not flip frame - something wrong with input data')
    frameUtils.clipframe(indata)

    messageLog.put('Creating interlaced wavelength reference frame')
    outdata = frameUtils.createframe(object="Interlaced wavelength reference")
    outdata[0].data = indata[0].data
    outdata[0].header.add_history('Data originates from %s' % config.options['interlacedwavedef'].value)

    bias = pyfits.open(config.options['masterbias'].value)
    outdata[0].data = outdata[0].data - bias[0].data

    fileUtils.removefiles(config.options['interlacedwaveref'].value)
    outdata.writeto(config.options['interlacedwaveref'].value)

    bias.close()
    indata.close()
    outdata.close()

    messageLog.put('Extracting orders from interlaced wavelength reference frame')
    irafWrappers.extractwave(config.options['interlacedwaveref'].value,
                             tempframe1,
			     config.options['masterinterlacedorderdef'].value)


    os.rename(tempframe1, config.options['interlacedwaveref'].value)      

    messageLog.put('Determining the interlaced wavelength solution')
#   irafWrappers.findwavesol(config.options['interlacedwaveref'].value)
    irafWrappers.refindwavesol(config.options['interlacedwaveref'].value, config.options['waveref'].value)

    fileUtils.removefiles(tempframe1)

