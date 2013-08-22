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
# frameUtils.py
#
# Last update 30/03/2009
#
################################################################################

"""
   Frame-related utility routine(s) that make life easier.
"""

# Import external modules

import os
import numpy
import pyfits
import config
import messageLog

import irafDbNames

################################################################################

# Helper functions :

def _flipaxis(arr, axis=0):

  # Might be done more efficiently with numpy.fliplr() and flipud()
  return arr.take(numpy.arange(arr.shape[axis], 0, -1) - 1, axis=axis)


def _chebev(a, b, c, x):
  """
     Calculate Chebyshev polynomial expansion with coefficients c for points
     (or vector of points) x.
     
     Adopted from Numerical Recipes.
  """

  # a = minimum of range
  # b = maximum of range
  # c = Chebyshev coefs
  # x = vector

  m = len(c)
    
  d  = 0.0
  dd = 0.0
  
  y = (2.0 * x - a - b) / (b - a)
  y2 = 2.0 * y
  
  rng = range(1, m)
  rng.reverse()
  
  for j in rng :
  
    sv = d
    d = y2 * d - dd + c[j]
    dd = sv
  
  # Note that the NR routine evaluates the Chebyshev polynomial sum - c[0]/2
  # (but... why???!!), so the following line is different from the NR code.
  # This gives correct IRAF order definitions!
  result = y * d - dd + c[0]
  
  return result


def sigclip(xarray, yarray, poly_order=5, sigma=3, max_iter=5):

  # Perform sigma-clipping on an array
  if xarray.size != yarray.size: raise ValueError;

  niter = 0
  while 1:

    coefs = numpy.polyfit(xarray, yarray, poly_order)
    polyfitfunc = numpy.poly1d(coefs)
    yfit = polyfitfunc(xarray)

    sigma = (yarray - yfit) / numpy.std(yarray - yfit)
    rejected = numpy.extract(sigma > 3, xarray)
    xarray   = numpy.extract(sigma < 3, xarray)
    yarray   = numpy.extract(sigma < 3, yarray)

    niter = niter + 1
    if niter == max_iter or rejected.size == 0: break

  return (xarray, coefs, polyfitfunc)

################################################################################


def WCSvector(frame, extension=0):

  """
     Return vectors of the WCS (World Coordinate System) x and y pixel
     coordinates in this frame. Particularly useful for binned or windowed
     frames.
  """

  # 'frame' is assumed to be an object created by pyfits.open(). Failure
  # to read any of the header keywords generates a KeyError exception, which
  # is supposed to be handled by the calling routine. However, for properly
  # defined frames, exceptions should not need to occur

  # Get the number of pixels for the two axes from the FITS headers
  # No try/except needed, because these headers MUST be defined.
  try:
    xsize = frame[extension].header['NAXIS1']
    ysize = frame[extension].header['NAXIS2']
  except KeyError:
    messageLog.put_error('NAXIS1 or NAXIS2 keyword not defined. Are you working with the correct extension?')
    return

  # Get the remaining vector values, in a truly failsafe way.

  try:
    x0 = frame[extension].header['CRVAL1']
  except KeyError:
    x0 = 1

  try:
    y0 = frame[extension].header['CRVAL2']
  except KeyError:
    y0 = 1
  
  try:
    xref = frame[extension].header['CRPIX1']
  except KeyError:
    xref = 1

  try:
    yref = frame[extension].header['CRPIX2']
  except KeyError:
    yref = 1

  try:
    xstep, dummy = frame[extension].header['CCDSUM'].split()
    xstep = int(xstep)
  except:
    raise
    try:
      xstep = frame[extension].header['CDELT1']
    except KeyError:
      xstep = 1

  try:
    dummy, ystep = frame[extension].header['CCDSUM'].split()
    ystep = int(ystep)
  except:
    raise
    try:
      ystep = frame[extension].header['CDELT2']
    except KeyError:
      ystep = 1
  

  # Construct the actual vectors with x and y pixel coordinates.
  # numpy.arange(n) gives a vector with n numbered elements
  xvec = (x0 + ( (numpy.arange(xsize) - xref + 1) * xstep))
  yvec = (y0 + ( (numpy.arange(ysize) - yref + 1) * ystep))

  # Return a tuple with the two vectors
  return (xvec, yvec)


################################################################################


def createframe(object="", xbin=1, ybin=1, nodata=0, extension=0):

  """
     Create a new frame from scratch using the default frame size and a
     given binning. Returns a PyFITS object.
  """

  # Determine the default frame size
  xmin = config.options['x_start'].value
  xmax = config.options['x_end'  ].value
  ymin = config.options['y_start'].value
  ymax = config.options['y_end'  ].value

  xsize = (xmax - xmin)/xbin + 1
  ysize = (ymax - ymin)/ybin + 1

  messageLog.put('Creating new frame: %s (%i, %i) pixels' % (object, xsize, ysize), 7)


  # Create a list of HDUs (Header Data Units).
  newframe = pyfits.HDUList()

  # ...and create the primary HDU
  primhdu  = pyfits.PrimaryHDU()

  # Set OBJECT header
  primhdu.header.update('OBJECT', object)

  # Attach the primary HDU to the frame
  newframe.append(primhdu)

  # Return if frame should not contain data block
  if nodata: return newframe

  # Also create an image HDU
  if (extension != 0):

    for extno in range(1, extension+1) :

      imhdu = pyfits.ImageHDU()

      # Attach the image HDU to the frame
      newframe.append(imhdu)

      # Create EXTNAME header (to satisfy IRAF)
      newframe[extno].header.update('EXTNAME', 'im%i' % extno)


  # Fill the data unit with an empty frame of requested size
  newframe[extension].data = numpy.zeros((ysize, xsize), dtype=numpy.float64)

  # Update header values accordingly
  newframe[extension].header.update('CRVAL1', xmin)
  newframe[extension].header.update('CRPIX1', 1)
  newframe[extension].header.update('CDELT1', xbin)
  newframe[extension].header.update('CTYPE1', 'PIXEL')

  newframe[extension].header.update('CRVAL2', ymin)
  newframe[extension].header.update('CRPIX2', 1)
  newframe[extension].header.update('CDELT2', ybin)
  newframe[extension].header.update('CTYPE2', 'PIXEL')

  newframe[extension].update_header()


  # Return the frame object
  return newframe


################################################################################


def flipframe(frame, direction=0, extension=0):

  # Retrieved the x and y coordinate vectors
  (xvec, yvec) = WCSvector(frame, extension)


  # MMM... MAYBE NOT TRUE... 
  # Comment : flipping the xvec and yvec below is not really
  # necessary, because the xnewmin and ynewmin will still be the
  # lowest pixel number in the (flipped) xvec and yvec. However,
  # exchanging the xvec and yvec is _very_ relevant.

  # NB: xaxis is axis nr 1, and yaxis is axis nr 0

  flippeddata = frame[extension].data

  # Rotate 'n' times 90 degrees
  # x ->  x, y ->  y
  if   direction == 0:
    flippeddata = flippeddata

  # x -> -y, y ->  x
  elif direction == 1:
    flippeddata = _flipaxis(flippeddata, axis=1)
    flippeddata.transpose()
    yvec = _flipaxis(yvec, axis=0)

  # x -> -x, y -> -y
  elif direction == 2:
    flippeddata = _flipaxis(flippeddata, axis=1)
    flippeddata = _flipaxis(flippeddata, axis=0)
    xvec = _flipaxis(xvec, axis=0)
    yvec = _flipaxis(yvec, axis=0)

  # x ->  y, y -> -x
  elif direction == 3:
    flippeddata = _flipaxis(flippeddata, axis=0)
    flippeddata.transpose()
    xvecnew = yvec
    yvec = _flipaxis(xvec, axis=0)
    xvec = xvecnew
    del(xvecnew)

  # Rotate 'n' times 90 degrees, and transpose
  # x ->  y, y ->  x
  elif direction == 4:
    flippeddata.transpose()
    xvecnew = yvec
    yvec = xvec
    xvec = xvecnew
    del(xvecnew)

  # x -> -x, y ->  y
  elif direction == 5:
    flippeddata = _flipaxis(flippeddata, axis=1)
    xvec = _flipaxis(xvec, axis=0)

  # x -> -y, y -> -x
  elif direction == 6:
    flippeddata = _flipaxis(flippeddata, axis=1)
    flippeddata = _flipaxis(flippeddata, axis=0)
    flippeddata.transpose()

  # x ->  x, y -> -y
  elif direction == 7:
    flippeddata = _flipaxis(flippeddata, axis=0)
    yvec = _flipaxis(yvec, axis=0)

  else:
    messageLog.put_warning('Frame orientation parameter out of range - ignored')


  # Determine 'new' vector minima and stepsize
  xnewbin = abs(xvec[2] - xvec[1])
  ynewbin = abs(yvec[2] - yvec[1])

  xnewmin = xvec.min()
  ynewmin = yvec.min()


  # Update header values accordingly
  frame[extension].header.update('CRVAL1', xnewmin)
  frame[extension].header.update('CRPIX1', 1)
  frame[extension].header.update('CDELT1', xnewbin)
  frame[extension].header.update('CTYPE1', 'PIXEL')

  frame[extension].header.update('CRVAL2', ynewmin)
  frame[extension].header.update('CRPIX2', 1)
  frame[extension].header.update('CDELT2', ynewbin)
  frame[extension].header.update('CTPYE2', 'PIXEL')

  # Store the rotated and flipped data in fits object
  frame[extension].data = flippeddata

  # And adjust the image headers to reflect the data
  frame[extension].update_header()


################################################################################


def clipframe(frame, extension=0):

  """
     Remove under- and overscan regions from a frame using the WCS
     pixel coordinates. If a frame is already within the limits,
     leave it unchanged (no padding).
  """

  # 'frame' is assumed to be an object created by pyfits.open().

  # Get the maximum allowed frame size from the current configuration
  xabsmin = config.options['x_start'].value
  xabsmax = config.options['x_end'  ].value
  yabsmin = config.options['y_start'].value
  yabsmax = config.options['y_end'  ].value
 
 
  # Retrieved the x and y coordinate vectors
  (xvec, yvec) = WCSvector(frame, extension)


  # Create new x and y vectors of the pixels coordinates that should be
  # in the output frame. Use boolean logic to select the pixels within bounds
  # Here, 'newxvec' and 'newyvec' contain only 0 or 1.
  newxvec = (xvec >= xabsmin) * (xvec <= xabsmax)
  newyvec = (yvec >= yabsmin) * (yvec <= yabsmax)
 

  # Get a subset of the data for the non-zero elements of the new x and y
  # vectors. Then, remove all non-zero elements of the new x and y vectors.
  #
  # The [0] after nonzero() is there because nonzero() returns a tuple
  # that contains an array (don't know why - bug?)
  #
  # Do not change the order of the statements below! Slicing is done twice.

  # Get the subset of the data along the x axis
  subset  = frame[extension].data.take(newxvec.nonzero()[0], axis=1)
  # Make 'newxvec' contain the actual pixel positions, rather than 0 or 1
  newxvec = xvec.take(newxvec.nonzero()[0])
  # Force the new vector to be one-dimensional (take returns 2D)
  newxvec.shape = (-1,)

  # Idem for y
  subset  = subset.take(newyvec.nonzero()[0], axis=0)
  newyvec = yvec.take(  newyvec.nonzero()[0])
  newyvec.shape = (-1,)

  # 'subset' only contains pointers (!) to (part of) the data in 'data'.
  # Therefore, make a real copy of the data, that can be stored in an output
  # file
  newdata = subset.copy()

  # Adjust the 2D shape of newdata
  newdata.shape = (len(newyvec), len(newxvec))

  # Determine if there is binning from the pixel positions
  newxbin = newxvec[2] - newxvec[1]
  newybin = newyvec[2] - newyvec[1]

  # Store the sliced data in the original frame
  frame[extension].data = newdata

  # And update the headers accordingly
  frame[extension].header.update('CRVAL1', newxvec[0])
  frame[extension].header.update('CRPIX1', 1)
  frame[extension].header.update('CDELT1', newxbin)
  frame[extension].header.update('CTYPE1', 'PIXEL')

  frame[extension].header.update('CRVAL2', newyvec[0])
  frame[extension].header.update('CRPIX2', 1)
  frame[extension].header.update('CDELT2', newybin)
  frame[extension].header.update('CTYPE2', 'PIXEL')

  frame[extension].update_header()


################################################################################


def padframe(frame, extension=0):

  """
     Extend an existing frame to the default image size. Returns 0 if the frame
     is already larger or equal to the default size, and 1 if the frame was
     padded.
  """

  # Get pixel numbering of this frame (from WCS)
  (xvec, yvec) = WCSvector(frame, extension)

  # Get current binning factor
  xbin = xvec[2] - xvec[1]
  ybin = yvec[2] - yvec[1]

  # Determine current frame size
  xmin = min(xvec)
  xmax = max(xvec)
  ymin = min(yvec)
  ymax = max(yvec)
  xsize = (xmax - xmin)/xbin + 1
  ysize = (ymax - ymin)/ybin + 1

  # Determine the maximum allowed frame size
  xnewmin = config.options['x_start'].value
  xnewmax = config.options['x_end'  ].value
  ynewmin = config.options['y_start'].value
  ynewmax = config.options['y_end'  ].value

  # And calculate the expected frame size from these numbers
  # (Use same algorithm as in createframe and createcube)
  xnewsize = (xnewmax - xnewmin)/xbin + 1
  ynewsize = (ynewmax - ynewmin)/ybin + 1


  # Check if padding is really necessary. If not, return with exit code 0
  if ( xsize == xnewsize and ysize == ynewsize ):
    return 0


  # Create new x and y index vectors
  # Note: arange starts with 0
  newxvec = numpy.arange(xnewsize)*xbin + xnewmin
  newyvec = numpy.arange(ynewsize)*ybin + ynewmin


  # Find which part of the new frame should contain the data
  # from the new frame (because old frame is smaller than new frame)
  #
  # The [0] after nonzero is there because nonzero return a tuple
  # that contains an array (don't know why - bug?)

  selxvec = (newxvec >= xmin) * (newxvec <= xmax)
  selxvec = selxvec.nonzero()[0]

  selyvec = (newyvec >= ymin) * (newyvec <= ymax)
  selyvec = selyvec.nonzero()[0]


  # Create a new zero-filled full-sized frame
  newdata = numpy.zeros((ynewsize, xnewsize), dtype=numpy.float64)

  # Transfer the data into the new frame. Note that the upper boundary of a
  # slice in a Python array is exclusive! Therefore +1...
  newdata[min(selyvec):max(selyvec)+1,
          min(selxvec):max(selxvec)+1] = frame[extension].data


  # Store the sliced data in the original frame
  frame[extension].data = newdata

  # Update header values accordingly
  frame[extension].header.update('CRVAL1', xnewmin)
  frame[extension].header.update('CRPIX1', 1)
  frame[extension].header.update('CDELT1', xbin)
  frame[extension].header.update('CTYPE1', 'PIXEL')

  frame[extension].header.update('CRVAL2', ynewmin)
  frame[extension].header.update('CRPIX2', 1)
  frame[extension].header.update('CDELT2', ybin)
  frame[extension].header.update('CTPYE2', 'PIXEL')

  frame[extension].update_header()

  # Return that padding was succesful
  return 1


################################################################################


def createdatacube(framelist, extension=0):

  """
     Read data for a list of images, and store this in a 3-dimensional array
     (numpy type). A 3-D array is practical for performing efficient
     calculations on a large number of images (for example, averaging).
  """

  # Determine the number of frames in the list (will be z-size of cube)
  nframes = len(framelist)

  # Counter to keep track of position in the datacube
  n = 0

  # Read first frame in list to determine frame size (x- and y-size of cube)
  messageLog.put('Reading frame %i of %i (%s)' % (n+1, nframes, os.path.basename(framelist[0])))
  try:
    image = pyfits.open(framelist[0])
  except IOError:
    messageLog.put_error('Cannot open %s' % framelist[0])
    return None

  image = extractMEF(image, extension=extension)

  # THIS IS NOT THE BEST PLACE TO PUT THESE, BUT IT SHOULD GO SOMEWHERE
  flipframe(image, config.options['frameorientation'].value)
  clipframe(image)
  data = image[0].data
  framesize = data.shape

  # Create am empty 3D dataset and insert first frame
  # (Use .copy() to keep the data in memory after the image is closed)
  datacube = numpy.zeros((nframes, framesize[0], framesize[1]), data.dtype)
  datacube[n, :, :] = data.copy()

  # Close the file
  image.close()

  # Next frame...
  n = n + 1

  # Loop over remaining frames (start from index 1, not 0)
  for frame in framelist[1:]:

    messageLog.put('Reading frame %i of %i (%s)' % (n+1, nframes, os.path.basename(frame)))

    # Read a frame
    try:
      image = pyfits.open(frame)
    except IOError:
      messageLog.put_error('Cannot open %s' % frame)
      return None


    # THIS IS NOT THE BEST PLACE TO PUT THESE, BUT IT SHOULD GO SOMEWHERE
    image = extractMEF(image, extension=extension)
    flipframe(image, config.options['frameorientation'].value)
    clipframe(image)

    data = image[0].data

    # Check that the frame size is compatible with first frame
    if data.shape != framesize:
      messageLog.put_error('Not all frames have identical sizes')
      return None

    # Insert the frame in the cube
    # (Use .copy() to keep the data in memory after the image is closed)
    datacube[n, :, :] = data.copy()

    # Close the file
    image.close()

    # Next frame...
    n = n + 1


  # Return the cube
  return datacube


################################################################################


def getorderdef(frame):

  """
     Not so nice routine to read the x and y positions of spectral orders from
     an IRAF definition file, and store these in a dictionary object.
  """

  # Separate the directory and the file name
#  base, frame = os.path.split(frame)

  # Split off the filename extension
#  frame, ext = os.path.splitext(frame)

  # Construct the filename of the IRAF order definition data file that
  # belongs to 'frame'
  databasefilename = irafDbNames.apname(frame)

  # Initialize parameters
  norders = 0
  orderarray = []

  # Open the definition file for reading
  fh = open(databasefilename)


  try:

    # Loop until broken
    while 1 :

      # Read a single line
      line = fh.readline()
      
      # If this was unsuccessful (EOF) then break the loop
      if not line: break

      # Split the line into words (using whitespace)
      words = line.split()

      # Empty line?
      if len(words) == 0 : continue

      # Does this line define the start of a new order definition?
      if words[0] == 'begin'  : 
				# Increase total number of orders by 1
				norders = norders + 1
				# Append empty order definition dictionary object
				orderarray.append( {} )

      # Does this line define the center position of the order?
      if words[0] == 'center' : 
    				# Store position of the center pixel and it's
				# corresponding (x-)value in keywords in the
				# dictionary 'orderarray'. Python array numbering
				# starts from 0, so therefore 'norders-1'.
				orderarray[norders-1]['centerpix'] = float(words[2])
				orderarray[norders-1]['centerval'] = float(words[1])

      # Is this line the first of several lines that define the coefficients
      # of the order location curve?
      if words[0] == 'curve' :  
				# Get the number of coefficients of this order
				ncoefs = int(words[1])
				orderarray[norders-1]['ncoefs'] = ncoefs
				
				# Create an array of coefficients
    				coefs = []

    				# Read the coefficients from the file
				for j in range(ncoefs):
			          value = fh.readline()
				  # Store in the coefficient array
			          coefs.append(float(value.strip()))

				# Store the coefficient array in 'orderarray'
				orderarray[norders-1]['coefs'] = coefs


  # Cope with i/o-errors (not implemented)
  except IOError: raise

  # Close the input file
  fh.close
  
  # Return the set of order definition parameters
  return orderarray

################################################################################

def getorderxy(order, npoints=50):

  """
     Simple hack to convert the IRAF order definition to x and y positions
     on the detector. 'npoints' determines how coarse the sampling along the
     orders is allowed to be.
  """

  # This routine can only handle order definitions using Chebychev polynoms!
  if order['coefs'][0] != 1.0 : 
    raise TypeError, 'Wrong type of polynom used in order fitting (not Chebychev)'

  # Get min and max pixel values
  ymin = order['coefs'][2]
  ymax = order['coefs'][3]

  # Determine the stepsize needed to put in 'npoints' points.
  stride = (ymax - ymin) / (npoints-1)

  # Get the central pixel offset
  centerpix = int((order['centerpix'] - ymin)/stride)
  centerval = order['centerval']

  # Take the coefficients...
  chebcoef = order['coefs'][4:]

  # ...and use these coefficients to determine the x-positions from the
  # chebychev polynom.
  yvec = numpy.arange(ymin, ymax, stride)
  xvec = _chebev(ymin, ymax, chebcoef, yvec)

  # Move the x-vector to the correct offset of the central pixel
  try:
    xvec = xvec - xvec[centerpix] + centerval
  except IndexError:
    print "Error in order definition?"
    pass

  # return two sets of arrays that can be plotted
  return (xvec, yvec)

################################################################################

def extractMEF(frame, extension=1):

  """
     Extract a FITS Header Data Unit (HDU) from a Multi-Extension FITS (MEF)
     file and return non-MEF data
  """

  if (extension == 0) : return frame

  try: frame[extension]
  except IndexError: return None

  # (NB: frame is an existing PyFITS object)

  # Create a list of HDUs (Header Data Units).
  newframe = pyfits.HDUList()

  # Create and attach a primary HDU to the frame
  newframe.append(pyfits.PrimaryHDU())

  # Copy headers from primary data unit
  for (headername, headervalue) in frame[0].header.items() :
    if (headername == "COMMENT") :
      newframe[0].header.add_comment(headervalue)
      continue
    if (headername == "HISTORY") :
      newframe[0].header.add_history(headervalue)
      continue
    newframe[0].header.update(headername, headervalue)

  # Copy headers from selected image data unit
  for (headername, headervalue) in frame[extension].header.items() :
    if (headername == "COMMENT") :
      newframe[0].header.add_comment(headervalue)
      continue
    if (headername == "HISTORY") :
      newframe[0].header.add_history(headervalue)
      continue
    newframe[0].header.update(headername, headervalue)

  # Remove old headers that may contain incorrect values or that are
  # not allowed in the primary HDU
  del(newframe[0].header['NAXIS'])
  del(newframe[0].header['NAXIS1'])
  del(newframe[0].header['NAXIS2'])
  del(newframe[0].header['BZERO'])
  del(newframe[0].header['BSCALE'])
  del(newframe[0].header['BITPIX'])
  del(newframe[0].header['EXTNAME'])
  del(newframe[0].header['XTENSION'])
  del(newframe[0].header['PCOUNT'])
  del(newframe[0].header['GCOUNT'])

  # Attach the new dataset to the primary HDU
  newframe[0].data = frame[extension].data.copy()
  
  # And fix the FITS output. This means reintroducing some of the keywords that
  # were just deleted, but now with values that agree with the new dataset
  newframe.verify('silentfix')
  newframe[0].update_header()

  return newframe

################################################################################

def fixheaders(frame):

  # Make sure header cards GAIN and RDNOISE are defined
  # Gain = 1.0 and RDNOISE = 0.0 are the only feasible assumptions one can
  # make if nothing else is known.

  data = pyfits.open(frame, 'update')
  try:
    gain = data[0].header['GAIN']
  except KeyError:
    data[0].header.update('GAIN', 1.0)
    messageLog.put_warning('Could not find GAIN header: assumed GAIN = 1.0')
    data[0].header.add_history('Inserted missing GAIN header card, assumed GAIN = 1.0')
  try:
    rdnoise = data[0].header['RDNOISE']
  except KeyError:
    data[0].header.update('RDNOISE', 0.0)
    messageLog.put_warning('Could not find RDNOISE header: assumed RDNOISE = 0.0')
    data[0].header.add_history('Inserted missing RDNOISE header card, assumed RDNOISE = 0.0')
  data.close()


################################################################################

def remove_bg(frame, extension=0, poly_order_x=10, poly_order_y=8, nsteps=100, max_gradient=5, orderdef=None):

  medianfiltersize = 5
  minfiltersize  = 10

  if nsteps <= poly_order_y: poly_order_y = int(nsteps / 2)

  im = frame[extension].data
  (ysize, xsize) = im.shape

  bgim    = numpy.zeros((ysize, xsize),  dtype=numpy.float64)
  xfitarr = numpy.zeros((nsteps, xsize), dtype=numpy.float64)

  xx = numpy.arange(xsize, dtype=numpy.float64)
  yy = numpy.arange(ysize, dtype=numpy.float64)
  
  xx_mean = xx.mean()
  yy_mean = yy.mean()

  ystep = int(ysize / (nsteps - 1))

  yvals = (numpy.arange(nsteps) * ystep).round()

  ycount = 0


  
  for yind in yvals:

    #medianvec = numpy.zeros(xsize, dtype=numpy.float64)
    minvec  = numpy.zeros(xsize, dtype=numpy.float64)

    ymin_ind = numpy.max([yind-medianfiltersize, 0])
    ymax_ind = numpy.min([yind+medianfiltersize, ysize-1])

    naver = ymax_ind - ymin_ind + 1

    meanvec = numpy.average(im[ymin_ind:ymax_ind, :], axis=0)
    kernel = numpy.repeat(0.2, 5)
    d_meanvec = numpy.gradient(numpy.convolve(meanvec, kernel, mode='sample'))
    d_meanvec[0:10] = 0
    d_meanvec[xsize-11:xsize-1] = 0


    order_throughs = []
    order_throughs2 = []
    order_peaks = []
    all_peaks = []



    # First try to detect the peaks myself...
    for j in numpy.arange(xsize-1) + 1:
     if (d_meanvec[j-1] > max_gradient and d_meanvec[j] < max_gradient): all_peaks.append(j)


    # Now make use of existing orderdef if it exists
    if orderdef is not None:

      for order in orderdef:
        centery  = order['centerpix']
        centerx  = order['centerval']

        cheb_min = order['coefs'][2]
        cheb_max = order['coefs'][3]
        chebcoef = order['coefs'][4:]

	orderpos = _chebev(cheb_min, cheb_max, chebcoef, numpy.array([yind+1]))
        order_peaks.append(int(orderpos[0] + centerx) - 1)

      # But also append 'own' peaks below and above defined orders
      min_order_peak = numpy.min(order_peaks)
      max_order_peak = numpy.max(order_peaks)
      for peak in all_peaks:
        if peak < min_order_peak: order_peaks.append(peak)
        if peak > max_order_peak: order_peaks.append(peak)

    else:
      # Use whatever we have
      order_peaks = all_peaks

    order_peaks = numpy.array(order_peaks)
    npeaks = order_peaks.size
    order_throughs = (order_peaks[1:] + order_peaks[:npeaks-1]) / 2

    # Correct the max gradient for frames with very low pixel counts
    max_gradient = numpy.min((5, numpy.mean(meanvec.take(order_peaks)) / 4))


    # Now start searching for points that sample the background
    # (scattered light level)

    # First, find low points before first order, with
    # at least one point every 20 pixels
    for xpix in numpy.arange(0, order_peaks.min(), 20):
      shortvec = numpy.arange(20) + xpix
      shortvec = numpy.extract(shortvec < order_peaks.min(), shortvec)

      if (shortvec.size > 0):
        order_throughs2.append(shortvec[numpy.argmin(meanvec.take(shortvec))])

# The following is fancy, but it does not work since derivative is put to
# zero at last pixels and hence it cannot be used to find peaks or troughs.
#
#    (inds,) = numpy.where(short_dd > -max_gradient)
#    if inds.size == 0: inds = short_dd.size
#    min_ind = numpy.min(inds)
#    max_ind = short_dd.size
#
#    for j in numpy.arange(min_ind, max_ind - min_ind):
#      order_throughs2.append(shortvec[j])
#      nfound = nfound + 1
#    if nfound == 0:
#      order_throughs2.append(shortvec[numpy.argmin(meanvec.take(shortvec))])



    # Then loop over orders and find the troughs between the orders
    for xpix in order_throughs:
      shortvec = numpy.arange(11) - 5 + xpix
      shortvec = numpy.extract(shortvec > 1,	   shortvec)
      shortvec = numpy.extract(shortvec < xsize-1, shortvec)
      short_dd = d_meanvec.take(shortvec)

      # Use gradient to find region between peaks
      (inds,) = numpy.where(short_dd < -max_gradient)
      if inds.size == 0: inds = 0
      min_ind = numpy.max(inds)
      (inds,) = numpy.where(short_dd >  max_gradient)
      if inds.size == 0: inds = short_dd.size
      max_ind = numpy.min(inds)

      # Add the found region to the array of troughs
      nfound = 0
      for j in numpy.arange(min_ind, max_ind - min_ind):
	order_throughs2.append(shortvec[j])
	nfound = nfound + 1
      if nfound == 0:
        # But, if no through found, try to take minimum value of meanvec at position of trough
	try: order_throughs2.append(shortvec[numpy.argmin(meanvec.take(shortvec))])
	except: pass

        # Old solution - middle point between the order peaks
	# Does not work well because of non-monotonous order spacing
	# order_throughs2.append(xpix)

      # Safety catch: remove any incorrectly selected order peaks...
      try: order_throughs2.remove(shortvec[j])
      except: pass



    # Finally, find lowest point after last order
    shortvec = numpy.arange(50) + order_peaks.max()
    shortvec = numpy.extract(shortvec < xsize-1, shortvec)
    short_dd = d_meanvec.take(shortvec)

    if shortvec.size > 0:
      order_throughs2.append(shortvec[numpy.argmin(meanvec.take(shortvec))])


    # Finally, make sure each point is only selected once
    order_throughs = numpy.unique(order_throughs2)


    niter = 0
    # Perform fitting with sigma-clipping
    while 1:

      coefs = numpy.polyfit(order_throughs - xx_mean, meanvec.take(order_throughs),
                            poly_order_x)
      xfit_poly = numpy.poly1d(coefs)
      xfit = xfit_poly(order_throughs - xx_mean)

      sigma = numpy.abs(meanvec.take(order_throughs) - xfit) / numpy.std(meanvec.take(order_throughs) - xfit)
      rejected = numpy.extract(sigma > 3, order_throughs)
      order_throughs = numpy.extract(sigma < 3, order_throughs)
      
      niter = niter + 1
      if niter == 5 or rejected.size == 0: break
    
    xfit = xfit_poly(xx - xx_mean)
    xfitarr[ycount, :] = xfit

    ycount = ycount + 1


    # Debugging plot in x-direction (across the orders)
    if 0:
    #if ycount % 10 == 0:

      import biggles
      biggles.configure('screen', 'width',  900)
      biggles.configure('screen', 'height', 450)
      fr = biggles.FramedPlot()
      fr.yrange = [-50, meanvec.take(order_throughs).max()+50]
      fr.add(biggles.Curve(xx, meanvec, color='green'))
      fr.add(biggles.Curve(xx, xfit, color='blue'))
      fr.add(biggles.Points(order_throughs, meanvec.take(order_throughs)))
      fr.add(biggles.Points(order_peaks, meanvec.take(order_peaks), color='red'))
      fr.title = str(yind)
      fr.show()


  # Now fit in the y-direction
  for xind in numpy.arange(xsize):

    # Perform fitting with sigma-clipping
    niter = 0
    goodind = numpy.arange(nsteps)

    while 1:

      coefs = numpy.polyfit(yvals.take(goodind) - yy_mean, xfitarr[goodind, xind], poly_order_y)

      yfit_poly = numpy.poly1d(coefs)
      yfit = yfit_poly(yvals.take(goodind) - yy_mean)

      sigma = (xfitarr[goodind, xind] - yfit) / numpy.std(xfitarr[goodind, xind] - yfit)
      rejected = numpy.extract(sigma > 3, goodind)
      goodind  = numpy.extract(sigma < 3, goodind)

      niter = niter + 1
      if niter == 3 or rejected.size == 0 or goodind.size == 0: break

    if goodind.size == 0: 
      print "Error: no points left when y-fitting the background"
      coefs=numpy.polyfit(xfitarr[:, xind])
    bgim[:, xind] = yfit_poly(yy - yy_mean)


    # Debugging plot in y-direction (along the orders)
    if 0:
    #if xind % 250 == 0:

      import biggles
      biggles.configure('screen', 'width',  900)
      biggles.configure('screen', 'height', 450)
      fr = biggles.FramedPlot()
      fr.yrange = [bgim[:, xind].min()-50, bgim[:, xind].max()+50]
      fr.add(biggles.Curve(yy, bgim[:, xind], color='blue'))
      fr.add(biggles.Points(yvals, xfitarr[:, xind], color='red'))
      fr.add(biggles.Points(yvals.take(goodind), xfitarr[goodind, xind], color='green'))
      fr.title = str(xind)
      fr.show()


  frame[extension].data = frame[extension].data - bgim

################################################################################
