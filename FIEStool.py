#!/usr/bin/env python

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
# FIEStool.py
#
# Last update 02/03/2009
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.		       #
    #         The body of this program is at the end of this file      #
    #								       #
    ####################################################################


"""
   Main routines for the automatic data reduction interface.
"""

_version     = "1.3.2"
_date        = "27/07/2009"
_author      = "Eric Stempels"


_minversion_numpy   = "1.0.1"
_minversion_pyfits  = "1.1"
_minversion_pyraf   = "1.4"
_minversion_biggles = "1.6.4"


################################################################################

# Import necessary modules

import checkImport
import config
import configFrame
import calibFrame
import messageLog
import pipeFrame
import plotFrame
import autoQueue
import autoLoader
import Queue
import popUp
import os
import tkMessageBox
import getopt

from Tkinter import *
from Dialog import Dialog

################################################################################

class AutoScrollbar(Scrollbar):
    # a scrollbar that hides itself if it's not needed.  only
    # works if you use the grid geometry manager.
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            # grid_remove is currently missing from Tkinter!
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)
    def pack(self, **kw):
        raise TclError, "cannot use pack with this widget"
    def place(self, **kw):
        raise TclError, "cannot use place with this widget"

################################################################################


class mainWindow(Frame):

  "The main control window that generates and controls all other windows"

  def __init__(self, parent, hasscrollbar=False):

    ###############################################
    #     First prepare and create this frame     #
    ###############################################
    self.hasscrollbar = hasscrollbar

    # GUIupdater' contains a Queue.Queue object into
    # which calls for updating the GUI may be placed (useful when running
    # as a thread). These variables are global in this class
    self.GUIupdater = Queue.Queue()

    # Was an additional scrollbar requested?
    if self.hasscrollbar:

      # Add a scrollbar
      vscrollbar = AutoScrollbar(parent)
      vscrollbar.grid(row=0, column=1, sticky=N+S)

      # Create a Canvas object that can contain the whole frame
      self.canvas = Canvas(parent, yscrollcommand=vscrollbar.set)
      self.canvas.columnconfigure(0, weight=1)
      self.canvas.rowconfigure(0, weight=1)

      # Attach the method to position the canvas to the scrollbar
      vscrollbar.config(command=self.canvas.yview)
      
      # Inialize as if this were an instance of the Frame class
      Frame.__init__(self, self.canvas)

    else:

      # Inialize as if this were an instance of the Frame class
      Frame.__init__(self, parent)

    # And make 'self' frame flexible in width and heigth
    self.columnconfigure(0, weight=1)
    self.rowconfigure(0, weight=1)

    # Build this frame
    self.makeFrame()

    ###############################################
    # Now prepare and create all the child frames #
    ###############################################

    # These are 'all' the editable options that should appear in the main configuration frame
    alloptions =      ['inpath', 'outpath', 'filename_filter', 'biaslist', 'flatlist',
                       'orderdef', 'wavedef', 'interlacedorderdef', 'interlacedwavedef',
		       'refpath', 'pixelmask', 'masterbias', 'masterflat',
		       'masternormflat', 'blazeshape', 'masterorderdef', 'waveref', 'masterwaveref',
		       'masterinterlacedorderdef', 'interlacedwaveref',
		       'fitsheaders', 'mef_dataext', 'frameorientation',
		       'config_dir', 'iraf_taskconfig_dir', 'x_start', 'x_end', 'y_start', 'y_end',
		       ]

    # Create the main configuration frame
    self.mainConfig = configFrame.window(self, optionlist=alloptions, title='FIEStool - Configure all settings')

    # Create the frame for calculating calibration (reference) data
    self.CF = calibFrame.window(self, self.GUIupdater, title='FIEStool - Calculate reference data')

    # Create the autoLoader frame
    self.AL = autoLoader.window(self)

    # Attach a menubar to this frame. This can only be done after the other frames
    # have been created, because the menubar makes calls to these frames
    self.MB = MenuBar(self)
    parent.config(menu = self.MB)

    # Start automatic updating of this frame
    self.updateGUI()

    if self.hasscrollbar:
      # Store it in a canvas object (to allow scrolling)
      self.canvas.create_window(0,0, anchor=NW, window=self)


  def makeFrame(self):

    "Construct the frame from scratch"

    # Put all data in a root frame that is embedded in the 'self' Frame object.
    # In this way, it is possible to destroy and rebuid the frame object,
    # without destroying the 'self' frame object.
    self.root = Frame(self)

    # The root frame should be flexible in width
    self.root.columnconfigure(0, weight=1)

    # Row counter
    i = 0


    # Add title
    Label(self.root, text='FIES automated data reduction interface').grid(row=i, sticky=E+W)

    # Next row
    i = i + 1


    # Insert an embedded configuration window so one can easily see and change th
    # location of the input and output frames
    self.CF = configFrame.window(self.root, optionlist=['inpath', 'outpath'], embed=1, width=200)

    # Display this embedded frame
    self.CF.grid(row=i, sticky=E+W)

    # Next row
    i = i + 1


    # Insert the autoQueue frame, controlling the queue of files to be processed
    self.AQ = autoQueue.window(self.root)

    # Display this embedded frame
    self.AQ.grid(row=i, sticky=E+W)

    # Next row
    i = i + 1


    # Create and display the button to start processing
    self.startbutton = Button(self.root, text='START processing',
                             command=self.startProcess)
    self.startbutton.grid(sticky=E+W)

    # Save the default fg and bg colors
    self.startbutton.defaultbgcolor     = self.startbutton.cget('bg')
    self.startbutton.defaultactivecolor = self.startbutton.cget('activebackground')
    
    # Add a popup to this button
    startPopup = popUp.popUp(self.root, self.startbutton, title='Start/Stop processing',
                          text="""Start (or stop) the processing of files in the queue. Only
			          one file at a time will be processed. Stopping the queue
				  processing will not interrupt any active process.
			       """)

    # Next row
    i = i + 1


    # Create and display the window containing the pipeline tasks
    self.PL = pipeFrame.window(self.root, self.GUIupdater)
    self.PL.grid(row=i, sticky=E+W)

    # Next row
    i = i + 1


    # Create and display the window controlling all plotting
    self.PF = plotFrame.window(self.root)
    self.PF.grid(row=i, sticky=E+W)

    # Next row
    i = i + 1


    # Create and display the logging window
    self.ML = messageLog.window(self.root)
    self.ML.grid(row=i, sticky=N+E+S+W)
    self.root.rowconfigure(i, weight=1)
    # Next row
    i = i + 1


    # Now display the entire frame
    self.root.grid(sticky=N+E+S+W)


  def updateGUI(self):

    """
       Checks if any GUI events have been queued in the GUIupdater Queue
       and executes them
    """

    while not self.GUIupdater.empty():
      command, args, kwds = self.GUIupdater.get()
      command(*args, **kwds)

    # Reschedule itself
    self.after(100, self.updateGUI)


  def startProcess(self):

    """
       Prepares this frame for automatic processing. Change the color
       and text of the start button and start the automatic queue checking.
    """

    # Make sure there is no other processing active
    try:
      self.after_cancel(self.pipeCheckID)
    except:
      pass

    messageLog.put('Automatic processing started', 5)

    # Change the color and text of the start/stop button
    self.startbutton.config(background = '#bb6666',
                            activebackground = '#dd6666',
			    text = 'STOP processing',
                            command = self.stopProcess)

    # Disable the pull-down menu. You don't want the user to change reduction mode
    # or configuration parameters while processing!
    self.MB.disable()

    # Start the periodic processing    
    self.Process()


  def Process(self):

    """
       Check for non-processed frames in the queue. As soon as a new frame
       appears (by user selection, or by automatic queuing) the reduction
       pipeline will be started.
    """

    messageLog.put('Checking queue for files to process', 9)

    # Update and check the number of frames in the queue
    self.AQ.updateQueueSize()
    nfiles = self.AQ.queuesize.get()

    messageLog.put('Found %i files in queue' % nfiles, 9)

    # Are there any pending frames?
    if (nfiles > 0): 

      # Yes! Pop a name from the queue
      nextframe = self.AQ.newfiles[0]

      # Start the pipelined processing (which may be using threads!)
      try: 
        self.PL.process_frame(nextframe)
        self.AQ.newfiles.remove(nextframe)
	# Run a process to check is the pipeline (thread) has produced any output
	self.checkPipeOut()
      except pipeFrame.Busy: 
        # Don't start reducing a frame if another frame is being processed.
        messageLog.put('Pipeline is busy - waiting...', 9)

      # Update the queue size
      self.AQ.updateQueueSize()

    else:

      # No... do nothing...
      messageLog.put('No files to process - waiting...', 9)

    # Reschedule this action
    self.pipeCheckID = self.after(1000, self.Process)


  def stopProcess(self):

    """
       Stop the automatic processing of frames and return the start button to
       its original text and color.
    """

    # Cancel the automatic checking for new frames (failsafe)
    try:
      self.after_cancel(self.pipeCheckID)
    except:
      pass

    messageLog.put('Automatic processing stopped', 5)

    # Revert the button to its original color and text
    self.startbutton.config(bg = self.startbutton.defaultbgcolor, 
                            activebackground = self.startbutton.defaultactivecolor,
			    text='START processing',
                            command = self.startProcess)

    # Re-enable the menu
    self.MB.enable()


  def checkPipeOut(self):

    """
       Do periodic checking for output from the pipeline. Because the pipeline
       may run in a background thread, the program must check if the thread has
       finished (i.e. it produced output) before starting autoplotting of the
       result and adjust the file queues.
    """

    try:
      # Retrieve output from the queue, not waiting for an answer if the
      # queue was empty, in which case the program raises the exception Queue.Empty
      (inframe, outframe) = self.PL.outQueue.get_nowait()
      messageLog.put('Got "%s" from output queue' % outframe, 9)

      # Add the file to the list of reduced files
      if inframe not in self.AQ.reducedfiles:
	self.AQ.reducedfiles.append(inframe)
      if inframe in self.AQ.newfiles:
	self.AQ.newfiles.remove(inframe)

      # Update all queue sizes
      self.AQ.updateQueueSize()

      # If output was produced and it really exists, set the default plotting name
      if outframe and os.path.exists(outframe):
	self.PF.lastplotname = outframe

    except Queue.Empty:

      # Pipeline has not finished yet, reschedule!
      self.after(1000, self.checkPipeOut)


  def about(self):

    "Display a box with information about this program"
  
    tkMessageBox.showinfo(message="FIEStool - version %s\nDate: %s\nAuthor: %s\nCopyright (c) 2005 NOTSA\nAll rights reserved" % (_version, _date, _author))


  def help(self):
  
    "Display a box with very short help information"

    tkMessageBox.showinfo(message="Right-click on any item for more info. Or read the manual...")


  def really_quit(self):
  
    "Double-check with the user if this frame (and thus the program) really should terminate."
    answer = Dialog(self, title='', bitmap='',
           strings=('Yes', 'No'), text='Really quit?', default=1)
    if answer.num == 0: self.quit()



################################################################################

class MenuBar(Menu):

  """
     Provides a pull-down menu to control the interface. It also contains all the
     hooks to the different windows and associated actions
  """

  def __init__(self, parent, relief=GROOVE, bd=2):

    # Create a counter for the number of menu items (useful for enable and disable)
    self.nmenu = 0

    # Create a menu object
    Menu.__init__(self, parent, relief=GROOVE, bd=2)


    # Create first pull-down menu ('file')
    filemenu = Menu(self, tearoff=0)
    # Add the 'quit' function
    filemenu.add_command(label='Quit',				command=parent.really_quit)


    # Create another pull-down menu ('edit') and add items/functions
    settmenu = Menu(self, tearoff=0)
    settmenu.add_command(label='Edit all settings',		command=parent.mainConfig.show)
    settmenu.add_command(label='Load settings',			command=parent.mainConfig.loadConfig)
    settmenu.add_command(label='Save settings',			command=parent.mainConfig.saveConfig)
    settmenu.add_separator()
    settmenu.add_command(label='Restore default settings',	command=lambda: config.load(config.options['default_savefile'].value,
                                                                                            ignore_options=['currentmode']))
    settmenu.add_command(label='Save current settings as default',
    								command=lambda: config.save(config.options['default_savefile'].value,
								                            ignore_options=['currentmode']))

    automenu = Menu(self, tearoff=0)
    automenu.add_command(label='Configure autoloading',		command=parent.AL.show)
    automenu.add_command(label='Load autoloader settings',	command=parent.AL.loadConfig)
    automenu.add_command(label='Save autoloader settings',	command=parent.AL.saveConfig)
    automenu.add_separator()
    automenu.add_command(label='Restore default autoloader settings',	command=lambda: config.load_raw(config.options['default_autoload_savefile'].value))
    automenu.add_command(label='Save current autoloader settings as default',
    								command=lambda: config.save_raw(config.options['default_autoload_savefile'].value,
								                                optionlist=['autoloadrules']))

    # Create another pull-down menu ('calib') and add items/functions
    calibmenu = Menu(self, tearoff=0)
    # This is a submenu in the 'calib' menu
    calibmenu2 = Menu(calibmenu, tearoff=0)
    calibmenu2.add_command(label="Yes, I know what I'm doing", command=parent.CF.show)
    # This is the actual calib menu entry
    calibmenu.add_cascade(label='Calculate calibration frames', menu=calibmenu2)

    # Create another pull-down menu ('log') and add items/functions
    logmenu = Menu(self, tearoff=0)
    # This is a submenu in the 'log' menu
    loglevelmenu = Menu(logmenu, tearoff=0)
    # Create a range of numbers from 1 to 9 to select log-level detail
    for i in range(9): loglevelmenu.add_radiobutton(label=i+1, value=i+1,
                                             variable=parent.ML.loglevel )
    # Add other items/functions
    logmenu.add_command(label='Clear log',	command=parent.ML.clear)
    logmenu.add_command(label='Dump config to logfile',	command=config.dump)
    logmenu.add_command(label='Save to file',	command=parent.ML.save)
    logmenu.add_separator()
    logmenu.add_cascade(label='Set log level',	menu=loglevelmenu)

    # Create another pull-down menu ('mode') and add items/functions
    modemenu = Menu(self, tearoff=0)
    for mode in config.options['availablemodes'].value:
      modemenu.add_radiobutton(label=mode, 	value=mode, variable=parent.PL.mode,
    						command=parent.PL.reset)

    # Create another pull-down menu ('help') and add items/functions
    helpmenu = Menu(self, tearoff=0)
    helpmenu.add_command(label='Help',		command=parent.help)
    helpmenu.add_command(label='About',		command=parent.about)


    # Add, one by one the pull-down menu objects to the menubar, keeping track of
    # the number of menus added
    self.add_cascade(label='File',		menu=filemenu)
    self.nmenu = self.nmenu + 1

    self.add_cascade(label='Settings',		menu=settmenu)
    self.nmenu = self.nmenu + 1

    self.add_cascade(label='AutoLoader',	menu=automenu)
    self.nmenu = self.nmenu + 1

    # The if-statement is now deprecated
    #if parent.enablecalib:
    self.add_cascade(label='Calibs',		menu=calibmenu)
    self.nmenu = self.nmenu + 1

    self.add_cascade(label='Log',		menu=logmenu)
    self.nmenu = self.nmenu + 1

    self.add_cascade(label='Mode',		menu=modemenu)
    self.nmenu = self.nmenu + 1

    self.add_cascade(label='Help',		menu=helpmenu)
    self.nmenu = self.nmenu + 1



  def disable(self):

    # Disable all items in the menubar
    i = 1
    while i <= self.nmenu:
      try:
        self.entryconfig(i, state=DISABLED)
	i = i + 1
      except:
        break


  def enable(self):
  
    # Enable all items in the menubar
    i = 1
    while i <= self.nmenu:
      try:
        self.entryconfig(i, state=NORMAL)
	i = i + 1
      except:
        break


################################################################################


def checkSetup():

  try:
    checkImport.checkImport('numpy',   _minversion_numpy)
    checkImport.checkImport('pyfits',  _minversion_pyfits)
    checkImport.checkImport('pyraf',   _minversion_pyraf)
    checkImport.checkImport('biggles', _minversion_biggles)
  except ImportError:
    sys.exit(1)

  if not os.path.exists('login.cl'):
    print "No 'login.cl' found in the FIEStool directory!"
    print "You need to run 'mkiraf' first to create this file."
    print "(choice of terminal does not matter)"
    sys.exit(1)


# This 'run()' module method in principle allows other programs to import this
# file as a module and start the reduction interface. Maybe never needed in
# this way.

def run():

  """
     Main program entry, initializing and starting the main window and managing
     its 'mainloop'
  """

  # Get and check command-line options
  #
  # NB: The -c (enable-calib) option is now deprecated
  #
  try:
    opts, args = getopt.getopt(sys.argv[1:], "cs", ["calib", "scrollbar"])
  except getopt.GetoptError:
    # print help information and exit:
    print "Unknown command line parameter. Allowed optional parameters are : "
    #print "-c [file]            Load calibration file at startup"
    print "-s / --scrollbar	Add vertical scrollbar is screen does not fit"
    sys.exit(2)

  #enablecalib  = False
  hasscrollbar = False

  for option, parameter in opts:
    if option in ("-c", "--calib"):
      print "The -c option is deprecated. You can now always access the calibs submenu."
    #  enablecalib = True
    if option in ("-s", "--scrollbar"):
      hasscrollbar = True

  # Display a welcome message in the message log
  messageLog.put('='*65)
  messageLog.put("Welcome to FIEStool - the automatic data processing tool")
  messageLog.put("for the FIES spectrograph at the Nordic Optical Telescope")
  messageLog.put("Copyright (c) 2005 NOTSA; All Rights Reserved")
  messageLog.put("This is version %s (%s), written by %s" % (_version, _date, _author))
  messageLog.put('='*65)
  messageLog.put("Right-click on any item to obtain instantaneous help")
  messageLog.put('='*65)

  # And tell the user what warnings and errors look like
  messageLog.put_warning("This is a sample warning")
  messageLog.put_error("This is a sample error")
  messageLog.put('='*65)

  # Now, check if a default config file exists; if yes, load it.
  try:
    config.load(config.options['default_savefile'].value)
  except IOError:
    messageLog.put_warning('No default configuration file found!')
    pass

  # Now, try to load the autoLoader config file.
  try:
    config.load_raw(config.options['default_autoload_savefile'].value)
  except IOError:
    messageLog.put_warning('No default autoloader configuration file found!')
    pass

  # Create the absolute top-level root window. Destroying this kills 'everything'!
  root = Tk()

  # Make the root window variable in width and height
  root.rowconfigure(0, weight=1)
  root.columnconfigure(0, weight=1)

  # Give it a title
  root.title('FIEStool - Data reduction for the FIES spectrograph')

  # Initialize the main window, with 'root' as its parent
  main = mainWindow(root, hasscrollbar)

  if hasscrollbar:

    # This puts the focus right, and will give the canvas child real sizes
    root.deiconify()
    root.focus_set()
    root.grab_set()
    root.wait_visibility()

    # Show the canvas
    main.canvas.grid(row=0, column=0, sticky=N+S+E+W)

    # Determine width and heigth of the contained frame
    canvaschild  = main.canvas.children.values()[0]

    neededwidth  = canvaschild.winfo_width()
    neededheight = canvaschild.winfo_height()

    maxheight = main.canvas.winfo_screenheight()
    maxwidth  = main.canvas.winfo_screenwidth()

    neededheight = min(maxheight, neededheight)
    neededwidth  = min(maxwidth,  neededwidth)

    # Set the canvas size and enable scrolling
    main.canvas.config(width=neededwidth, height=neededheight)
    main.canvas.config(scrollregion=main.canvas.bbox("all"))

  else:

    # Display the main window
    main.grid(sticky=N+E+S+W)

  # Catch the 'delete window' signal that may be triggered by the window
  # manager if the user clicks on its 'close' button
  root.protocol("WM_DELETE_WINDOW", main.really_quit)


  # Start the mainloop. The program executing pauses here, only to execute events
  # generated by the main window (which in essence is everything)

  ###############
  main.mainloop()
  ###############


  # Apparently an event triggered an exit out of the mainloop (see main.reallyquit)
  messageLog.put("Destroying widgets")

  # Destroy the root window and all its children
  root.destroy()


  # This message will never arrive, though...
  messageLog.put("Program finished - EXIT FIEStool")

  # Return to caller (main)
  return


# Here's where the program starts
if __name__ == "__main__": 

  # Check if the minimum required external packages are installed
  checkSetup()

  # Run it!
  run()

  # That's it - BYE!
