"""mxi ui: an interface to mxi, for processing data and for teaching.

The workflow's steps as forms that build the real mxi command, run it and show
what it wrote; the images, with spots, predictions and untrusted regions; the
results, as plots. tkinter for the windows and matplotlib for the plots, the
images and results read with h5py -- nothing scraped from logs.

Import-light on purpose: the window's modules import tkinter, and only they
do, so that the rest -- steps, runner, the directory's state, the dependency
check -- can be used and tested where there is no display or no tkinter.
"""
