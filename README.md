# DBSpace
A library to bridge clinical electrophysiology with dynamical/control systems modeling in deep brain stimulation.


## Overview
Deep brain stimulation (DBS) is growing rapidly as a therapy for neurological and psychiatric disorders, but we aren't sure how it works and what it does to the brain.
My PhD work was in 'reverse-engineering' the effects and efficacy of SCCwm-DBS for depression \cite{} and this library was a product of that work.

This library is meant to focus on the link between the messay world of clinical electrophysiology and the goals of dynamical neuroscience.
As such, it prioritizes an 'n=1' inference approach where our goal is to see how well individual patients fit the clinical models being used to treat them.
While the results of this work may be generalizable to broad neuropsychiatric care, this is not the primary goal.
Indeed, this library and its toolset are meant to embrace false-positive errors in the effort to generate models for further focused study in parallel with clinical care.

## Running Spot Check

### Choosing files
spot_check.py is meant to be a quick option to view an LFP recording. It can be run (F5) without any arguments to open up a file chooser and select an LFP recording. Once chosen, the file chooser adds the file to a running list of files to visualize. The file chooser then opens again to enable further selection of files.

If the user is done with selecting the files of interest, then the user presses cancel in the file chooser which triggers the remainder of the script.

### Pre-selected files
A file list can be chosen *a priori* for spot_check.py to run through. This is currently implemented as a list being placed into the script, but will be implemented as an argument that can be passed to the script.

## Outputs and Visualizations
The key visualizations done for each of the selected files is (1) a spectrogram of the left and right channels and (2) a PSD of a particular interval in the recording for the left and right channels.

## History

This library was developed slowly over the course of my PhD work under Dr. Helen Mayberg and Dr. Robert Butera.
