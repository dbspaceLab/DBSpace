# DBSpace
A model-centric Python library to study the dynamics and control of deep brain stimulation.

## Overview

## Structure

## Example

### Choosing files
spot_check.py is meant to be a quick option to view an LFP recording. It can be run (F5) without any arguments to open up a file chooser and select an LFP recording. Once chosen, the file chooser adds the file to a running list of files to visualize. The file chooser then opens again to enable further selection of files.

If the user is done with selecting the files of interest, then the user presses cancel in the file chooser which triggers the remainder of the script.

### Pre-selected files
A file list can be chosen *a priori* for spot_check.py to run through. This is currently implemented as a list being placed into the script, but will be implemented as an argument that can be passed to the script.

## Outputs and Visualizations
The key visualizations done for each of the selected files is (1) a spectrogram of the left and right channels and (2) a PSD of a particular interval in the recording for the left and right channels.

## Historical Note
This library was developed slowly over the course of my PhD work under Dr. Helen Mayberg and Dr. Robert Butera.
The PhD was focused on analysing timeseries from intracranial and extracranial recordings, particularly problematic ones.
The importance of establishing a model _a priori_ in analysing neural recordings became clear - especially because we had high-dimensional models that were being used ever-day in the management of the patient.

The goal of this library was to recented DBS around the models we were using to treat patients, and enable a way to check in with data to shape those models.

