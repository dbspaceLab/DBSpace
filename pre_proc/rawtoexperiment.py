"""
The objective of this script will be to create an experiment folder with duplicate raw files renamed
to actual experiments.


Created on Mon June 4th 11:42AM 2018

@author: Dolu Obatusin

"""

# import required modules / libraries
import os
import numpy as np
import json
import re
from datetime import datetime
import shutil
import json
import sys, getopt

# Note the location of the raw files output from Sensing tablet

def experimentdir(inputfile, outputfile):
    # Create temporary folder as experiment folder
    if not(os.path.exists(outputfile)):
        os.makedirs(outputfile)

# # Function to extract Patient ID from path
#
#
# def extract_patientid(basedir):
#     paths = []
#     pattern = re.compile("DBS")
#     i = 0
#     for path, directories, files in os.walk(basedir, topdown=True):
#         paths.append(path)
#         if files:
#             for match in re.finditer(pattern, path):
#                 break
#         #             print(match)
#         #             print('Found on line %s' % (i+1))
#         i += 1
#     patient_id = paths[1][match.span()[0]:]
#     return patient_id
#
# # Function to get timestamps from log file
#
#
# def ins_enable_timestamps(log_file):
#     pattern = re.compile("Enable Sensing")
#     for i, line in enumerate(open(log_file)):
#         for match in re.finditer(pattern, line):
#             t_sensing = line[15:]
#             datetime_object1 = datetime.strptime(t_sensing.strip(), '%m/%d/%Y %I:%M:%S %p')
#     #             print("Enable Sensing:",datetime_object2)
#
#     pattern = re.compile("INS TimeStamp")
#     for i, line in enumerate(open(log_file)):
#         for match in re.finditer(pattern, line):
#             t_ins = line[14:]
#             datetime_object2 = datetime.strptime(t_ins.strip(), '%m/%d/%Y %I:%M:%S %p')
#     #             print("INS TimeStamp:",datetime_object1)
#     #         print('Found on line %s' % (i+1))
#
#     # Convert Enable Sensing to datetime
#
#     return datetime_object1, datetime_object2
#
# # Function to create log dictionary
#
#
# def create_logdict(root):
#     log_files = []
#     logfile_dict = {}
#     for path, directories, files in os.walk(root, topdown=True):
#         if files:
#             for f in files:
#                 if f.endswith('LOG.txt'):
#                     log_files.append(f)
#                     logfile_dict[path] = log_files
#             log_files = []
#     return logfile_dict
#
#
# # Create a data dictionary of path and log files - only include paths that have log files
#
#
# def createpathfiles():
#     logfile_dict = create_logdict(basedir)
#
#     pathfiles_dict = {}
#     path_files = []
#
#     # Operate according to each path
#     for key in logfile_dict:
#         #     print(key)
#
#         # create a folder in path based on log file
#         for k in logfile_dict[key]:
#             #         print(k)
#             log_file = key + '/'+ k
#             if k.endswith('LOG.txt'):
#                 #             pass
#                 new_folder = log_file.strip('_LOG.txt')
#                 if not(os.path.exists(new_folder)):
#                     os.makedirs(new_folder)
#             # Get INS & Sensing start time from log file
#             #         print(log_file)
#             datetime_object1, datetime_object2 = ins_enable_timestamps(log_file)
#             #         print(datetime_object1)
#             #         print(datetime_object2)
#             # Loop through path to get all files
#             #     print( ' ')
#             for path, directories, files in os.walk(key, topdown=True):
#                 if files:
#                     #             print(files)
#                     for f in files:
#                         if f.endswith('.xml') or f.endswith('.txt'):
#                             #                 print(f)
#                             datetime_object3 = datetime.strptime(f[6:27], '_%Y_%m_%d_%H_%M_%S_')
#                         #                     print(f,'--',datetime_object3)
#
#                         else:
#                             continue
#
#                         # compare each file's timestamp to INS & Sensing starttime and add to dictionary of {new folder: files}
#                         if datetime_object1 <= datetime_object3 <= datetime_object2:
#                             path_files.append( key + '/'+ f)
#                             pathfiles_dict[new_folder] = path_files
#                     path_files = []
#     return pathfiles_dict
#
#
#
# # Finally move files into their respective session folders
#
#
# def movefiles(pathfiles_dict):
#     for key in pathfiles_dict:
#         print(key)
#         print('-'*40)
#         for i in pathfiles_dict[key]:
#             #         print(i)
#             if os.path.exists(i):
#                 shutil.move(i, key)
#             else:
#                 continue
#
# # Get the file name for the new file to write
#
#
# def getfilename(pathfiles_dict):
#     patient_id = extract_patientid(basedir)
#     dst = basedir + '/' + patient_id + '_' + 'pathfiles.txt'
#     # Save map to text file
#     # If the file name exists, write a JSON string into the file.
#     if dst:
#         # Writing JSON data
#         with open(dst, 'w') as f:
#             json.dump(pathfiles_dict, f) # use `json.loads` to do the reverse


"""Collect command-line options in a dictionary"""


# def getopts(argv):
#     opts = {}  # Empty dictionary to store key-value pairs.
#     while argv:  # While there are arguments left to parse...
#         if argv[0][0] == '-':  # Found a "-name value" pair.
#             opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
#         argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
#     return opts
#
#
# if __name__ == '__main__':
#     from sys import argv
#     myargs = getopts(argv)
#     if '-i' in myargs:  # Example usage.
#         print(myargs['-i'])
#     print(myargs)


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print ('rawtoexperiment.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('rawtoexperiment.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print ('Input file is "', inputfile)
    print ('Output file is "', outputfile)




if __name__ == "__main__":
    main(sys.argv[1:])

