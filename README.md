APEX SNAKE files directory
==========================

Overview
---------

This git repo is for the files needed by SNAKE for APEX experiment

Directory structure:
*APEX_doc : files that document the APEX septum setup.
*infiles  : input files used by SNAKE
*trajfiles: trajectory files used by SNAKE
*mapfiles : Magnet map files used by SNAKE
*macros   : Root C code/scripts
*matrix-files : fitted matrix elements

Each directory has a README.md to document the directory.

Setup
------

In addition to the files in the git repo, one should 
setup the following:

*Need to create the file sn-rand.dat which will store the initial seed for the
random number generator. The SNAKE code will overwrite this file with new seed at end of running.
*Symbolic link __snake__ to the location of the snake excutable.
*Symbolic link __DIR.DAT__ to the SNAKE input file.
*Symbolic link __TRAJ.DAT__ to the SNAKE trajectory file.
*Create directory __rootfiles__ to store root files created by the root macros.
*Create directory __plots__  to store pdf files created by the root macros.