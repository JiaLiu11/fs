#!/usr/bin/env python 
# 
#  file: lm_cmdline.py
#  
#  Author:  Jia liu    <liu.2053@osu.edu>
#
#  Revision history:
#	  May-19-2013  Original version
#  
#*************************************************************************
from subprocess import call    # make the "call" function available
import numpy

# Define the inputs at the top, so they are easy to find and change
value_list = numpy.linspace(1,3,num=3)

# Ok, let's do it . . .
print "\nTry to run the lm.e for several times:"

for nevents in value_list:    # don't forget the colon!
  my_command = './lm.e ' + str(nevents)# convert nevents to a string
  retcode = call(my_command, shell=True)    # pass "my_command" to be executed
  print "This run completes!\n\n"

print "The last return code was", retcode

#*************************************************************************
