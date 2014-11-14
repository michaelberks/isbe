#!/bin/bash

#
# == Set SGE options:
#

#$ -S /bin/bash
#$ -cwd
  #
  # -- ensure BASH is used, no matter what the queue config;
  # -- run the job in the current-working-directory;
  # -- submit the job to the "serial.q" queue;
  #

#
# == The job itself :
#

echo " "
echo "Job ran on and started at: "
/bin/hostname
/bin/date
echo " "

export PATH=$PATH:cxx/vxl/bin/isbe_apm/tools/
  #
  # ...ensure VXL RF executable is on our path...
  #

mrfr_build_vec_trees

echo " "
echo "Job finished at: "
/bin/date