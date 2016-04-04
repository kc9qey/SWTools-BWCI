#!/bin/ksh

############################## standard interface to /sw tools
# Input:
#   Environment variables
#     SW_BLDDIR    current directory (PWD) minus /autofs/na1_ stuff
#     SW_ENVFILE   file to be sourced which has alternate prog environment
#                     only to be used in special circumstances
#     SW_WORKDIR   work dir that local script can use
# Output:
#   Return code of 0=success or 1=failure   or 2=job submitted
#
# Notes:
#   If this script is called from swtest, then swtest requires 
#   SW_WORKDIR to be set.  Then swtest adds a unique path to what 
#   user gave swtest (action+timestamp+build) and provides this
#   script with a uniquely valued SW_WORKDIR.  swtest will
#   automatically remove this unique workspace when retest is done.
##################################################################

# exit 3 is a signal to the sw infrastructure that this template has not 
# been updated; please delete it when ready

if [ -z $SW_BLDDIR ]; then
  echo "Error: SW_BLDDIR not set!"
  exit 1
else
  cd $SW_BLDDIR
fi

if [ -z $SW_ENVFILE ]; then
  ### Set Environment (do not remove this line only change what is in between)
  . ${MODULESHOME}/init/ksh
  . ${SW_BLDDIR}/remodule
  ### End Environment (do not remove this line only change what is in between)
else
  . $SW_ENVFILE
fi

############################## app specific section
#  

#clear out status file since re-testing
rm -f status test.log

cd ${SW_WORKDIR}
tar -xf /sw/testcases/${PACKAGE}/apoa1.tar.gz 
mv apoa1/* .
cat > namd.pbs << EOF
#!/bin/bash
#PBS -N namd.apoa1
#PBS -j oe
#PBS -l size=32,walltime=01:00:00
#PBS -A UT-SUPPORT
#PBS -W group_list=install
#PBS -W umask=002

set -o verbose
cd \$PBS_O_WORKDIR

NAMD_BIN=${SW_BLDDIR}/bin/namd2
aprun -n 4 \$NAMD_BIN apoa1.namd > apoa1.namd.out
total=\`grep "ENERGY:     500" apoa1.namd.out | cut -d" " -f73| cut -d. -f1\`
echo \$total

if [ "\$total" == "-222495" ]; then
  echo "verified" > ${SW_BLDDIR}/status 
else
 echo "unverified" > ${SW_BLDDIR}/status
fi  

JOBID=\`echo \$PBS_JOBID | cut -d "." -f1 \`
chmod 775 ${SW_BLDDIR}/status
rm ${SW_BLDDIR}/.running
cp namd.apoa1.o\$JOBID ${SW_BLDDIR}/test.log
cat apoa1.namd.out >> ${SW_BLDDIR}/test.log
chmod 664 ${SW_BLDDIR}/test.log
EOF

#submit job and touch .running file - marker to infrastructure that
qsub namd.pbs 
touch ${SW_BLDDIR}/.running

# qsub returns 0 on successful job launch, so if failure return 1
if [ $? -ne 0 ]; then
  rm -f ${SW_BLDDIR}/.running
  exit 1
else
  exit 2
fi

cd ../

############################### if this far, return 0
exit 0
