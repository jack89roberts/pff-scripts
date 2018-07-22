#!/bin/bash
#
# This script creates in the current directory soft links to necessary files of the model.
# In this way one uses the official model and doesn't need to copy it and keep it uptodate.
# NOTE: of course you will try to edit the official model! pay attention at what you do.
#
# Jan 2014
# davide.gamba@cern.ch
#

#Jack.Roberts@cern.ch - check username to determine whether control room environment or my laptop
scriptUser=$(whoami)

# update this variable if you move the model!
if [ "$scriptUser" = "ctf3op" ]
then
	echo "User is ctf3op, assuming control room environment"
	MODEL_PATH='/ctf/data/ctfmod/elements'
	MADX_EXECUTABLE='/ctf/data/ctfmod/madx_pro64'
elif [ "$scriptUser" = "jack" ]
then
	echo "User is jack, assuming environment on jack's laptop"	
	MODEL_PATH='/home/jack/Documents/ctf3model/elements'
	MADX_EXECUTABLE='/usr/local/bin/madx'
else
	echo "User name not jack or ctf3op. Don't know where model is."
	exit 1
fi

# Usage and check that we have an argument given.
if [ -z $1 ]; then
  echo " ";
  echo "   This script creates soft links to the seqx files of CTF3 model.";
  echo "   It also copies some commonly used files, e.g. the file with nominal optics currents.";
  echo "   Currently the madx model of CTF3 is stored in $MODEL_PATH."
  echo " ";
  echo "      Usage: $0 [directory where to create links]";
  echo " ";
  exit 1
fi



# lines sequences
ln -s $MODEL_PATH/common.seqx $1/
ln -s $MODEL_PATH/linac.seqx $1/
ln -s $MODEL_PATH/ct.seqx $1/
ln -s $MODEL_PATH/dl.seqx $1/
ln -s $MODEL_PATH/dl.wiggler.matching.madx $1/
ln -s $MODEL_PATH/tl1.seqx $1/
ln -s $MODEL_PATH/cr.seqx $1/
ln -s $MODEL_PATH/cr.wiggler.matching.madx $1/
ln -s $MODEL_PATH/cr.SInjection.seqx $1/
ln -s $MODEL_PATH/tl2.seqx $1/
ln -s $MODEL_PATH/tbts.seqx $1/
ln -s $MODEL_PATH/tbl.seqx $1/

ln -s $MODEL_PATH/califes.seqx $1/


# this will be anyway modified all the time, so better copy it.
cp $MODEL_PATH/ctf3.seqx $1/
# for the time being I just do a symlink...
#ln -s $MODEL_PATH/ctf3.seqx $1/

# same for nominal currents and a general starting script
cp $MODEL_PATH/nominalcurrents.txt $1/
cp $MODEL_PATH/generalScript.madx $1/

# create a link to current madx executable
ln -s $MADX_EXECUTABLE $1/madx
