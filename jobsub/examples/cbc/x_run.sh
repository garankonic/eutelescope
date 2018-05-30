#!/bin/sh

#################################
#                               #
# CMS CBC Analysis - EUDAQ Data #
#                               #
# GBL alignment                 #
#                               #
#################################

# usage: sh x_run.sh runnumber

#jobsub.py -c config.cfg -csv runlist.csv telescope-converter $1
#jobsub.py -c config.cfg -csv runlist.csv combined-clustering $1
#jobsub.py -c config.cfg -csv runlist.csv combined-filter $1
#jobsub.py -c config.cfg -csv runlist.csv combined-hitmaker $1
#jobsub.py -c config.cfg -csv runlist.csv stub $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-1 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-2 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-3 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-4 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-5 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-6 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-7 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-8 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-9 $1
jobsub.py -c config_align.cfg -csv runlist.csv alignment-gbl-10 $1
jobsub.py -c config.cfg -csv runlist.csv tracking-gbl $1
jobsub.py -c config.cfg -csv runlist.csv cbc-recovering $1

################################
# event-viewer displays tracks
################################
################################
# Uncomment this for event display server:
################################

#	# check if glced is running, if not start it
#	SERVICE='glced'
#	if ps ax | grep -v grep | grep $SERVICE > /dev/null
#	then
#		echo "$SERVICE is running!"
#	else
#		echo "$SERVICE is not running, will start it!"
#	$SERVICE &
#	fi
#
#	jobsub.py -c config.cfg -csv runlist.csv event-viewer $1

#
