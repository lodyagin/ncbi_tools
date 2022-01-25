#! /bin/sh
#
# $Id: wrapper.sh,v 6.10 1999/04/19 17:32:30 kimelman Exp $
#
# this is CGI handler wrapper. It works as a membrane between httpd and actual
# cgi program and allow to run new technological version of such program in 
# parallel to old one and compare results in "real life" condition. 
#
# $Log: wrapper.sh,v $
# Revision 6.10  1999/04/19 17:32:30  kimelman
# new wrapper
#
# Revision 6.9  1998/10/07 21:58:16  kimelman
# code cleaned & madeeasy to read.
# check for directory access permission added
# binaries expected names changed
# this version of script is expected to work on public servers
#
# Revision 6.8  1998/06/09 19:16:07  kimelman
# extra bracket removed.
#
# Revision 6.7  1998/06/09 18:50:35  kimelman
# permissions problems workaround added to wrapper
#
#

progname=$0
options="$*"
progdir=`dirname $progname`

. $progdir/wrapper_lib.sh

basic_settings # most of environment - default binary name etc

# here you can override default values by your lovely ones
# EXECs_to_try="./mmdbsrv.REAL ./mmdbsrv.NEW new/mmdbsrv.REAL"
# EXECs="./mmdbsrv.REAL ./mmdbsrv.NEW new/mmdbsrv.REAL"
# and here you can define the victim
# THEvictim="kimelman zimmerma"

# comment out the next line if you don't want timing info archiving
do_timing

prestart_checks # check env , binaries , sybase settings 

run_all # run EXECs in parallel - diff results , gather time statistics, report deltas to 

#
# The end
#
exit 0
