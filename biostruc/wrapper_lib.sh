#! /bin/sh
#
# $Id: wrapper_lib.sh,v 6.2 1999/04/29 23:08:31 kimelman Exp $
#
# this is CGI handler wrapper library. It works as a membrane between httpd and actual
# cgi program and allow to run new technological version of such program in
# parallel to old one and compare results in "real life" condition.
#

atexit() {
    rm -rf $TMPdir
}

basic_settings() {
    basename=`basename $progname | sed 's/[.][^.]*//g'`
    EXECs_to_try="./$basename.REAL ./$basename.NEW ./$basename.OLD"
    TMPtop=/tmp/$basename
    STATs_to_try="./stats.$basename log/stats.$basename"
    
    if [ "x$THEvictim" = x ]; then
        THEvictim="kimelman"
#        THEvictim="kimelman zimmerma"
    fi
    while [ -d $TMPtop -a ! -w $TMPtop ]; 
    do
        TMPtop="${TMPtop}.A"
    done
    if [ ! -d $TMPtop ]; then
        mkdir -p $TMPtop
        chmod a+rw $TMPtop
    fi
    TMPdir=$TMPtop/$$
    awk_prg=$TMPtop/awkp
    if [ ! -d $TMPdir ]; then
        mkdir -p $TMPdir 
        chmod a+rw $TMPdir
    fi

    trap atexit 0,1
}

res_name() {
  res_fname=${TMPdir}/`echo "$1" | tr '/.' '__'`
  dir=`dirname $res_fname`
  [ -d  $dir ] || mkdir -p $dir
}

run_cgi() {
  res_name $1
  case $1 in 
  /* ) 
      prg="$1"
      ;;
  *) 
      prg="$progdir/$1"
      ;;
  esac
  [ x$stats = x ] || timep="/usr/bin/time -p"
  (
     cd `dirname $prg`
     $timep ./`basename $prg` $options <$input_file >$res_fname 2>$res_fname.2
  )  &
  waitedID="$waitedID $!"
}

proc_res() {
  res_name $1
  fn=`basename $res_fname`
  # echo "proc_res: $1 --> $fn"
  if [ "x$rcode" != x -a $rcode -eq 0 -a x$stats != x ]; then
    timing
  fi
}

prestart_checks() {
  httpdenv_check_and_fix
  check_executables
  # check sybase settings
  if [ -f ${progdir}/st_configure.sh ] ; then
    if [ x$stats != x ] ; then
        sybase_settings_fname=`dirname $stats`/.syb_set
    fi
    . ${progdir}/st_configure.sh
  fi
}

run_all() {
    # stage A : store input stream
    res_name input
    input_file=$res_fname
    cat >$input_file

    #
    # run all cgi programs in parallel
    #

    for texec in $EXECs ; do
        run_cgi $texec
    done

    #
    # wait for termination of all of them
    #
    main_result=

    for execname in $EXECs ; do # for every runned version of sgi program
        # get it's process id
        set -- $waitedID
        wid=$1 ; shift ; waitedID="$*"
        
        wait $wid            # wait for it to complete
        rcode=$?             # end get it return code and 
        proc_res $execname   # and process the result
        if [ $rcode -ne 0 ] ; then 
            # report mmdbsrv.REAL failure
            report_failure $execname $rcode $res_fname $res_fname.2 
            continue
        fi
        if [ x$main_result = x ] ; then
            main_result=$res_fname
            # output result of first program
            cat $main_result
            continue
        fi
        # here we have: a: 0 result code & more than 1 output 
        # compare results now
        diff -c $main_result $res_fname >$res_fname.diff
        if [ $? -gt 0 ] ; then
            report_failure $execname $rcode $res_fname.diff $res_fname.2
        fi
    done
}

#
# gather timing statistics.
#

do_timing() {
    if [ x$stats = x ] ; then
        for stats_file in $STATs_to_try ; do
            if [ -w `dirname $stats_file` ]; then
                stats=$stats_file
                break
            fi
        done
    fi
}

create_stats_prog() {
cat >$awk_prg <<EOF
BEGIN   { real=0 ; user=0 ; sys=0 ; cnt=0 ; }
/real/  { real=\$2 }
/user/  { user=\$2 }
/sys/   { sys=\$2  }
        { if ( \$1 == fn ) { cnt=\$5 ; real+=\$2 * cnt ;  user+=\$3 *cnt ; sys +=\$4 * cnt ;  } }
END     { cnt+=1 ; print fn, real/cnt, user/cnt, sys/cnt, cnt }
EOF
}

timing() {
    dd=$$
    while [ -f $stats.lock ] ; do sleep 1 ; done
    echo $dd >$stats.lock
    if [ -f $stats.lock -a "$dd" != "`cat $stats.lock`" ] ; then
        report_lock_problem
    fi
    [ -f $stats ] || echo '! Stats file' >$stats
    grep -v "^$fn " $stats > $res_fname.3
    grep "^$fn " $stats | cat $res_fname.2 -  >$res_fname.4
    if [ ! -r $awk_prg ] ; then
        create_stats_prog
    fi
    nawk -v fn=$fn -f $awk_prg $res_fname.4 >>$res_fname.3
    rm -f $stats
    cat $res_fname.3 >$stats
    chmod a+w $stats
    rm $res_fname.3
    [ ! -f ${stats}.lock ] || rm ${stats}.lock
}

#
#  CHECK & REPORT functions  
#
#

# get aname of person who should got bug reports
get_victim() { # (filename)
    victim=
    [ "x$1" = x ] || victim="`ls -l $1 | awk '  { print $3 ; }'`"
    [ "x$victim" != xpubmed ] || victim=""
    [ "x$victim" != x ] || [ ! -r $1.recepient ] || victim="`cat $1.recepients`"
    [ "x$victim" != x ] || [ ! -r .mail_recepients ] || victim="`cat .mail_recepients`"
    [ "x$victim" != x ] || victim="$THEvictim"
}

httpdenv_check_and_fix() { # check environment & fix required but unset env vars.
    # especially useful in case of debug environment
    if [ "x$REMOTE_ADDR" = x ] ; then
        REMOTE_ADDR=127.0.0.1
        export REMOTE_PORT
    fi
    if [ "x$REMOTE_HOST" = x ] ; then
        REMOTE_HOST=localhost
        export REMOTE_HOST
    fi
    if [ "x$REMOTE_PORT" = x ] ; then
        REMOTE_PORT=0
        export REMOTE_ADDR
    fi
    if [ "x$HTTP_USER_AGENT" = x ] ; then
        HTTP_USER_AGENT="pseudo agent : $progname"
        export HTTP_USER_AGENT
    fi
    if [ "x$REQUEST_URI" = x ] ; then
        REQUEST_URI=$progname
        export REQUEST_URI
    fi
    if [ "x$QUERY_STRING" = x ] ; then
        QUERY_STRING=''
        export QUERY_STRING
    fi
}

check_executables() { # check & find executable from the list
    EXECs_given=$EXECs
    EXECs=
    if [ "x$EXECs_given" != x ] ; then
        for fexec in $EXECs_given ; do
            [ ! -x $fexec ] || EXECs="$EXECs $fexec"
        done
    else
        for fexec in $EXECs_to_try ; do
            [ ! -x $fexec ] || EXECs="$EXECs $fexec"
        done
    fi
    if [ "x$EXECs" = x ]; then
        get_victim
        mail $victim <<EOF
Subject ${progname} : can find binaries to run

`ls -al`

EOF
        exit 1
    fi
}

report_lock_problem() {
     get_victim
     opid=`cat $stats.lock`
     mail $victim <<EOF
Subject: mmdbsrv.wrapper locking problem

dd="$dd"
$stats.lock=$opid

`ps -ef | grep $opid`

`ls -al`

`ps -ef`

`set`

EOF
     exit 1
}

report_failure() { # $execname $rcode $res_fname $res_fname.2
    rf_execn=$1
    rf_code=$2 
    rf_rn=$3
    rf_rn2=$4
    get_victim $rf_execn
    mail $victim <<EOF
Subject: $rf_execn failed with status $rf_code

HTTP_HOST=${HTTP_HOST}
HTTP_REFERER=${HTTP_REFERER}
QUERY_STRING=${QUERY_STRING}
REMOTE_HOST=${REMOTE_HOST}
REQUEST_METHOD=${REQUEST_METHOD}
REQUEST_URI=${REQUEST_URI}
SYBASE=${SYBASE}
COMMAND_LINE=${options}

====================================================================
Input file:
`cat $input_file`

====================================================================
Log.2:
`cat $rf_rn2`

Log:
`cat $rf_rn`

====================================================================
----------- Environment:----------------

`env`

EOF
}
