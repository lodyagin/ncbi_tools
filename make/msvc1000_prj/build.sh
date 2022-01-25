#! /bin/sh
# $Id: build.sh,v 1.8 2012/04/18 15:12:13 ivanov Exp $
# Author:  Vladimir Ivanov (ivanov@ncbi.nlm.nih.gov)
#
# Build C Toolkit on MSVC10 32/64.


#---------------- Arguments ----------------

script="$0"
cfgs="${1:-Debug Release DebugDLL ReleaseDLL}"
arch="$2"
target="${3:-all_ncbi}"


#---------------- Global variables ----------------

timer="date +'%H:%M'"
out=".build.$$"

# TRUE if parallell project build system is enabled in Visual Studio
is_ppb=false
need_ppb_check=true


#---------------- Functions ----------------

error()
{
    echo "[`basename $script`] ERROR:  $1"
    exit 1

}

generate_msvc10_error_check_file() {
    cat <<-EOF >$1
	/(| : |The source )([fatal ]*error [A-Z]*[0-9]* *: |The .* are both configured to produce|.*: error [0-9]*:|Error executing |ERROR: This project depends)/ {
	  print \$0
	  exit
	}
	EOF
}

generate_simple_log()
{
    echo Parallel project build detected! Creating simplified log.
    echo

    log=$1
    sol=$2
    cfg=$3

    # All 64-bit builds go under x64 subdirectories
    test "$arch" = 64  &&  cfg="x64/$cfg"

    # Get built projects
    projects=`grep '.*--* Build started:' $log | awk '{ sub(/^.* started:/, ""); gsub(/ /,"#"); print $0}'`

    for p in $projects ; do
        echo "------$p" | awk '{gsub(/[#]/," "); print}'
        prj_name=`echo $p | awk '{gsub(/[#,]/," "); print $2}'`

        # Get path for specified project name from solution
        s=`grep \"$prj_name\" $sol | awk '{gsub(/,/," "); print $4}' | sed -e 's%"%%g' -e 's%\\\%/%g' -e 's%.vcxproj%%'`

        target_dir=`echo $s | sed 's%/[^/]*$%%'`
        test $target_dir = $s  &&  target_dir=''
        target_name=`echo $s | sed 's%^.*/%%'`

        # Path to regular logfile for current project
        if [ -z "$target_dir" ]; then
           prj_log="${arch_dir}$cfg/$target_name.log"
        else
           prj_log="$target_dir/$cfg/$target_name.log"
        fi

        # Add it to new combined log
        if test ! -f "$prj_log" ; then
            echo "BUILD_SYSTEM_ERROR: Cannot find log file for this project: $prj_log"
            echo
            continue
        fi
        # Remove 3 first bytes from logfile (EF BB BF)
        cat $prj_log | tr -d '\357\273\277'
        echo
    done
    grep '.*========== Build:' $log
    echo
}



#---------------- Main ----------------

# Get build dir
build_dir=`dirname $script`
build_dir=`(cd "$build_dir"; pwd)`

if [ ! -d $build_dir ] ; then
   error "Build directory $build_dir not found"
   exit 1
fi
cd $build_dir


# Generate errors check script

check_awk=$build_dir/build_check.awk
generate_msvc10_error_check_file $check_awk


# Build

for cfg in $cfgs ; do
    start=`eval $timer`
    echo
    echo Start time: $start
    echo "INFO: Building \"$cfg\""
    echo
    $build_dir/build_exec.bat "ncbi.sln" build "$arch" "$cfg" "$target" $out >/dev/null
    status=$?
    if $need_ppb_check; then
       need_ppb_check=false
       grep '^1>------ Build started:' $out >/dev/null 2>&1  &&  is_ppb=true
    fi
    if $is_ppb; then
       generate_simple_log $out ncbi.sln $cfg > $out.simple
       mv $out $cfg.log
       mv $out.simple $out
    fi 
    cat $out
    echo "Build time: $start - `eval $timer`"
    echo STATUS = $status

    if [ $status -ne 0 ] ; then
       # Check on errors
       failed="1"
       awk -f $check_awk $out >$out.res 2>/dev/null  &&  test ! -s $out.res  &&  failed="0"
       rm -f $out $out.res >/dev/null 2>&1
       if [ "$failed" = "1" ]; then
          echo FAILED: Build $cfg
          exit 4
       fi
    fi
    rm -f $out >/dev/null 2>&1
done


exit 0
