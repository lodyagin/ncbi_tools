#! /bin/sh
# $Id: install.sh,v 1.26 2003/01/31 14:45:06 lavr Exp $
# Authors:  Denis Vakatov    (vakatov@ncbi.nlm.nih.gov)
#           Vladimir Ivanov  (ivanov@ncbi.nlm.nih.gov)
#           Anton Lavrentiev (lavr@ncbi.nlm.nih.gov)
#
# Deploy sources, headers, libraries and executables for the further use
# by the "external" users' projects

# Cmd.-line args  -- source and destination
script="$0"
builddir="${1:-//u/coremake/ncbi}"
target="${2:-//u/coremake/public/ncbi}"

if test -n "$3" ; then
  echo "USAGE:  `basename $script` [build_dir] [install_dir]"
fi


error()
{
  echo "[`basename $script`] ERROR:  $1"
  exit 1
}


makedir()
{
  test -d "$1"  ||  mkdir $2 "$1"  ||  error "Cannot create \"$1\""
}


echo "[`basename $script`] NCBI C:  \"$builddir\" to \"$target\"..."
sleep 2


# Derive the destination dirs
incdir="$target"/include
srcdir="$target"/src
dbgdir="$target"/Debug
libdir="$target"/Release
bindir="$target"/bin
datdir="$target"/data

# Alternate dirs (mirrors)
srcdir_a="$target"/altsrc
dbgdir_a="$target"/dbglib
libdir_a="$target"/lib


# Check
test -d "$builddir"  ||  error "Absent build dir \"$builddir\""


# Reset the public directory
test -d "$target"  &&  find "$target" -type f -exec rm -f {} \; >/dev/null 2>&1
test -d "$target"  ||  mkdir -p "$target"
test -d "$target"  ||  error "Cannot create target dir \"$target\""


# Make all dirs
makedir "$incdir" -p
makedir "$incdir"/connect
makedir "$incdir"/ctools
makedir "$srcdir" -p
makedir "$srcdir"/access
makedir "$srcdir"/api
makedir "$srcdir"/asnlib
makedir "$srcdir"/biostruc
makedir "$srcdir"/biostruc/cdd
makedir "$srcdir"/biostruc/cn3d
makedir "$srcdir"/cdromlib
makedir "$srcdir"/cn3d
makedir "$srcdir"/corelib
makedir "$srcdir"/connect
makedir "$srcdir"/connect/test
makedir "$srcdir"/ctools
makedir "$srcdir"/ddv
makedir "$srcdir"/demo
makedir "$srcdir"/desktop
makedir "$srcdir"/gif
makedir "$srcdir"/network/blast3/client  -p
makedir "$srcdir"/network/entrez/client  -p
makedir "$srcdir"/network/id1arch        -p
makedir "$srcdir"/network/medarch/client -p
makedir "$srcdir"/network/nsclilib       -p
makedir "$srcdir"/network/nsdemocl       -p
makedir "$srcdir"/network/taxon1/common  -p
makedir "$srcdir"/network/taxon1/taxon2  -p
makedir "$srcdir"/network/vibnet         -p
makedir "$srcdir"/object
makedir "$srcdir"/sequin
makedir "$srcdir"/tools
makedir "$srcdir"/vibrant
for i in '' 'DLL' ; do
  makedir "$dbgdir$i" -p
  makedir "$libdir$i" -p
done
makedir "$bindir"   -p
makedir "$datdir"   -p
makedir "$srcdir_a" -p
makedir "$dbgdir_a" -p
makedir "$libdir_a" -p


# Copy files

bd="$builddir"

cp -p "$bd"/corelib/tsprintf.h            "$incdir"
cp -p "$bd"/corelib/gifgen.h              "$incdir"
cp -p "$bd"/corelib/ncbi*.h               "$incdir"
cp -p "$bd"/corelib/tree*.h               "$incdir"
cp -p "$bd"/corelib/matrix.h              "$incdir"
cp -p "$bd"/corelib/binary.h              "$incdir"
cp -p "$bd"/corelib/*.c                   "$srcdir"/corelib
cp -p "$bd"/corelib/core*.h               "$srcdir"/corelib
cp -p "$bd"/corelib/regex.h               "$incdir"
cp -p "$bd"/asnlib/*.h                    "$srcdir"/asnlib
mv "$srcdir"/asnlib/asn.h                 "$incdir"
cp -p "$bd"/asnlib/*.c                    "$srcdir"/asnlib
cp -p "$bd"/connect/*.c                   "$srcdir"/connect
cp -p "$bd"/connect/ncbi_priv.h           "$srcdir"/connect
cp -p "$bd"/connect/ncbi_comm.h           "$srcdir"/connect
cp -p "$bd"/connect/ncbi_host_infop.h     "$srcdir"/connect
cp -p "$bd"/connect/ncbi_server_infop.h   "$srcdir"/connect
cp -p "$bd"/connect/ncbi_servicep.h       "$srcdir"/connect
cp -p "$bd"/connect/ncbi_dispd.h          "$srcdir"/connect
cp -p "$bd"/connect/ncbi_lbsmd.h          "$srcdir"/connect
cp -p "$bd"/connect/*.h                   "$incdir"
mv "$incdir"/ncbi_*.h                     "$incdir"/connect
mv "$incdir"/connect_export.h             "$incdir"/connect
cp -p "$bd"/connect/test/*.[ch]           "$srcdir"/connect/test
cp -p "$bd"/ctools/*.c                    "$srcdir"/ctools
cp -p "$bd"/ctools/*.h                    "$incdir"/ctools
cp -p "$bd"/object/*.c                    "$srcdir"/object
cp -p "$bd"/object/*.h                    "$incdir"
cp -p "$bd"/access/*.c                    "$srcdir"/access
cp -p "$bd"/access/*.h                    "$incdir"
cp -p "$bd"/asnstat/*.h                   "$incdir"
cp -p "$bd"/api/*.c                       "$srcdir"/api
cp -p "$bd"/api/*.h                       "$incdir"
cp -p "$bd"/cdromlib/*.c                  "$srcdir"/cdromlib
cp -p "$bd"/cdromlib/*.h                  "$incdir"
cp -p "$bd"/biostruc/*.c                  "$srcdir"/biostruc
cp -p "$bd"/biostruc/*.h                  "$incdir"
cp -p "$bd"/biostruc/cdd/*.c              "$srcdir"/biostruc/cdd
cp -p "$bd"/biostruc/cdd/*.h              "$incdir"
cp -p "$bd"/biostruc/cn3d/*.c             "$srcdir"/biostruc/cn3d
cp -p "$bd"/biostruc/cn3d/*.h             "$incdir"
cp -p "$bd"/tools/*.c                     "$srcdir"/tools
cp -p "$bd"/tools/*.h                     "$incdir"
cp -p "$bd"/link/mswin/ncbirc.h           "$srcdir"
cp -p "$bd"/vibrant/*.c                   "$srcdir"/vibrant
cp -p "$bd"/vibrant/*.h                   "$incdir"
cp -p "$bd"/desktop/*.c                   "$srcdir"/desktop
cp -p "$bd"/desktop/*.h                   "$incdir"
cp -p "$bd"/gif/*.c                       "$srcdir"/gif
cp -p "$bd"/gif/*.h                       "$incdir"
cp -p "$bd"/cn3d/*.c                      "$srcdir"/cn3d
cp -p "$bd"/cn3d/*.h                      "$incdir"
cp -p "$bd"/ddv/*.c                       "$srcdir"/ddv
cp -p "$bd"/ddv/*.h                       "$incdir"


# Copy network files

cp -p "$bd"/network/entrez/client/*.c     "$srcdir"/network/entrez/client
cp -p "$bd"/network/entrez/client/*.h     "$incdir"
cp -p "$bd"/network/nsclilib/*.c          "$srcdir"/network/nsclilib
cp -p "$bd"/network/nsclilib/*.h          "$incdir"
cp -p "$bd"/network/medarch/client/*.c    "$srcdir"/network/medarch/client
cp -p "$bd"/network/medarch/client/*.h    "$incdir"
cp -p "$bd"/network/taxon1/common/*.c     "$srcdir"/network/taxon1/common
cp -p "$bd"/network/taxon1/common/*.h     "$incdir"
cp -p "$bd"/network/taxon1/taxon2/*.c     "$srcdir"/network/taxon1/taxon2
cp -p "$bd"/network/taxon1/taxon2/*.h     "$incdir"
cp -p "$bd"/network/blast3/client/*.c     "$srcdir"/network/blast3/client
cp -p "$bd"/network/blast3/client/*.h     "$incdir"
cp -p "$bd"/network/id1arch/*.c           "$srcdir"/network/id1arch
cp -p "$bd"/network/id1arch/*.h           "$incdir"
cp -p "$bd"/network/nsdemocl/*.[hc]       "$srcdir"/network/nsdemocl
cp -p "$bd"/demo/entrez.c                 "$srcdir"/demo
cp -p "$bd"/demo/entrezcf.c               "$srcdir"/demo
cp -p "$bd"/demo/netentcf.c               "$srcdir"/demo
cp -p "$bd"/demo/ccpv.c                   "$srcdir"/demo
cp -p "$bd"/demo/ccp.c                    "$srcdir"/demo
cp -p "$bd"/demo/dustv.c                  "$srcdir"/demo
cp -p "$bd"/demo/dst.c                    "$srcdir"/demo
cp -p "$bd"/demo/epiv.c                   "$srcdir"/demo
cp -p "$bd"/demo/epi.c                    "$srcdir"/demo
cp -p "$bd"/demo/sigmev.c                 "$srcdir"/demo
cp -p "$bd"/demo/sigme.c                  "$srcdir"/demo
cp -p "$bd"/demo/searchv.c                "$srcdir"/demo
cp -p "$bd"/demo/srchaa.c                 "$srcdir"/demo
cp -p "$bd"/demo/srchnt.c                 "$srcdir"/demo
cp -p "$bd"/demo/twopv.c                  "$srcdir"/demo
cp -p "$bd"/demo/twop.c                   "$srcdir"/demo
cp -p "$bd"/demo/cnsrtv.c                 "$srcdir"/demo
cp -p "$bd"/demo/cnsrt.c                  "$srcdir"/demo
cp -p "$bd"/demo/cnsgnv.c                 "$srcdir"/demo
cp -p "$bd"/demo/cnsgn.c                  "$srcdir"/demo
cp -p "$bd"/demo/udvmain.c                "$srcdir"/demo
cp -p "$bd"/sequin/*.c                    "$srcdir"/sequin
cp -p "$bd"/sequin/*.h                    "$incdir"
cp -p "$bd"/network/vibnet/*.c            "$srcdir"/network/vibnet
cp -p "$bd"/network/vibnet/*.h            "$incdir"


# Object files
for i in '' 'DLL' ; do
  cp -p "$bd"/make/msvc_prj/corelib/ncbimain/"Debug$i"/ncbimain.obj   "$dbgdir$i"
  cp -p "$bd"/make/msvc_prj/corelib/ncbi/"Debug$i"/ncbithr.obj        "$dbgdir$i"
  cp -p "$bd"/make/msvc_prj/corelib/ncbimain/"Release$i"/ncbimain.obj "$libdir$i"
  cp -p "$bd"/make/msvc_prj/corelib/ncbi/"Release$i"/ncbithr.obj      "$libdir$i"
done

for i in '' 'DLL' ; do
  # Debug libs
  cp -p `find $buiddir -name '*.lib' | grep "Debug$i/"` "$dbgdir$i"

  # Release libs
  cp -p `find $buiddir -name '*.lib' | grep "Release$i/"` "$libdir$i"
done


# Executables
cp -p `find $buiddir -name '*.exe' | grep "ReleaseDLL/"` "$bindir"


# Data
cp -pr "$builddir"/data/* "$datdir"


# Fill alt source tree
cd "$srcdir"
x_dirs=`find . -type d -print`
for ddd in $x_dirs ; do
  cd "$srcdir/$ddd"
  cp -p *.[ch] "$srcdir_a"
done


# Fill mirror dirs
cp -p "${dbgdir}DLL"/*.* "$dbgdir_a"
cp -p "${libdir}DLL"/*.* "$libdir_a"


exit 0
