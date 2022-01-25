#! /bin/sh
# $Id: install.sh,v 1.15 2002/07/26 15:13:41 lavr Exp $
# Author:  Denis Vakatov (vakatov@ncbi.nlm.nih.gov)
#          Vladimir Ivanov (ivanov@ncbi.nlm.nih.gov)
#
# Deploy sources, headers, libraries and executables for the further use
# by the "external" users' projects

# Cmd.-line args  -- source and destination

builddir="$1"
target="$2"

if test -z "$target" ; then
  echo "USAGE:  `basename $0` [build_dir] [install_dir]"
fi

builddir="${builddir:-//u/coremake/ncbi}"
target="${target:-//u/coremake/public/ncbi}"

echo "[`basename $0`]  NCBI C:   \"$builddir\" to \"$target\"..."


# Derive the destination dirs
incdir="$target/include"
srcdir="$target/src"
dbgdir="$target/Debug"
libdir="$target/Release"
bindir="$target/bin"
datdir="$target/data"

# Alternate dirs (mirrors)
srcdir_a="$target/altsrc"
dbgdir_a="$target/dbglib"
libdir_a="$target/lib"


# Check
if test ! -d "$builddir" ; then
  echo "[$0] ERROR:  Absent build dir \"$builddir\""
  exit 1
fi

# Reset the public directory
if test -d "$target" ; then
  rm -rf "$target"
fi
mkdir -p "$target"

# Make dirs without mirrors
mkdir -p "$incdir"
mkdir -p "$incdir/connect"
mkdir -p "$incdir/ctools"
mkdir -p "$srcdir"
for i in "" "MT"; do
  mkdir -p "$dbgdir$i"
  mkdir -p "$libdir$i"
done
mkdir -p "$bindir"
mkdir -p "$datdir"


# Copy files

bd="$builddir"

cp -p "$bd"/corelib/tsprintf.h            "$incdir"
cp -p "$bd"/corelib/gifgen.h              "$incdir"
cp -p "$bd"/corelib/ncbi*.h               "$incdir"
cp -p "$bd"/corelib/tree*.h               "$incdir"
cp -p "$bd"/corelib/matrix.h              "$incdir"
cp -p "$bd"/corelib/binary.h              "$incdir"
mkdir "$srcdir/corelib"
cp -p "$bd"/corelib/*.c                   "$srcdir/corelib"
cp -p "$bd"/corelib/core*.h               "$srcdir/corelib"
cp -p "$bd"/corelib/regex.h               "$incdir"
mkdir "$srcdir/asnlib"
cp -p "$bd"/asnlib/*.h                    "$srcdir/asnlib"
mv "$srcdir/asnlib/asn.h"                 "$incdir"
cp -p "$bd"/asnlib/*.c                    "$srcdir/asnlib"
mkdir "$srcdir/connect"
cp -p "$bd"/connect/*.c                   "$srcdir/connect"
cp -p "$bd"/connect/ncbi_priv.h           "$srcdir/connect"
cp -p "$bd"/connect/ncbi_comm.h           "$srcdir/connect"
cp -p "$bd"/connect/ncbi_server_infop.h   "$srcdir/connect"
cp -p "$bd"/connect/ncbi_servicep.h       "$srcdir/connect"
cp -p "$bd"/connect/ncbi_servicep_dispd.h "$srcdir/connect"
cp -p "$bd"/connect/ncbi_servicep_lbsmd.h "$srcdir/connect"
cp -p "$bd"/connect/*.h                   "$incdir"
mv "$incdir"/ncbi_*.h                     "$incdir/connect"
mkdir "$srcdir/connect/test"
cp -p "$bd"/connect/test/*.[ch]           "$srcdir/connect/test"
mkdir "$srcdir/ctools"
cp -p "$bd"/ctools/*.c                    "$srcdir/ctools"
cp -p "$bd"/ctools/*.h                    "$incdir/ctools"
mkdir "$srcdir/object"
cp -p "$bd"/object/*.c                    "$srcdir/object"
cp -p "$bd"/object/*.h                    "$incdir"
mkdir "$srcdir/access"
cp -p "$bd"/access/*.c                    "$srcdir/access"
cp -p "$bd"/access/*.h                    "$incdir"
cp -p "$bd"/asnstat/*.h                   "$incdir"
mkdir "$srcdir/api"
cp -p "$bd"/api/*.c                       "$srcdir/api"
cp -p "$bd"/api/*.h                       "$incdir"
mkdir "$srcdir/cdromlib"
cp -p "$bd"/cdromlib/*.c                  "$srcdir/cdromlib"
cp -p "$bd"/cdromlib/*.h                  "$incdir"
mkdir "$srcdir/biostruc"
cp -p "$bd"/biostruc/*.c                  "$srcdir/biostruc"
cp -p "$bd"/biostruc/*.h                  "$incdir"
mkdir "$srcdir/biostruc/cdd"
cp -p "$bd"/biostruc/cdd/*.c              "$srcdir/biostruc/cdd"
cp -p "$bd"/biostruc/cdd/*.h              "$incdir"
mkdir "$srcdir/biostruc/cn3d"
cp -p "$bd"/biostruc/cn3d/*.c             "$srcdir/biostruc/cn3d"
cp -p "$bd"/biostruc/cn3d/*.h             "$incdir"
mkdir "$srcdir/tools"
cp -p "$bd"/tools/*.c                     "$srcdir/tools"
cp -p "$bd"/tools/*.h                     "$incdir"
cp -p "$bd"/link/mswin/ncbirc.h           "$srcdir"
mkdir "$srcdir/vibrant"
cp -p "$bd"/vibrant/*.c                   "$srcdir/vibrant"
cp -p "$bd"/vibrant/*.h                   "$incdir"
mkdir "$srcdir/desktop"
cp -p "$bd"/desktop/*.c                   "$srcdir/desktop"
cp -p "$bd"/desktop/*.h                   "$incdir"
mkdir "$srcdir/gif"
cp -p "$bd"/gif/*.c                       "$srcdir/gif"
cp -p "$bd"/gif/*.h                       "$incdir"
mkdir "$srcdir/cn3d"
cp -p "$bd"/cn3d/*.c                      "$srcdir/cn3d"
cp -p "$bd"/cn3d/*.h                      "$incdir"
mkdir "$srcdir/ddv"
cp -p "$bd"/ddv/*.c                       "$srcdir/ddv"
cp -p "$bd"/ddv/*.h                       "$incdir"


# Copy network files

mkdir -p "$srcdir/network/entrez/client"
cp -p "$bd"/network/entrez/client/*.c     "$srcdir/network/entrez/client"
cp -p "$bd"/network/entrez/client/*.h     "$incdir"
mkdir -p "$srcdir/network/nsclilib"
cp -p "$bd"/network/nsclilib/*.c          "$srcdir/network/nsclilib"
cp -p "$bd"/network/nsclilib/*.h          "$incdir"
mkdir -p "$srcdir/network/medarch/client"
cp -p "$bd"/network/medarch/client/*.c    "$srcdir/network/medarch/client"
cp -p "$bd"/network/medarch/client/*.h    "$incdir"
mkdir -p "$srcdir/network/taxon1/common"
cp -p "$bd"/network/taxon1/common/*.c     "$srcdir/network/taxon1/common"
cp -p "$bd"/network/taxon1/common/*.h     "$incdir"
mkdir -p "$srcdir/network/taxon1/taxon2"
cp -p "$bd"/network/taxon1/taxon2/*.c     "$srcdir/network/taxon1/taxon2"
cp -p "$bd"/network/taxon1/taxon2/*.h     "$incdir"
mkdir -p "$srcdir/network/blast3/client"
cp -p "$bd"/network/blast3/client/*.c     "$srcdir/network/blast3/client"
cp -p "$bd"/network/blast3/client/*.h     "$incdir"
mkdir -p "$srcdir/network/id1arch"
cp -p "$bd"/network/id1arch/*.c           "$srcdir/network/id1arch"
cp -p "$bd"/network/id1arch/*.h           "$incdir"
mkdir -p "$srcdir/network/nsdemocl"
cp -p "$bd"/network/nsdemocl/*.[hc]       "$srcdir/network/nsdemocl"
mkdir "$srcdir/demo"
cp -p "$bd"/demo/entrez.c                 "$srcdir/demo"
cp -p "$bd"/demo/entrezcf.c               "$srcdir/demo"
cp -p "$bd"/demo/netentcf.c               "$srcdir/demo"
cp -p "$bd"/demo/ccpv.c                   "$srcdir/demo"
cp -p "$bd"/demo/ccp.c                    "$srcdir/demo"
cp -p "$bd"/demo/dustv.c                  "$srcdir/demo"
cp -p "$bd"/demo/dst.c                    "$srcdir/demo"
cp -p "$bd"/demo/epiv.c                   "$srcdir/demo"
cp -p "$bd"/demo/epi.c                    "$srcdir/demo"
cp -p "$bd"/demo/sigmev.c                 "$srcdir/demo"
cp -p "$bd"/demo/sigme.c                  "$srcdir/demo"
cp -p "$bd"/demo/searchv.c                "$srcdir/demo"
cp -p "$bd"/demo/srchaa.c                 "$srcdir/demo"
cp -p "$bd"/demo/srchnt.c                 "$srcdir/demo"
cp -p "$bd"/demo/twopv.c                  "$srcdir/demo"
cp -p "$bd"/demo/twop.c                   "$srcdir/demo"
cp -p "$bd"/demo/cnsrtv.c                 "$srcdir/demo"
cp -p "$bd"/demo/cnsrt.c                  "$srcdir/demo"
cp -p "$bd"/demo/cnsgnv.c                 "$srcdir/demo"
cp -p "$bd"/demo/cnsgn.c                  "$srcdir/demo"
cp -p "$bd"/demo/udvmain.c                "$srcdir/demo"
mkdir "$srcdir/sequin"
cp -p "$bd"/sequin/*.c                    "$srcdir/sequin"
cp -p "$bd"/sequin/*.h                    "$incdir"
mkdir -p "$srcdir/network/vibnet"
cp -p "$bd"/network/vibnet/*.c            "$srcdir/network/vibnet"
cp -p "$bd"/network/vibnet/*.h            "$incdir"


# Object files
for i in "" "MT"; do
  cp -p "$bd/make/msvc_prj/corelib/ncbimain/Debug$i/ncbimain.obj" "$dbgdir$i"
  cp -p "$bd/make/msvc_prj/corelib/ncbimain/Release$i/ncbimain.obj" "$libdir$i"
  cp -p "$bd/make/msvc_prj/corelib/ncbi/Debug$i/ncbithr.obj" "$dbgdir$i"
  cp -p "$bd/make/msvc_prj/corelib/ncbi/Release$i/ncbithr.obj" "$libdir$i"
done

for i in "" "MT"; do
  # Debug libs
  cp -p `find $buiddir -name '*.lib' | grep "Debug${i}/"` "$dbgdir$i"

  # Release libs
  cp -p `find $buiddir -name '*.lib' | grep "Release${i}/"` "$libdir$i"
done


# Executables
cp -p `find $buiddir -name '*.exe' | grep "ReleaseMT/"` "$bindir"


# Data
cp -pr "$builddir/data/" "$target"


# Make mirrors dirs
cp -pr "${dbgdir}MT" "$dbgdir_a"
cp -pr "${libdir}MT" "$libdir_a"


# Alt source tree
mkdir "$srcdir_a"
cd "$srcdir"
x_dirs=`find . -type d -print`
for ddd in $x_dirs ; do
  cd "$srcdir/$ddd"
  cp -p *.[ch] "$srcdir_a"
done
