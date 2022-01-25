#	UNIX shell file to source to make Motif version of Entrez using
#   cc compiler on Solaris.  Just modify this according to the header
#   info in makeall.unx for other systems.
#	
#   Substitute the path to your Motif includes and libraries.
#	

make -f makedemo.unx CC="cc -Xa" VIBLIBS="-R/usr/openwin/lib -R/opt/SUNWmotif/lib -lXm -lXmu -lXt -lX11 -lXext -lresolv -lsocket -lnsl -lgen" LIB4=-lvibrant LIB20=-lncbidesk VIBFLAG="-I/usr/openwin/include -I/opt/SUNWmotif/include -L/usr/openwin/lib -L/opt/SUNWmotif/lib -DWIN_MOTIF" VIB=entrez copy entrez OTHERLIBS="-lthread -lm"
