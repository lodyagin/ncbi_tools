#	UNIX shell file to source to make Vibrant version of tools using
#   cc compiler on Solaris.  Just modify this according to the header
#   info in makeall.unx for other systems.
#	
#   Substitute the path to your Motif includes and libraries.
#	

make LCL=sol CC="cc -Xa"  VIBLIBS="-R/usr/openwin/lib -R/opt/SUNWmotif/lib -lXm -lXmu -lXt -lX11 -lXext" LIB30=libncbicn3d.a LIB28=libvibgif.a LIB4=libvibrant.a LIB20=libncbidesk.a LIB20=libncbidesk.a VIBFLAG="-I/usr/openwin/include -I/opt/SUNWmotif/include -L/usr/openwin/lib -L/opt/SUNWmotif/lib -DWIN_MOTIF" OTHERLIBS="-lthread -lm"

