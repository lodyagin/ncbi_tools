  Summary on the new configuration resources for NCBI network clients
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      By Denis V. Vakatov NCBI/NLM/NIH (vakatov@ncbi.nlm.nih.gov)


Due to the recent works on the old NCBI dispatcher replacement
I introduced some new resources which affect the new NCBI network
client functionality.


Supported network interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SRV_CONN_MODE -- to specify which one of the now available network
                 interfaces should be used.
Value:
   WWW        -- Web-based dispatcher in the stateful connection mode [default]
   DISPATCHER -- old NCBI dispatcher
   FIREWALL   -- just like "WWW" but the data connection goes through the
                 NCBI firewall daemon that listens at a well-known host/port


Resources for the DISPATCHER interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See in "???".


Resources for WWW and FIREWALL interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(NOTE):  In most cases, the default values should be fine;  we do not recommend
         you to modify them unless you have real problems(e.g. caused by a
         slow or unstable network connection) connecting to the NCBI services.


SRV_DEBUG_PRINTOUT:  to printout messages sent by dispatcher when
                     establishing client-server connection
Value:
   "1", "yes", "true" (case-insensitive) -- do printout
   otherwise -- do not printout [default]


SRV_CONN_TIMEOUT:  to set timeout for the connection establishment
Value:
   a positive integer or floating-point value(in sec.)  [default = 30]


SRV_CONN_TRY:  maximum number of attempts to establish the connection
Value:
   a positive integer  [default = 3]


How to set the resource values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The above resource values can be set either by using the NCBI resource
(configuration) file or (on UNIX, MS Windows and VMS) by setting relevant
environment variable.  The environment variable overrides the value
set in the NCBI resource file.

[NET_SERV]
SRV_CONN_MODE
SRV_DEBUG_PRINTOUT
SRV_CONN_TIMEOUT
SRV_CONN_TRY

For example, in order to set the connection establishment timeout to
12.345 sec. you should set environment variable:

UNIX csh:
   setenv SRV_CONN_TIMEOUT 12.345
UNIX sh, bash:
   SRV_CONN_TIMEOUT=12.345
   export SRV_CONN_TIMEOUT
DOS prompt:
   set SRV_CONN_TIMEOUT=12.345

or, in the NCBI configuration file:

   [NET_SERV]
   SRV_CONN_TIMEOUT=12.345


How to access NCBI services from behind a firewall
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the FIREWALL network interface and ask your system administrator
to configure your firewall for hosts/ports advertised in "???".


$Date: 1998/09/09 22:38:01 $
