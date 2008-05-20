FUNCTION getdata, snapshot, wkdir=wkdir,nzeros=nzeros,_EXTRA=extra

on_error, 2
IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'
IF (N_ELEMENTS(nzeros) EQ 0) THEN nzeros=4

IF N_PARAMS() EQ 0 THEN BEGIN
    print, "Usage: result = getdata(snapnumber [,wkdir=<dir>, /empty | /rho, /temp, /vx ...])"
    RETURN, "Usage: result = getdata(snapnumber [,wkdir=<dir>, /empty | /rho, /temp, /vx ...])"
ENDIF
string1='("/",'
string2=',".cfd")'
string_desc=string1+string(nzeros,format="('I',I03)")+string2
file = wkdir + string(snapshot,format=string_desc)

RETURN, LoadCFDFile(file,_EXTRA=extra)

END
