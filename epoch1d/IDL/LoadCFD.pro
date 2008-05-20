FUNCTION LoadCFDFile, filename, Variables=requestv,request_classes=requestc, _extra=extra

COMMON BlockTypes,TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE
COMMON MeshTypes, MESH_CARTESIAN, MESH_PARTICLE

on_error, 2

Version=1L
Revision=0L

base_mesh_in_place=0

offset=0LL

IF N_PARAMS() EQ 0 THEN BEGIN
    print , "Usage: result = LoadCFDFile(Filename[,/variables])"
    RETURN, "Usage: result = LoadCFDFile(Filename[,/variables])"
ENDIF

; array of names of parameters
IF (N_ELEMENTS(extra) NE 0) THEN BEGIN
    NAME_ARR=TAG_NAMES(extra)
    ELEMENT_BLOCK=intarr(N_ELEMENTS(NAME_ARR))
    ELEMENT_BLOCK(*)=0
ENDIF ELSE BEGIN
    NAME_ARR=""
    ELEMENT_BLOCK=1
ENDELSE

header = assoc(1,{cfd:bytarr(3), headeroff:0L, blockheaderoff:0L, Version:0L, Revision:0L, MaxString:0L, nblocks:0L},0,/packed)
close,1
openr,1,filename
fileheader = header[0]
MaxStringLen=fileheader.MaxString
offset=fileheader.headeroff

;Whole load of boring tests

IF (STRING(fileheader.cfd) NE "CFD") THEN BEGIN
    PRINT,"The file ",filename," is not a valid CFD file"
    CLOSE,1
    RETURN,0
ENDIF

IF (fileheader.Version GT Version) THEN BEGIN
    PRINT,"The file ",filename," is of a version too high to be read by this program"
    PRINT,"Please contact the CFSA, University of Warwick to obtain a new version"
    CLOSE,1
    RETURN,0
ENDIF

IF (fileheader.Revision GT Revision) THEN BEGIN
    PRINT,"WARNING : The file ",filename," has a higher revision number than this reader"
    PRINT,"Not all data in the file will be available"
    PRINT,"Please contact the CFSA, University of Warwick to obtain a new version"
ENDIF

IF (fileheader.nblocks LE 0) THEN BEGIN
    PRINT,"The file ",filename," either contains no blocks or is corrupted"
    CLOSE,1
    RETURN,0
ENDIF

;The file seems valid, spool through blocks

f = {filename: filename}

IF (N_ELEMENTS(requestv) NE 0) THEN BEGIN
    PRINT,"Available elements are "
ENDIF
vBlock=0
FOR iBlock=0,fileheader.nblocks-1 DO BEGIN
    data=assoc(1,{Name:bytarr(MaxStringLen),Class:bytarr(MaxStringLen),Type:0L,BlockMDLen:0LL,BlockLen:0LL},offset,/packed)
    blockheader=data[0]
                                ;Read what we know of the block header, so skip the rest
    offset=offset+fileheader.blockheaderoff
    IF (N_ELEMENTS(requestv) NE 0) THEN BEGIN
        HandleBlock,fileheader,blockheader,f,offset,name_arr,/onlymd,md=md
        q=ReturnIDLUsable(blockheader,md)
        IF (q EQ 1) THEN BEGIN
          PRINT,STRTRIM(STRING(vBlock+1),2),") ",STRTRIM(STRING(blockheader.Name),2)," : " + ReturnFriendlyTypeName(blockheader,md)
          vBlock=vBlock+1
        ENDIF
	ELEMENT_BLOCK(*)=1
    ENDIF ELSE BEGIN
        HandleBlock,fileheader,blockheader,f,offset,name_arr,element_block
    ENDELSE
    offset=offset+blockheader.BlockLen
ENDFOR

;Check for any names given which have not been understood
Errcount=0
FOR iEl=0,N_ELEMENTS(name_arr)-1 DO BEGIN
	IF (ELEMENT_BLOCK(iEl) EQ 0) THEN BEGIN
		PRINT,"WARNING! Unrecognised variable requested (",name_arr(iEl),")"
                Errcount=Errcount+1
	ENDIF
ENDFOR
IF(Errcount NE 0) THEN PRINT,"You have specified nonexistant variables. To list available variables, use the '/variables' switch"
close,1

RETURN, f

END

PRO handleblock,fileheader,blockheader,outputobject,offset,name_arr,element_block,md=md,onlymd=onlymd
COMMON BlockTypes,TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE

NameMatch=CheckName(blockheader,Name_Arr,element_block)
IF (NameMatch EQ 1 || blockheader.Type EQ TYPE_MESH) THEN BEGIN
    IF (blockheader.Type EQ TYPE_MESH) THEN BEGIN
        GetMesh,fileheader,blockheader,outputobject,offset,onlymd=onlymd,md=md,byname=NameMatch
        RETURN
    ENDIF
    IF (blockheader.Type EQ TYPE_MESH_VARIABLE) THEN  BEGIN
        GetMeshVar,fileheader,blockheader,outputobject,offset,onlymd=onlymd,md=md
        RETURN
    ENDIF
ENDIF
END

PRO GetMesh,file_header,block_header,output_struct,offset,onlymd=onlymd,md=md,byname=byname
COMMON MeshTypes, MESH_CARTESIAN, MESH_PARTICLE
COMMON ParticleCoords, PARTICLE_CARTESIAN

mesh=assoc(1,{MeshType:0L,nd:0L,sof:0L},offset,/packed)
mesh_header=mesh[0]

mdonly_f=0
IF (n_elements(onlymd) NE 0) THEN mdonly_f=1
byname_f=0
IF (n_elements(byname) NE 0) THEN byname_f=byname

IF (mesh_header.MeshType EQ MESH_CARTESIAN) THEN BEGIN
                                ;Read in the cartesian metadata
;Read in the actual mesh

    names=["x","y","z","a","b","c","d","e","f"]

        mesh=assoc(1,{MeshType:0L,nd:0L,sof:0L,npts:lonarr(mesh_header.nd)},offset,/packed)
        mesh_header=mesh[0]
        IF (mdonly_f NE 1) THEN BEGIN
            IF (mesh_header.sof EQ 4) THEN BEGIN
               datastruct=CREATE_STRUCT(names[0],fltarr(mesh_header.npts[0]))
               FOR iDim=1,mesh_header.nd-1 DO BEGIN
                  datastruct=CREATE_STRUCT(datastruct,names[iDim],fltarr(mesh_header.npts[iDim]))
               END
            ENDIF
            IF (mesh_header.sof EQ 8) THEN BEGIN
               datastruct=CREATE_STRUCT(names[0],dblarr(mesh_header.npts[0]))	
               FOR iDim=1,mesh_header.nd-1 DO BEGIN
                 datastruct=CREATE_STRUCT(datastruct,names[iDim],dblarr(mesh_header.npts[iDim]))
               END
           ENDIF
        ENDIF

    IF(mdonly_f NE 1) THEN BEGIN
        data=assoc(1,datastruct,offset+block_header.BlockMDLen,/packed)
        d=data[0]
	d=CREATE_STRUCT(d,mesh_header)
    ENDIF
ENDIF ELSE IF (mesh_header.MeshType EQ MESH_PARTICLE AND byname_f EQ 1) THEN BEGIN
    mesh=assoc(1,{MeshType:0L,nd:0L,sof:0L,CoordType:0L,npart:0LL},offset,/packed)
    mesh_header=mesh[0]
    IF (mdonly_f NE 1) THEN BEGIN
;        IF (mesh_header.CoordType EQ PARTICLE_CARTESIAN) THEN BEGIN
            IF (mesh_header.sof EQ 4) THEN datastruct=CREATE_STRUCT("ParticlePositions",fltarr(mesh_header.npart,mesh_header.nd))
            IF (mesh_header.sof EQ 8) THEN datastruct=CREATE_STRUCT("ParticlePositions",dblarr(mesh_header.npart,mesh_header.nd))
            PRINT,"Warning, you have loaded a particle mesh. At present, IDL support for particle meshes is very limited"
            PRINT,"The raw data has simply been loaded, and you will have to decode it yourself"
            data=assoc(1,datastruct,offset+block_header.BlockMDLen,/packed)
            d=data[0]
	    d=CREATE_STRUCT("MetaData",mesh_header,d)
;        ENDIF
    ENDIF
ENDIF
md=mesh_header
IF(mdonly_f NE 1 AND N_ELEMENTS(d) NE 0) THEN BEGIN
    output_struct=CREATE_STRUCT(output_struct,STRTRIM(STRING(block_header.Name),2),d)
ENDIF
END 

PRO GetMeshVar,file_header,block_header,output_struct,offset,md=md,onlymd=onlymd
COMMON VarTypes, VAR_CARTESIAN, VAR_PARTICLE

var=assoc(1,{VarType:0L,nd:0L,sof:0L},offset,/packed)
var_header=var[0]
mdonly_f=1
IF (N_ELEMENTS(onlymd) EQ 0) THEN mdonly_f=0
IF (var_header.VarType EQ VAR_CARTESIAN) THEN BEGIN
;Read in the actual variable
        var=assoc(1,{VarType:0L,nd:0L,sof:0L,npts:lonarr(var_header.nd)},offset,/packed)
        var_header=var[0]
        IF (mdonly_f NE 1) THEN BEGIN
            IF (var_header.sof EQ 4) THEN datastruct=CREATE_STRUCT(STRTRIM(STRING(block_header.name),2),fltarr(var_header.npts))
            IF (var_header.sof EQ 8) THEN datastruct=CREATE_STRUCT(STRTRIM(STRING(block_header.name),2),dblarr(var_header.npts))
        ENDIF

    md=var_header
    IF (mdonly_f NE 1) THEN BEGIN
        data=assoc(1,datastruct,offset+block_header.BlockMDLen,/packed)
        d=data[0]
    ENDIF
ENDIF ELSE IF (var_header.VarType EQ VAR_PARTICLE) THEN BEGIN
	var=assoc(1,{VarType:0L,nd:0L,sof:0L,npart:0LL},offset,/packed)
	var_header=var[0]
        IF (mdonly_f NE 1) THEN BEGIN
            IF (var_header.sof EQ 4) THEN datastruct=CREATE_STRUCT(STRTRIM(STRING(block_header.name),2),fltarr(var_header.npart))
            IF (var_header.sof EQ 8) THEN datastruct=CREATE_STRUCT(STRTRIM(STRING(block_header.name),2),dblarr(var_header.npart))
        ENDIF
	md=var_header
	IF (mdonly_f NE 1) THEN BEGIN
        	data=assoc(1,datastruct,offset+block_header.BlockMDLen,/packed)
	        d=data[0]
	ENDIF

ENDIF
IF(mdonly_f NE 1 AND N_ELEMENTS(d) NE 0) THEN BEGIN
    output_struct=CREATE_STRUCT(output_struct,d)
ENDIF
END 
