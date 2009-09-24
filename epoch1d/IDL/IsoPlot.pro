;************************
;*       IsoPlot.       *
;*----------------------*
;*  N J Sircombe, 2001  *
;************************


;This procedure updates the tranceformation matrix using Trackball
;which in turn reponds to mouse movements
pro TrackEx_Event, sEvent 
    WIDGET_CONTROL, sEvent.Top, GET_UVALUE=sState, /NO_COPY
    bHaveXform = sState.oTrackball->Update(sEvent,TRANSFORM=TrackXform)
    IF (bHaveXform) THEN BEGIN
    	sState.oModel->GetProperty,TRANSFORM=ModelXform
    	sState.oModel->SetProperty, TRANSFORM =ModelXform # TrackXform	;Update model
    	sState.oWindow->SetProperty, QUALITY = 0 ;Draws in low quality whilst moving
	sState.oWindow ->Draw,sState.oView  ;Redraw
    ENDIF ELSE BEGIN
    	sState.oWindow->SetProperty, QUALITY = 2 ;Draws in high quality when stationary
	sState.oWindow ->Draw,sState.oView  ;Redraw
    ENDELSE
    
    WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
END

;*******************************************************************************

PRO IsoPlot, volData, level,ViewPt = viewPt, Print=print, Win=win, IsoColor=isocolor, $
    	XC=xc, YC=yc, ZC=zc, PrCol=PrCol, XG=xg,YG=yg, ZG=zg
; NAME:
;   	IsoPlot
; PURPOSE:
;	Draws an isosurface for given data at a given value.
; CALLING SEQUENCE:
;	IsoPlot, VolData, Level
; INPUTS:
;	VolData = 3D scalar array containing raw data for isosurface plot.
;   	Level = Numerical value at which the isosurface is to be drawn
; KEYWORD PARAMETERS:
;	ViewPt = Three element array which sets the initial viewing angle.
;   	    e.g. ViewPt = [a,b,c]  will rotate the model about the x-axis
;   	    by 'a' degrees, about the y-axis by 'b' degrees and about the z-axis
;   	    by 'c' degrees.
;   	Print = Setting '/Print' at the command line will output the INITIAL state
;   	    of the model to a printer or .eps file. A dialog box allows the user
;   	    to control the output. Default output file is 'xprinter.eps'.
;   	PrCol = Setting '/prcol' will output the initial state in a lower-quality
;   	    color format to the file 'iso.ps'
;   	Win = Window size (in pixels), default is 600
;   	IsoColor = Three element array which controls the colour of the isosurface
;   	    by changing the lighting colour.
;   	    e.g. IsoColor = [r,g,b]  where r,g,b represent the proportion of
;   	    red, green and blue light respectively. 
;   	    Values must be in the range 0-255.
;   	XC,YC,ZC = Three two element arrays which, if given, control contour plots
;   	    in the x,y & z planes. They take the form xc, yc, zc=[a,b] where 'a' is the
;   	    coordinate of the plane to be contoured and 'b' is the coordinate of the
;   	    plane onto which the contour will be drawn.
;   	XG,YG,ZG = Three arrays containing the coordinates of the grid points. The first
;   	    and last element of each array will be used to calulate the correct scaling
;   	    allong the x,y and z axis.
; ADDITIONAL:
;	The view-point can be changed dynamically whilst the program is running
;   	    by clicking and dragging the mouse across the window.


;Find the size of input array - volData
    volsize=SIZE(volData, /DIMENSIONS) 
    xn=volsize[0] - 1
    yn=volsize[1] - 1
    zn=volsize[2] - 1

    IF  (N_Elements(XG) NE 0) OR Arg_Present(XG) THEN BEGIN 
    	coConX = (ABS(xg[xn]-xg[0]))
    ENDIF ELSE BEGIN
	coConX = 1
    ENDELSE
    IF  (N_Elements(YG) NE 0) OR Arg_Present(YG) THEN BEGIN 
    	coConY = (ABS(yg[yn]-yg[0]))
    ENDIF ELSE BEGIN
	coConY = 1
    ENDELSE
    IF  (N_Elements(ZG) NE 0) OR Arg_Present(ZG) THEN BEGIN 
    	coConZ = (ABS(zg[zn]-zg[0]))
    ENDIF ELSE BEGIN
	coConZ = 1
    ENDELSE
    
    maxn=((xn*coConX) > (yn*coConY) > (zn*coConZ))

;If no initial viewpoint is given, defalt values are used
    IF  (N_Elements(ViewPt) NE 0) OR Arg_Present(ViewPt) THEN BEGIN 
   	xR = viewPt[0]
     	yR = viewPt[1]
    	zR = viewPt[2]
    ENDIF ELSE BEGIN
	xR = -60     ;default values
	yR = 0
	zR = 30
    ENDELSE
;Sets the window size. If no value is given, a defalt value is used.
    IF  (N_Elements(Win) NE 0) OR Arg_Present(Win) THEN BEGIN 
    	xDim=win
	yDim=win
    ENDIF ELSE BEGIN
	xdim=600    ;default values
	ydim=600
    ENDELSE
;Sets the lighting color. If no value is given, defaluts to white light
    IF  NOT Keyword_Set(IsoColor) THEN BEGIN 
    	isoColor = [255,255,255] ;default value
    ENDIF
    shade_volume,volData,level,verts,polys,/low
    ;IsoSurface, volData, level, verts, polys	;Computes the isosurface
    
    wBase = WIDGET_BASE()
;initailise draw window    
    wDraw = WIDGET_DRAW(wBase,XSIZE = xdim,YSIZE = yDim, GRAPHICS_LEVEL=2,$
	    /BUTTON_EVENTS, /MOTION_EVENTS,/EXPOSE_EVENTS,RETAIN=0)
    WIDGET_CONTROL, wBase,/REALIZE
    WIDGET_CONTROL,wDraw,GET_VALUE=oWindow   

;initialise main objects
    oTrackball = Obj_New('Trackball', [xDim/2,yDim/2.], xDim/2.)

;Main model view. Relative sizes of the model and the window are controled by the
;VIEWPLANE_RECT keyword. Replacing '0.9' and '1.8' with larger factors will make the model
;appear smaller with respect to the viewing window.
;e.g. VIEWPLANE_RECT=[-3*maxn,-3*maxn,6*maxn,6*maxn]
;Zclip and EYE should also be changed to match. 
    oView = Obj_New('IDLgrView', ZClip = [1.8*maxn,-1.8*maxn],$
    	VIEWPLANE_RECT = [-0.9*maxn,-0.9*maxn, 1.8*maxn,1.8*maxn], EYE = 1.8*(maxn)+10)
    oModel  = Obj_New('IDLgrModel', lighting=0)
    oModelIso = Obj_New('IDLgrModel')
    oModelRX = Obj_New('IDLgrModel')
    oModelRY = Obj_New('IDLgrModel')
    oIso = Obj_New('IDLgrPolygon', Shading = 1, COLOR = [125,125,125]) 

    colS=IntArr(3,22)
;generate colors for contour plots
     FOR i=0,10 DO BEGIN
     	colS[*,i]=[(25*i),0,250]
     ENDFOR
     FOR i=1,11 DO BEGIN
     	colS[*,i+10]=[250,(25*i),250]
     ENDFOR

    IF  (N_Elements(XC) NE 0) OR Arg_Present(XC) THEN BEGIN
    	oModelXC=Obj_New('IDLgrModel', lighting=0)
       	oXContour=Obj_New('IDLgrContour', REFORM(volData[xc[0],*,*]), /Fill,C_COLOR=cols, $
    	    	GEOMZ=xc[1]*coConX,/Planar,N_LEVELS=22, YCOORD_CONV=[0,coConZ],XCOORD_CONV=[0,coConY])
    	oModelXC->Add, oXContour
    	oModelXC->Rotate, [1,0,0], 90
	oModelXC->Rotate, [0,0,1], 90
	oModel->Add, oModelXC
    ENDIF
    IF  (N_Elements(YC) NE 0) OR Arg_Present(YC) THEN BEGIN
    	oModelYC=Obj_New('IDLgrModel', lighting=0)
       	oYContour=Obj_New('IDLgrContour', REFORM(volData[*,yc[0],*]), /Fill,C_COLOR=cols, $
    	    	GEOMZ=yc[1]*coConY,/Planar, N_LEVELS=22, YCOORD_CONV=[0,coConZ],XCOORD_CONV=[0,coConX])
    	oModelYC->Add, oYContour
    	oModelYC->Rotate, [1,0,0], 90
 	oModelYC->Translate, 0,2*yc[1]*coConY,0
    	oModel->Add, oModelYC
    ENDIF
    IF  (N_Elements(ZC) NE 0) OR Arg_Present(ZC) THEN BEGIN
    	oModelZC=Obj_New('IDLgrModel', lighting=0)
       	oZContour=Obj_New('IDLgrContour', REFORM(volData[*,*,zc[0]]), /Fill,C_COLOR=cols, $
    	    	GEOMZ=zc[1]*coConZ,/Planar,N_LEVELS=22, YCOORD_CONV=[0,coConY],XCOORD_CONV=[0,coConX])
    	oModelZC->Add, oZContour
    	oModel->Add, oModelZC
    ENDIF

;initialise light sources
    oLight = OBJ_NEW('IDLgrLight', LOCATION=[0,0,0], TYPE=1, COLOR=isoColor)
    oLight2 = OBJ_NEW('IDLgrLight', LOCATION=[xn*coConX,yn*coConY,zn*coConZ], TYPE=1, COLOR=isoColor)
    oLight3= Obj_New('IDLgrLight', LOCATION=[0,0,0], TYPE=0, INTENSITY = 0.3)
;Initialise the axis objects used to make up the bounding box
    oXaxis = Obj_New('IDLgrAxis',0)
    oYaxis = Obj_New('IDLgrAxis',1)
    oZaxis = Obj_New('IDLgrAxis',2)
    oXaxis2 = Obj_New('IDLgrAxis',0)
    oYaxis2 = Obj_New('IDLgrAxis',1)
    oZaxis2 = Obj_New('IDLgrAxis',2)
    oXaxis3 = Obj_New('IDLgrAxis',0)
    oYaxis3 = Obj_New('IDLgrAxis',1)
    oZaxis3 = Obj_New('IDLgrAxis',2)
    oXaxis4 = Obj_New('IDLgrAxis',0)
    oYaxis4 = Obj_New('IDLgrAxis',1)
    oZaxis4 = Obj_New('IDLgrAxis',2)
;Pass isosurface data over to IDLgrSurface which biulds the 3d model    
    oIso->SetProperty, data=verts,Polygons = polys,XCOORD_CONV=[0,coConX],YCOORD_CONV=[0,coConY],ZCOORD_CONV=[0,coConZ]
;Setup axis, removing labels and ticks    
    oXaxis ->SetProperty, TickLen = 0, Location = [0,0,0], /NoText,/Exact, Range =  [0,xn], XCOORD_CONV=[0,coConX]
    oYaxis ->SetProperty, TickLen = 0, Location = [0,0,0], /NoText,/Exact, Range =  [0,yn], YCOORD_CONV=[0,coConY]
    oZaxis ->SetProperty, TickLen = 0, Location = [0,0,0], /NoText,/Exact, Range =  [0,zn], ZCOORD_CONV=[0,coConZ]
    oXaxis2->SetProperty, TickLen = 0, Location = [0,yn*coConY,0], /NoText,/Exact, Range =  [0,xn], XCOORD_CONV=[0,coConX]
    oYaxis2->SetProperty, TickLen = 0, Location = [xn*coConX,0,0], /NoText,/Exact, Range =  [0,yn], YCOORD_CONV=[0,coConY]
    oZaxis2->SetProperty, TickLen = 0, Location = [xn*coConX,0,0], /NoText,/Exact, Range =  [0,zn], ZCOORD_CONV=[0,coConZ]
    oXaxis3->SetProperty, TickLen = 0, Location = [0,yn*coConY,zn*coConZ],/NoText,/Exact, Range =  [0,xn], XCOORD_CONV=[0,coConX]
    oYaxis3->SetProperty, TickLen = 0, Location = [xn*coConX,0,zn*coConZ],/NoText,/Exact, Range =  [0,yn], YCOORD_CONV=[0,coConY]
    oZaxis3->SetProperty, TickLen = 0, Location = [xn*coConX,yn*coConY,0],/NoText,/Exact, Range =  [0,zn], ZCOORD_CONV=[0,coConZ]
    oXaxis4->SetProperty, TickLen = 0, Location = [0,0,zn*coConZ], /NoText,/Exact, Range =  [0,xn], XCOORD_CONV=[0,coConX]
    oYaxis4->SetProperty, TickLen = 0, Location = [0,0,zn*coConZ], /NoText,/Exact, Range =  [0,yn], YCOORD_CONV=[0,coConY]
    oZaxis4->SetProperty, TickLen = 0, Location = [0,yn*coConY,0], /NoText,/Exact, Range =  [0,zn], ZCOORD_CONV=[0,coConZ]
;Add axis, lighting, contours and isosurface to the model
    oView->Add, oModelRX
    oModel->Add,oModelIso
    oModelIso->Add,oIso
    oModelRY->Add,oModel    ;RX and RY models are used to setup the initial viewpt relative
    oModelRX->Add,oModelRY  ;to absolute x,y,z axis and not the models own axis (which will
    oModel->Add,oXaxis	    ;always rotate with the models.
    oModel->Add,oYaxis
    oModel->Add,oZAxis
    oModel->Add,oXaxis2
    oModel->Add,oYaxis2
    oModel->Add,oZAxis2
    oModel->Add,oXaxis3
    oModel->Add,oYaxis3
    oModel->Add,oZAxis3
    oModel->Add,oXaxis4
    oModel->Add,oYaxis4
    oModel->Add,oZAxis4
    oModel->Add,oLight
    oModel->Add,oLight2
    oModel->Add,oLight3
;Setup the initial viewpoint
    oModel->Translate, -((xn*coConX)/2.0),-((yn*coConY)/2.0),-((zn*coConZ)/2.0)
    oModel->Rotate, [0,0,1], zR
    oModelRY->Rotate, [0,1,0], yR
    oModelRX->Rotate, [1,0,0], xR


    XMANAGER, 'TrackEx', wBase,/NO_BLOCK
    oWindow->Draw, oView    ;Draw isosurface on screen

;If /print keyword is set then send inital viewpoint to the printer/eps file
    IF Keyword_Set(print) THEN BEGIN	
 	oPrint=Obj_New('IDLgrPrinter')
	;prnOK1=DIALOG_PRINTERSETUP(oPrint)   	;activate for more printer control
	;IF prnOK1 EQ 1 THEN BEGIN		;   "      "   "      "       "
	    prnOK2=DIALOG_PRINTJOB(oPrint)	     ;Best results with compression on and scale=3
	    IF prnOK2 EQ 1 THEN BEGIN
	    	oPrint->Draw,oView, Vector=0
    	    	oPrint->NewDocument	
	    ENDIF
	;ENDIF	;activate for more printer control
    ENDIF
    IF Keyword_Set(prCol) THEN BEGIN	;color print output
	oPrintC=Obj_New('IDLgrClipboard',  DIMENSIONS=[xDim,yDim]) ;print using clipboard object 
	oPrintC->Draw,oView, Vector=0,FILENAME='iso.ps',POSTSCRIPT=1
    ENDIF


    sState = {oView:oView, $
	oTrackBall:oTrackBall, $
	oModel:oModel,$
        oWindow:oWindow}  

    WIDGET_CONTROL, wBase, SET_UVALUE=sState,/NO_COPY
END
