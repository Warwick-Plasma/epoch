; ************************
; *       IsoPlot.       *
; *----------------------*
; *  N J Sircombe, 2001  *
; ************************

; This procedure updates the transformation matrix using Trackball
; which in turn reponds to mouse movements

PRO TrackEx_Event, sEvent
  COMPILE_OPT idl2, hidden
  WIDGET_CONTROL, sEvent.Top, GET_UVALUE=sState, /NO_COPY
  bHaveXform = sState.oTrackball->Update(sEvent, TRANSFORM=TrackXform)

  IF (bHaveXform) THEN BEGIN
    sState.oModel->GetProperty, TRANSFORM=ModelXform
    ; Update model
    sState.oModel->SetProperty, TRANSFORM=ModelXform # TrackXform
    ; Draws in low quality whilst moving
    sState.oWindow->SetProperty, QUALITY=0
    ; Redraw
    sState.oWindow->Draw, sState.oView
  ENDIF ELSE BEGIN
    ; Draws in high quality when stationary
    sState.oWindow->SetProperty, QUALITY=2
    ; Redraw
    sState.oWindow->Draw, sState.oView
  ENDELSE

  WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
END

; --------------------------------------------------------------------------

PRO IsoPlot, volData, level, ViewPt=viewPt, Print=print, Win=win, $
    IsoColor=isocolor, XC=xc, YC=yc, ZC=zc, PrCol=PrCol, XG=xg, YG=yg, ZG=zg
  ; NAME:
  ;   IsoPlot
  ; PURPOSE:
  ;   Draws an isosurface for given data at a given value.
  ; CALLING SEQUENCE:
  ;   IsoPlot, VolData, Level
  ; INPUTS:
  ;   VolData = 3D scalar array containing raw data for isosurface plot.
  ;   Level = Numerical value at which the isosurface is to be drawn
  ; KEYWORD PARAMETERS:
  ;   ViewPt = Three element array which sets the initial viewing angle.
  ;       e.g. ViewPt = [a, b, c]  will rotate the model about the x-axis
  ;       by 'a' degrees, about the y-axis by 'b' degrees and about the z-axis
  ;       by 'c' degrees.
  ;   Print = Setting '/Print' at the command line will output the INITIAL
  ;       state of the model to a printer or .eps file. A dialog box allows
  ;       the user to control the output. Default output file is 'xprinter.eps'
  ;   PrCol = Setting '/prcol' will output the initial state in a lower-quality
  ;       color format to the file 'iso.ps'
  ;   Win = Window size (in pixels), default is 600
  ;   IsoColor = Three element array which controls the colour of the
  ;       isosurface by changing the lighting colour.
  ;       e.g. IsoColor = [r, g, b]  where r, g, b represent the proportion of
  ;       red, green and blue light respectively.
  ;       Values must be in the range 0-255.
  ;   XC, YC, ZC = Three two element arrays which, if given, control contour
  ;       plots in the x, y & z planes. They take the form xc, yc, zc = [a, b]
  ;       where 'a' is the coordinate of the plane to be contoured and 'b' is
  ;       the coordinate of the plane onto which the contour will be drawn.
  ;   XG, YG, ZG = Three arrays containing the coordinates of the grid points.
  ;       The first and last element of each array will be used to calulate
  ;       the correct scaling allong the x, y and z axis.
  ; ADDITIONAL:
  ;   The view-point can be changed dynamically whilst the program is running
  ;   by clicking and dragging the mouse across the window.

  COMPILE_OPT idl2

  ; Find the size of input array - volData
  volsize = SIZE(volData, /DIMENSIONS)
  xn = volsize[0] - 1
  yn = volsize[1] - 1
  zn = volsize[2] - 1

  IF (N_ELEMENTS(XG) NE 0) OR ARG_PRESENT(XG) THEN BEGIN
    coConX = ABS(xg[xn] - xg[0])
  ENDIF ELSE BEGIN
    coConX = 1
  ENDELSE

  IF (N_ELEMENTS(YG) NE 0) OR ARG_PRESENT(YG) THEN BEGIN
    coConY = ABS(yg[yn] - yg[0])
  ENDIF ELSE BEGIN
    coConY = 1
  ENDELSE

  IF (N_ELEMENTS(ZG) NE 0) OR ARG_PRESENT(ZG) THEN BEGIN
    coConZ = ABS(zg[zn] - zg[0])
  ENDIF ELSE BEGIN
    coConZ = 1
  ENDELSE

  maxn = ((xn*coConX) > (yn*coConY) > (zn*coConZ))

  ; If no initial viewpoint is given, defalt values are used
  IF (N_ELEMENTS(ViewPt) NE 0) OR ARG_PRESENT(ViewPt) THEN BEGIN
    xR = viewPt[0]
    yR = viewPt[1]
    zR = viewPt[2]
  ENDIF ELSE BEGIN
    xR = -60     ; default values
    yR = 0
    zR = 30
  ENDELSE

  ; Sets the window size. If no value is given, a defalt value is used.
  IF (N_ELEMENTS(Win) NE 0) OR ARG_PRESENT(Win) THEN BEGIN
    xDim = win
    yDim = win
  ENDIF ELSE BEGIN
    xdim = 600    ; default values
    ydim = 600
  ENDELSE

  ; Sets the lighting color. If no value is given, defaluts to white light
  IF NOT KEYWORD_SET(IsoColor) THEN BEGIN
    isoColor = [255, 255, 255] ; default value
  ENDIF

  SHADE_VOLUME, volData, level, verts, polys, /LOW
  ; IsoSurface, volData, level, verts, polys ;Computes the isosurface

  wBase = WIDGET_BASE()
  ; initailise draw window
  wDraw = WIDGET_DRAW(wBase, XSIZE=xdim, YSIZE=yDim, GRAPHICS_LEVEL=2, $
      /BUTTON_EVENTS, /MOTION_EVENTS, /EXPOSE_EVENTS, RETAIN=0)
  WIDGET_CONTROL, wBase, /REALIZE
  WIDGET_CONTROL, wDraw, GET_VALUE=oWindow

  ; initialise main objects
  oTrackball = OBJ_NEW('Trackball', [xDim/2, yDim/2.], xDim/2.)

  ; Main model view. Relative sizes of the model and the window are controled
  ; by the VIEWPLANE_RECT keyword. Replacing '0.9' and '1.8' with larger
  ; factors will make the model appear smaller with respect to the viewing
  ; window.
  ; e.g. VIEWPLANE_RECT = [-3*maxn, -3*maxn, 6*maxn, 6*maxn]
  ; Zclip and EYE should also be changed to match.
  oView = OBJ_NEW('IDLgrView', ZCLIP=[1.8*maxn, -1.8*maxn], $
      VIEWPLANE_RECT=[-0.9*maxn, -0.9*maxn, 1.8*maxn, 1.8*maxn], $
      EYE=1.8*(maxn)+10)
  oModel  = OBJ_NEW('IDLgrModel', LIGHTING=0)
  oModelIso = OBJ_NEW('IDLgrModel')
  oModelRX = OBJ_NEW('IDLgrModel')
  oModelRY = OBJ_NEW('IDLgrModel')
  oIso = OBJ_NEW('IDLgrPolygon', SHADING=1, COLOR=[125, 125, 125])

  colS = INTARR(3, 22)
  ; generate colors for contour plots
  FOR i = 0, 10 DO BEGIN
    colS[*, i] = [(25*i), 0, 250]
  ENDFOR
  FOR i = 1, 11 DO BEGIN
    colS[*, i+10] = [250, (25*i), 250]
  ENDFOR

  IF  (N_ELEMENTS(XC) NE 0) OR ARG_PRESENT(XC) THEN BEGIN
    oModelXC = OBJ_NEW('IDLgrModel', LIGHTING=0)
    oXContour = OBJ_NEW('IDLgrContour', REFORM(volData[xc[0], *, *]), /FILL, $
        C_COLOR=cols, GEOMZ=xc[1]*coConX, /PLANAR, N_LEVELS=22, $
        YCOORD_CONV=[0, coConZ], XCOORD_CONV=[0, coConY])
    oModelXC->Add, oXContour
    oModelXC->Rotate, [1, 0, 0], 90
    oModelXC->Rotate, [0, 0, 1], 90
    oModel->Add, oModelXC
  ENDIF
  IF  (N_ELEMENTS(YC) NE 0) OR ARG_PRESENT(YC) THEN BEGIN
    oModelYC = OBJ_NEW('IDLgrModel', LIGHTING=0)
    oYContour = OBJ_NEW('IDLgrContour', REFORM(volData[*, yc[0], *]), /FILL, $
        C_COLOR=cols, GEOMZ=yc[1]*coConY, /PLANAR, N_LEVELS=22, $
        YCOORD_CONV=[0, coConZ], XCOORD_CONV=[0, coConX])
    oModelYC->Add, oYContour
    oModelYC->Rotate, [1, 0, 0], 90
    oModelYC->Translate, 0, 2*yc[1]*coConY, 0
    oModel->Add, oModelYC
  ENDIF
  IF  (N_ELEMENTS(ZC) NE 0) OR ARG_PRESENT(ZC) THEN BEGIN
    oModelZC = OBJ_NEW('IDLgrModel', LIGHTING=0)
    oZContour = OBJ_NEW('IDLgrContour', REFORM(volData[*, *, zc[0]]), /FILL, $
        C_COLOR=cols, GEOMZ=zc[1]*coConZ, /PLANAR, N_LEVELS=22, $
        YCOORD_CONV=[0, coConY], XCOORD_CONV=[0, coConX])
    oModelZC->Add, oZContour
    oModel->Add, oModelZC
  ENDIF

  ; initialise light sources
  oLight = OBJ_NEW('IDLgrLight', $
      LOCATION=[0, 0, 0], TYPE=1, COLOR=isoColor)
  oLight2 = OBJ_NEW('IDLgrLight', $
      LOCATION=[xn*coConX, yn*coConY, zn*coConZ], TYPE=1, COLOR=isoColor)
  oLight3 = OBJ_NEW('IDLgrLight', $
      LOCATION=[0, 0, 0], TYPE=0, INTENSITY=0.3)

  ; Initialise the axis objects used to make up the bounding box
  oXaxis = OBJ_NEW('IDLgrAxis', 0)
  oYaxis = OBJ_NEW('IDLgrAxis', 1)
  oZaxis = OBJ_NEW('IDLgrAxis', 2)
  oXaxis2 = OBJ_NEW('IDLgrAxis', 0)
  oYaxis2 = OBJ_NEW('IDLgrAxis', 1)
  oZaxis2 = OBJ_NEW('IDLgrAxis', 2)
  oXaxis3 = OBJ_NEW('IDLgrAxis', 0)
  oYaxis3 = OBJ_NEW('IDLgrAxis', 1)
  oZaxis3 = OBJ_NEW('IDLgrAxis', 2)
  oXaxis4 = OBJ_NEW('IDLgrAxis', 0)
  oYaxis4 = OBJ_NEW('IDLgrAxis', 1)
  oZaxis4 = OBJ_NEW('IDLgrAxis', 2)

  ; Pass isosurface data over to IDLgrSurface which biulds the 3d model
  oIso->SetProperty, DATA=verts, POLYGONS=polys, $
      XCOORD_CONV=[0, coConX], YCOORD_CONV=[0, coConY], ZCOORD_CONV=[0, coConZ]

  ; Setup axis, removing labels and ticks
  oXaxis ->SetProperty, TICKLEN=0, LOCATION=[0, 0, 0], $
      /NOTEXT, /EXACT, RANGE=[0, xn], XCOORD_CONV=[0, coConX]
  oYaxis ->SetProperty, TICKLEN=0, LOCATION=[0, 0, 0], $
      /NOTEXT, /EXACT, RANGE=[0, yn], YCOORD_CONV=[0, coConY]
  oZaxis ->SetProperty, TICKLEN=0, LOCATION=[0, 0, 0], $
      /NOTEXT, /EXACT, RANGE=[0, zn], ZCOORD_CONV=[0, coConZ]
  oXaxis2->SetProperty, TICKLEN=0, LOCATION=[0, yn*coConY, 0], $
      /NOTEXT, /EXACT, RANGE=[0, xn], XCOORD_CONV=[0, coConX]
  oYaxis2->SetProperty, TICKLEN=0, LOCATION=[xn*coConX, 0, 0], $
      /NOTEXT, /EXACT, RANGE=[0, yn], YCOORD_CONV=[0, coConY]
  oZaxis2->SetProperty, TICKLEN=0, LOCATION=[xn*coConX, 0, 0], $
      /NOTEXT, /EXACT, RANGE=[0, zn], ZCOORD_CONV=[0, coConZ]
  oXaxis3->SetProperty, TICKLEN=0, LOCATION=[0, yn*coConY, zn*coConZ], $
      /NOTEXT, /EXACT, RANGE=[0, xn], XCOORD_CONV=[0, coConX]
  oYaxis3->SetProperty, TICKLEN=0, LOCATION=[xn*coConX, 0, zn*coConZ], $
      /NOTEXT, /EXACT, RANGE=[0, yn], YCOORD_CONV=[0, coConY]
  oZaxis3->SetProperty, TICKLEN=0, LOCATION=[xn*coConX, yn*coConY, 0], $
      /NOTEXT, /EXACT, RANGE=[0, zn], ZCOORD_CONV=[0, coConZ]
  oXaxis4->SetProperty, TICKLEN=0, LOCATION=[0, 0, zn*coConZ], $
      /NOTEXT, /EXACT, RANGE=[0, xn], XCOORD_CONV=[0, coConX]
  oYaxis4->SetProperty, TICKLEN=0, LOCATION=[0, 0, zn*coConZ], $
      /NOTEXT, /EXACT, RANGE=[0, yn], YCOORD_CONV=[0, coConY]
  oZaxis4->SetProperty, TICKLEN=0, LOCATION=[0, yn*coConY, 0], $
      /NOTEXT, /EXACT, RANGE=[0, zn], ZCOORD_CONV=[0, coConZ]

  ; Add axis, lighting, contours and isosurface to the model
  oView->Add, oModelRX
  oModel->Add, oModelIso
  oModelIso->Add, oIso
  ; RX and RY models are used to setup the initial viewpt relative to
  ; absolute x, y, z axis and not the models own axis (which will always
  ; rotate with the models.
  oModelRY->Add, oModel
  oModelRX->Add, oModelRY
  oModel->Add, oXaxis
  oModel->Add, oYaxis
  oModel->Add, oZAxis
  oModel->Add, oXaxis2
  oModel->Add, oYaxis2
  oModel->Add, oZAxis2
  oModel->Add, oXaxis3
  oModel->Add, oYaxis3
  oModel->Add, oZAxis3
  oModel->Add, oXaxis4
  oModel->Add, oYaxis4
  oModel->Add, oZAxis4
  oModel->Add, oLight
  oModel->Add, oLight2
  oModel->Add, oLight3
  ; Setup the initial viewpoint
  oModel->Translate, -((xn*coConX)/2.0), -((yn*coConY)/2.0), -((zn*coConZ)/2.0)
  oModel->Rotate, [0, 0, 1], zR
  oModelRY->Rotate, [0, 1, 0], yR
  oModelRX->Rotate, [1, 0, 0], xR

  XMANAGER, 'TrackEx', wBase, /NO_BLOCK
  oWindow->Draw, oView    ; Draw isosurface on screen

  ; If /print keyword is set then send inital viewpoint to the printer/eps file
  IF KEYWORD_SET(print) THEN BEGIN
    oPrint = OBJ_NEW('IDLgrPrinter')
    ; prnOK1 = DIALOG_PRINTERSETUP(oPrint) ; activate for more printer control
    ; IF prnOK1 EQ 1 THEN BEGIN  ; activate for more printer control
    ; Best results with compression on and scale = 3
    prnOK2 = DIALOG_PRINTJOB(oPrint)
    IF prnOK2 EQ 1 THEN BEGIN
      oPrint->Draw, oView, VECTOR=0
      oPrint->NewDocument
    ENDIF
    ; ENDIF ;activate for more printer control
  ENDIF
  IF KEYWORD_SET(prCol) THEN BEGIN ; color print output
    ; print using clipboard object
    oPrintC = OBJ_NEW('IDLgrClipboard', DIMENSIONS=[xDim, yDim])
    oPrintC->Draw, oView, VECTOR=0, FILENAME='iso.ps', POSTSCRIPT=1
  ENDIF

  sState = {oView:oView, oTrackBall:oTrackBall, oModel:oModel, oWindow:oWindow}

  WIDGET_CONTROL, wBase, SET_UVALUE=sState, /NO_COPY
END
