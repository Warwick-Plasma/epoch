COMMON background, wkdir_global, retro_global
COMMON gdlset, gdl
DEFSYSV, '!GDL', EXISTS=gdl

retro_global = 0
wkdir_global = "Data"
@ IDL/StartCFD.pro

.r IDL/widget
.r IDL/StartPIC.pro

q0 = 1.602176565d-19 ; C
m0 = 9.10938291d-31  ; kg
v0 = 2.99792458d8    ; m/s^2
kb = 1.3806488d-23   ; J/K
mu0 = 4.0d-7 * !dpi  ; N/A^2
epsilon0 = 8.8541878176203899d-12 ; F/m
h_planck = 6.62606957d-34 ; J s

device, true_color=24
device, decompose=0
device, retain=2

!p.charsize = 2
loadct, 3
