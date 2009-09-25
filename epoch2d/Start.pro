COMMON background,wkdir_global
wkdir_global="Data"
@ IDL/StartCFD.pro
.r IDL/StartPIC.pro

Q0 = 1.60217646d-19 ;C
M0 = 9.10938188d-31 ;kg
V0 = 2.98d8         ;ms^(-2)
kb = 1.3806503d-23  ;m^2kgs(-2)K^(-1)
epsilon0 = 8.85418782d-12
h_planck=6.626068e-34

device,true_color=24
device,decompose=0
device,retain=2

!p.charsize=2
loadct,3
