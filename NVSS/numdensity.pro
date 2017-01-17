PRO numdensity,sensitivity
name = 'CATALOG.FIT'
;a=mrdfits(name)
a = mrdfits(name,1)
ra = a.RA_2000_
dec = a.DEC_2000_
intensity = a.PEAK_INT
major = a.MAJOR_AX
minor = a.MINOR_AX
pa = a.POSANGLE
qcenter = a.Q_CENTER
ucenter = a.U_CENTER
flux = a.P_FLUX
irms = a.I_RMS
polrms = a.POL_RMS
resrms = a.RES_RMS
respeak = a.RES_PEAK
resflux = a.RES_FLUX
centerx = a.CENTER_X
centery = a.CENTER_Y
jd = a.JD_PROCESSED

;good=where((intensity gt 0.1) and (intensity lt 5.0) and $
;        (sqrt(major^2+minor^2) lt 0.02) and $
;        (dec gt -15.0) and (dec lt 65.0),count)
;T=10.0
T=10.0
; For temp

;print, ra

good=where((intensity gt sensitivity) and $
        (dec gt 10) and (dec lt 20) $
        and (ra gt 10) and (ra lt 20),count)
; For FAST
;good=where((intensity gt T/(40000.0/1380.0)) and $
;        (sqrt(major^2+minor^2) lt 0.02) and $
;        (dec gt -15.0) and (dec lt 65.0),count)
; For Arecibo
;good=where((intensity gt T/(10000.0/1380.0)) and $
;        (dec gt -1.0) and (dec lt 37.0),count)
print,'num of sources:',count
ra1=ra(good)
dec1=dec(good)
intensity1=intensity(good)
flux1=flux(good)



set_plot,'PS'
;filename='calibsources.eps' ; set the file name of the output ps file
filename='numdensity.eps' ; set the file name of the output ps file
device,file=filename,/ENCAPSULATED,/COLOR, BITS=8;,xsize=xsize,ysize=ysize
map_set,0,0,/AITOFF,/ISOTROPIC,/HORIZON,/GRID

lons=indgen(360/45+1)*45-180
lonnames=[180, 135, 90, 45, 0, 315, 270, 225,180]
map_grid,latdel=20,londel=20,lonnames=lonnames,lons=lons,$
color=0.00*!d.n_colors,charthick=2,glinethick=3,/LABEL,/HORIZON

;openw,lun,'NVSS_GALFA.txt',/get_lun
openw,lun,'NVSS_temp.txt',/get_lun


for j=1L,n_elements(ra1) do begin
xyouts, 360.0-ra1[j-1], dec1[j-1], 'o',charsize=0.5
printf,lun,ra1[j-1],dec1[j-1],intensity1[j-1],intensity1[j-1]*(40000.0/1380.0)
endfor
print,n_elements(ra1)

free_lun,lun

device,/CLOSE


END
;** Structure <8699574>, 18 tags, length=88, data length=88, refs=1:
;   RA_2000_        DOUBLE       0.00039996055
;   DEC_2000_       DOUBLE          -34.119346
;   PEAK_INT        FLOAT        0.00248270
;   MAJOR_AX        FLOAT         0.0144693
;   MINOR_AX        FLOAT         0.0125015
;   POSANGLE        FLOAT          -12.5905
;   Q_CENTER        FLOAT      -0.000443626
;   U_CENTER        FLOAT      -0.000704582
;   P_FLUX          FLOAT       0.000732245
;   I_RMS           FLOAT       0.000480345
;   POL_RMS         FLOAT       0.000300825
;   RES_RMS         FLOAT       0.000347723
;   RES_PEAK        FLOAT       0.000882899
;   RES_FLUX        FLOAT      -0.000565786
;   CENTER_X        FLOAT           511.921
;   CENTER_Y        FLOAT           964.276
;   FIELD           STRING    'C0000M36'
;   JD_PROCESSED    LONG           2450823

