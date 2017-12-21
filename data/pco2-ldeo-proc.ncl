; In December 2010 our colleague, R. Wanninkhof pointed out 
; a problem with the flux data for the month of December.  
; The file of ice coverage for December was corrupted and
; showed zero ice for the entire month, worldwide.  This 
; has been corrected with the estimated percent of ice 
; and the flux recalculated.  Version "c" of the data files 
; contain this correction.
; 
; In October 2009 our colleague, N. Metzl pointed out a problem
; with one of the data sets used in the interpolation.  Specifically
; the OISO cruises in the Indian Ocean.  We corrected these data
; and re-ran the interpolation and the files on the current web
; site reflect those corrections.
; 
; In June 2009, we found a minor programming error in the
; interpolation algorithm used in Takahashi et al., Deep-Sea 
; Research II, vol. 56, 554-577, 2009. The error affects only 
; the CO2 and SST values in the box areas located at the prime 
; meridian.  This file contains revised values for the climatological
; mean monthly values in 4 degree latitude by 5 degree longitude 
; pixel areas.
; 
; LAT is the latitude of the center of the designated box.  
;   Unit is decimal degrees.
; 
; LON is the longitude of the center of the designated box.
;   Unit is decimal degrees all in EAST longitude.
; 
; MONTH is the calendar month, January = 1, February = 2, etc.
; 
; PCO2_SW is monthly mean pCO2 (uatm) in surface water at SST.
;   The individual pCO2 values measured in different year are
;   normalized to a reference year 2000 using 1.5 uatm per year
;   as the rate of pCO2 change.  The grid 4 degree latitude by
;   5 degree longitude, centered on 2.5, 7.5, etc longitude and
;   0, 4, 8, etc. latitude.
;  
; MOD_SST is based on the individual SST observations that were
;   made concurrently with each pCO2 measurement, and interpolated
;   using  the same interpolation scheme run for the pCO2 data.
;   El Nino months for the area from 10N - 10S and 170W to 100W
;   are excluded from the database before the interpolation is made. 
; 
; VCO2_AIR is the concentration of CO2 in dry air from the GlobalView
;   CO2 database for the week nearest the center date of each month,
;   in 2000.  Unit is parts per million.
; 
; PCO2_AIR is monthly pCO2 (uatm) in air for 2000.  It is computed using:
;   (a) the zonal mean CO2 conc. in dry air from GlobalView for 2000 above,
;   (b) climatological monthly mean barometric pressure in 2 by 2 degree 
;       boxes from NCAR, 
;   (c) SST referenced below and, 
;   (d) salinity referenced below.
; 
; PBARO is monthly climatological mean barometric pressure
;   in millibars interpolated onto our 4 degree by 5 degree grid,
;   from a 2 degree by 2 degree monthly data base from NOAA
;   Climate Diagnostic Center Objective Interpolation.
; 
; SST is climatological sea surface temperature, interpolated onto
;   the same 4 degree by 5 degree grid as the pCO2 and flux data. 
;   The source is monthly data on a 1 degree by  1 degree grid from
;   the NOAA Climate Diagnostic Center Objective Interpolation.
;   The years 1981-2006 are averaged with the El Nino months removed,
;   so that the values are consistent with the pCO2 data.
; 
; SAL is salinity, from the NODC World Ocean Database 1998.
; 
; DELTA_PCO2 is seawater - air pCO2 difference (uatm) for a 
;   reference year 2000 computed using PCO2_SW and PCO2_AIR
;   data cited above.
; 
; WINDSPEED is the climatological monthly mean scalar wind speed.
;   Unit is meters per second.  Derived from the 1979-2005 
;   NCEP-DOE AMIP-II Reanalysis (R-2) wind speed data.
; 
; P_VEL_SW06 is the piston velocity.  Unit is centimeter per hour (cm hr^-1).
; 
; GASCOEFSW is the Gas Exchange Coefficient.  
;   Evaluated as:  7.3 * 0.26 * windspeed2 * 
;   (Schmidt Number / 660)^-0.5 * K0_CO2. 
;   Unit is millimole m^-2 month^-1 microatmosphere^-1.
; 
; K0_CO2 is the solubility of CO2 in seawater.  Unit is
;   millimole liter^-1 atmosphere^-1.
; 
; AREA_MKM2 is the area of each box.  Unit is millions of square
;   kilometers, or 10^6 km^2.
; 
; ICE_PCT is the percent of ice coverage.  Where the percentage is
;   below 10% we assume open water.  Where the percentage is greater
;   than 95% we set it = to 95% to account for leads, polynyas, etc.
;   The flux is then multiplied by 1.0 - (ice_pct / 100).
;   Unit is percent.
; 
; FLUXM2SW06 is the Monthly CO2 flux in Moles Carbon per square meter.
;   Unit is Mole Carbon m^-2 month^-1.
; 
; FLUXGMSW06 is the Monthly CO2 flux in grams Carbon per square meter.
;   Unit is grams Carbon m^-2 month^-1.
; 
; TFLUXSW06 is the total mass of Carbon transferred in this
;   box each month, or the area integrated flux.  Unit is 
;   10^15 g Carbon month^-1.
; 
; The SW06 in the names of the variables refers to the parametrization 
;   of the wind speed to gas flux relationship.
; 
begin

; LAT is the latitude of the center of the designated box.  
;   Unit is decimal degrees.
; 
; LON is the longitude of the center of the designated box.
;   Unit is decimal degrees all in EAST longitude.
; 
; MONTH is the calendar month, January = 1, February = 2, etc.
; 
; PCO2_SW is monthly mean pCO2 (uatm) in surface water at SST.
;   The individual pCO2 values measured in different year are
;   normalized to a reference year 2000 using 1.5 uatm per year
;   as the rate of pCO2 change.  The grid 4 degree latitude by
;   5 degree longitude, centered on 2.5, 7.5, etc longitude and
;   0, 4, 8, etc. latitude.
;  
; MOD_SST is based on the individual SST observations that were
;   made concurrently with each pCO2 measurement, and interpolated
;   using  the same interpolation scheme run for the pCO2 data.
;   El Nino months for the area from 10N - 10S and 170W to 100W
;   are excluded from the database before the interpolation is made. 
; 
; VCO2_AIR is the concentration of CO2 in dry air from the GlobalView
;   CO2 database for the week nearest the center date of each month,
;   in 2000.  Unit is parts per million.
; 
; PCO2_AIR is monthly pCO2 (uatm) in air for 2000.  It is computed using:
;   (a) the zonal mean CO2 conc. in dry air from GlobalView for 2000 above,
;   (b) climatological monthly mean barometric pressure in 2 by 2 degree 
;       boxes from NCAR, 
;   (c) SST referenced below and, 
;   (d) salinity referenced below.
; 
; PBARO is monthly climatological mean barometric pressure
;   in millibars interpolated onto our 4 degree by 5 degree grid,
;   from a 2 degree by 2 degree monthly data base from NOAA
;   Climate Diagnostic Center Objective Interpolation.
; 
; SST is climatological sea surface temperature, interpolated onto
;   the same 4 degree by 5 degree grid as the pCO2 and flux data. 
;   The source is monthly data on a 1 degree by  1 degree grid from
;   the NOAA Climate Diagnostic Center Objective Interpolation.
;   The years 1981-2006 are averaged with the El Nino months removed,
;   so that the values are consistent with the pCO2 data.
; 
; SAL is salinity, from the NODC World Ocean Database 1998.
; 
; DELTA_PCO2 is seawater - air pCO2 difference (uatm) for a 
;   reference year 2000 computed using PCO2_SW and PCO2_AIR
;   data cited above.
; 
; WINDSPEED is the climatological monthly mean scalar wind speed.
;   Unit is meters per second.  Derived from the 1979-2005 
;   NCEP-DOE AMIP-II Reanalysis (R-2) wind speed data.
; 
; P_VEL_SW06 is the piston velocity.  Unit is centimeter per hour (cm hr^-1).
; 
; GASCOEFSW is the Gas Exchange Coefficient.  
;   Evaluated as:  7.3 * 0.26 * windspeed2 * 
;   (Schmidt Number / 660)^-0.5 * K0_CO2. 
;   Unit is millimole m^-2 month^-1 microatmosphere^-1.
; 
; K0_CO2 is the solubility of CO2 in seawater.  Unit is
;   millimole liter^-1 atmosphere^-1.
; 
; AREA_MKM2 is the area of each box.  Unit is millions of square
;   kilometers, or 10^6 km^2.
; 
; ICE_PCT is the percent of ice coverage.  Where the percentage is
;   below 10% we assume open water.  Where the percentage is greater
;   than 95% we set it = to 95% to account for leads, polynyas, etc.
;   The flux is then multiplied by 1.0 - (ice_pct / 100).
;   Unit is percent.
; 
; FLUXM2SW06 is the Monthly CO2 flux in Moles Carbon per square meter.
;   Unit is Mole Carbon m^-2 month^-1.
; 
; FLUXGMSW06 is the Monthly CO2 flux in grams Carbon per square meter.
;   Unit is grams Carbon m^-2 month^-1.
; 
; TFLUXSW06 is the total mass of Carbon transferred in this
;   box each month, or the area integrated flux.  Unit is 
;   10^15 g Carbon month^-1.
; 
  ni = 72
  nj = 45

  month = new(12,integer)
  month@long_name = "Month" ;
  month@units = " " ;

  lon = new(ni,double) ;
  lon@long_name = "Longitude" ;
  lon@units = "degrees_east" ;
  lon!0 = "lon"

  lat = new(nj,double)
  lat@long_name = "Latitude" ;
  lat@units = "degrees_north" ;
  lat!0 = "lat"

  PCO2_SW = new((/12,nj,ni/),double)
  PCO2_SW@long_name = "pCO2_sw" ;
  PCO2_SW@units = "\\muatm" ;
  PCO2_SW!0 = "time"
  PCO2_SW!1 = "lat"
  PCO2_SW!2 = "lon"
  

  MOD_SST = new((/12,nj,ni/),double)
  MOD_SST@long_name = "SST observations" ;
  MOD_SST@units = "^\\circC" ;
  MOD_SST!0 = "time"
  MOD_SST!1 = "lat"
  MOD_SST!2 = "lon"
  
  VCO2_AIR = new((/12,nj,ni/),double)
  VCO2_AIR@long_name = "CO_2 in dry air" ;
  VCO2_AIR@units = "ppm" ;
  VCO2_AIR!0 = "time"
  VCO2_AIR!1 = "lat"
  VCO2_AIR!2 = "lon"
  
  PCO2_AIR = new((/12,nj,ni/),double)
  PCO2_AIR@long_name = "pCO_2 Atmosphere" ;
  PCO2_AIR@units = "\\muatm" ;
  PCO2_AIR!0 = "time"
  PCO2_AIR!1 = "lat"
  PCO2_AIR!2 = "lon"
  
  PBARO = new((/12,nj,ni/),double)
  PBARO@long_name = "SLP" ;
  PBARO@units = "millibars" ;
  PBARO!0 = "time"
  PBARO!1 = "lat"
  PBARO!2 = "lon"

  SST = new((/12,nj,ni/),double)
  SST@long_name = "SST Climatology" ;
  SST@units = "^\\circC" ;
  SST!0 = "time"
  SST!1 = "lat"
  SST!2 = "lon"
  
  SAL = new((/12,nj,ni/),double)
  SAL@long_name = "Salinty" ;
  SAL@units = " " ;
  SAL!0 = "time"
  SAL!1 = "lat"
  SAL!2 = "lon"

  DELTA_PCO2 = new((/12,nj,ni/),double)
  DELTA_PCO2@long_name = "\\DeltapCO_2" ;
  DELTA_PCO2@units = "\\muatm" ;
  DELTA_PCO2!0 = "time"
  DELTA_PCO2!1 = "lat"
  DELTA_PCO2!2 = "lon"
  
  WINDSPEED = new((/12,nj,ni/),double)
  WINDSPEED@long_name = "NCEP-DOE wind speed" ;
  WINDSPEED@units = "m s^{-1}" ;
  WINDSPEED!0 = "time"
  WINDSPEED!1 = "lat"
  WINDSPEED!2 = "lon"
  
  P_VEL_SW06 = new((/12,nj,ni/),double)
  P_VEL_SW06@long_name = "Piston velocity" ;
  P_VEL_SW06@units = "m mon^{-1}" ;
  P_VEL_SW06!0 = "time"
  P_VEL_SW06!1 = "lat"
  P_VEL_SW06!2 = "lon"
  
  GASCOEFSW = new((/12,nj,ni/),double)
  GASCOEFSW@long_name = "Gas Exchange Coefficient" ;
  GASCOEFSW@units = "mmol m^{-2} month^{-1} \\muatm^{-1}" ;
  GASCOEFSW!0 = "time"
  GASCOEFSW!1 = "lat"
  GASCOEFSW!2 = "lon"
  
  K0_CO2 = new((/12,nj,ni/),double)
  K0_CO2@long_name = "mmol L^{-1} atm^{-1}" ;
  K0_CO2@units = "CO_2 solubility" ;
  K0_CO2!0 = "time"
  K0_CO2!1 = "lat"
  K0_CO2!2 = "lon"

  AREA_MKM2 = new((/12,nj,ni/),double)
  AREA_MKM2@long_name = "10^6 km^2" ;
  AREA_MKM2@units = "Box area" ;
  AREA_MKM2!0 = "time"
  AREA_MKM2!1 = "lat"
  AREA_MKM2!2 = "lon"
  
  ICE_PCT = new((/12,nj,ni/),double)
  ICE_PCT@long_name = "%" ;
  ICE_PCT@units = "Ice coverage" ;
  ICE_PCT!0 = "time"
  ICE_PCT!1 = "lat"
  ICE_PCT!2 = "lon"

  FLUXM2SW06 = new((/12,nj,ni/),double)
  FLUXM2SW06@long_name = "CO_2 flux" ;
  FLUXM2SW06@units = "mol C m^{-2} month^{-1}" ;
  FLUXM2SW06!0 = "time"
  FLUXM2SW06!1 = "lat"
  FLUXM2SW06!2 = "lon"
  
  FLUXGMSW06 = new((/12,nj,ni/),double)
  FLUXGMSW06@long_name = "CO_2 flux" ;
  FLUXGMSW06@units = "g C m^{-2} month^{-1}" ;
  FLUXGMSW06!0 = "time"
  FLUXGMSW06!1 = "lat"
  FLUXGMSW06!2 = "lon"

  TFLUXSW06 = new((/12,nj,ni/),double)
  TFLUXSW06@long_name = "Total CO2 flux" ;
  TFLUXSW06@units = "Tg C month^{-1}" ;
  TFLUXSW06!0 = "time"
  TFLUXSW06!1 = "lat"
  TFLUXSW06!2 = "lon"
  

  pth = "/glade/p/work/mclong/pco2-ldeo/"
  fnm = pth + "cdiac_orig/ldeo_monthly_clim.month_flux_2006c.txt"
  strhead = "    LAT      LON MONTH PCO2_SW MOD_SST VCO2_AIR PCO2_AIR   PBARO      SST     SAL DELTA_PCO2 WINDSPEED P_VEL_SW06    GASCOEFSW      K0_CO2   AREA_MKM2 ICE_PCT  FLUXM2SW06  FLUXGMSW06    TFLUXSW06"
  
  head = str_split(strhead," ")
  print(head)
  ncol = dimsizes(head)
  strdata = asciiread(fnm,-1,"string")
  nrow = dimsizes(strdata)
  i = -1
  do while (i .le. nrow)
      i = i + 1
      if (strdata(i)) .eq. strhead
          print("header:")
          print(""+strdata(i))
          break
      end if
  end do
  I = i + 1 ; first data row
  ndata = nrow - I
  print("N obs: " + ndata)

  lon = (/ispan(2,357,5)*1d0+0.5/)
  lat = (/ispan(-88,88,4)*1d0/)
  lon&lon = lon
  lat&lat = lat

  do r = 0,ndata-1
      row = stringtodouble(str_split(strdata(r+I)," "))
      l = tointeger(row(2) - 1)
      j = ind(row(0) .eq. lat)
      i = ind(row(1) .eq. lon)

      c = 3
      PCO2_SW(l,j,i)    = (/row(c)/)
      c = c + 1
      MOD_SST(l,j,i)    = (/row(c)/)
      c = c + 1
      VCO2_AIR(l,j,i)   = (/row(c)/)
      c = c + 1
      PCO2_AIR(l,j,i)   = (/row(c)/)
      c = c + 1
      PBARO(l,j,i)      = (/row(c)/)
      c = c + 1
      SST(l,j,i)        = (/row(c)/)
      c = c + 1
      SAL(l,j,i)        = (/row(c)/)
      c = c + 1
      DELTA_PCO2(l,j,i) = (/row(c)/)
      c = c + 1
      WINDSPEED(l,j,i)  = (/row(c)/)
      c = c + 1
      P_VEL_SW06(l,j,i) = (/row(c)/)
      c = c + 1
      GASCOEFSW(l,j,i)  = (/row(c)/)
      c = c + 1
      K0_CO2(l,j,i)     = (/row(c)/)
      c = c + 1
      AREA_MKM2(l,j,i)  = (/row(c)/)
      c = c + 1
      ICE_PCT(l,j,i)    = (/row(c)/)
      c = c + 1
      FLUXM2SW06(l,j,i) = (/row(c)/)
      c = c + 1
      FLUXGMSW06(l,j,i) = (/row(c)/)
      c = c + 1
      TFLUXSW06(l,j,i)  = (/row(c)/)
  end do

  system("rm -f " + pth + "ldeo_monthly_clim_2009_new.nc")
  f = addfile(pth + "ldeo_monthly_clim_2009_new.nc","c")
  
  dpm  = (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)*1d0
  eom = cumsum(dpm,0)
  bom = dpm
  bom(1:11) = cumsum(dpm(0:10),0)+1
  bom(0) = 1.
  time = dim_avg_n( (/bom,eom/),0 )
  time@units = "julian day"
  mmdd = round(2000*10000 + ispan(1,12,1) * 100. + dpm/2,3)
  mmdd!0 = "time"
  mmdd&time = time
  mmdd@format = "YYYYMMDD"

  f->PCO2_SW    = PCO2_SW
  f->MOD_SST    = MOD_SST
  f->VCO2_AIR   = VCO2_AIR
  f->PCO2_AIR   = PCO2_AIR
  f->PBARO      = PBARO
  f->SST        = SST
  f->SAL        = SAL
  f->DELTA_PCO2 = DELTA_PCO2
  f->WINDSPEED  = WINDSPEED
  f->P_VEL_SW06 = P_VEL_SW06
  f->GASCOEFSW  = GASCOEFSW
  f->K0_CO2     = K0_CO2
  f->AREA_MKM2  = AREA_MKM2
  f->ICE_PCT    = ICE_PCT
  f->FLUXM2SW06 = FLUXM2SW06
  f->FLUXGMSW06 = FLUXGMSW06
  f->TFLUXSW06  = TFLUXSW06


  fgco2 = FLUXM2SW06
  fgco2 = (/FLUXM2SW06 / conform(FLUXM2SW06,dpm,0) * 365./)
  fgco2@units = "mol C m^{-2} yr^{-1}"
  f->FLUXMOLM2Y = fgco2

  f->date = mmdd
  f->lat  = lat
  f->lon  = lon
end 