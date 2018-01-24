load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
setfileoption("nc","Format","LargeFile")

fv = "9.9692099683868690e+36"

begin

  pth = "/glade/p/work/mclong/orcas/fltdata-merge"

  file_list = (/"ORCASall.merge60"/) ;","ORCASall.mergeTOGA","ORCASall.merge60"/) ;,

  now = systemfunc("date +%Y%m%d")

;------------------------------------------------------------------
; read metadata file
;------------------------------------------------------------------

  strmeta := asciiread("gv-vars-metadata.txt",-1,"string")
  nvar_meta = dimsizes(strmeta)

  units = True
  long_name = True
  do i = 0,nvar_meta-1
      str = str_strip(str_split(strmeta(i),"|"))
      if str(1) .ne. ""
          units@$str(0)$ = str(1)
      end if
      if str(2) .ne. ""
          long_name@$str(0)$ = str(2)
      end if
  end do

  ;; some ammendments
  long_name@Year ="Year"
  units@Year = "Year"
  
  long_name@flt = "Flight number"
  
  long_name@DOY = "Day of year"
  units@DOY = "days"

  long_name@UTC = "seconds since midnight on UTC day of takeoff"
  units@UTC = "seconds"

  units@U_FMPD = "m/s"
  long_name@U_FMPD = "Zonal wind"

  units@V_FMPD = "m/s"
  long_name@V_FMPD = "Meridional wind"

  units@Q_FMPD = "kg/kg"
  long_name@Q_FMPD = "Specific humidity"

  units@PS_FMPD = "Pa"
  long_name@PS_FMPD = "Surface pressure"

  units@CO2_AO2 = "ppm"
  long_name@CO2_AO2 = "CO2 (AO2)"

  units@O2_AO2 = "per meg"
  long_name@O2_AO2 = "O2 (AO2)"

  units@APO_AO2 = "per meg"
  long_name@APO_AO2 = "APO (AO2"

  units@CWC = "g/m^3"
  long_name@CWC = "Condensed water concentration"

  units@O2_LND_FMPD = "per meg"
  long_name@O2_LND_FMPD = "O2_LND (CESM)"

  units@O2_FFF_FMPD = "per meg"
  long_name@O2_FFF_FMPD = "O2_FFF (CESM)"

  units@O2_FMPD = "per meg"
  long_name@O2_FMPD = "O2 (total,CESM)"

  units@APO_OCN_FMPD = "per meg"
  long_name@APO_OCN_FMPD = "APO_OCN (CESM)"

  units@APO_FFF_FMPD = "per meg"
  long_name@APO_FFF_FMPD = "APO_FFF (CESM)"

  units@APO_FMPD = "per meg"
  long_name@APO_FMPD = "APO (total,CESM)"

  units@DMS_FMPD = "ppb"
  long_name@DMS_FMPD = "DMS (model)"

  units@DMS = "ppt"
  long_name@DMS = "DMS (TOGA)"


  do file_i = 0,dimsizes(file_list)-1
      file_base = file_list(file_i)
      print("======================")
      print("working on: "+file_base)

      file_in = pth+"/"+file_base+".tbl"
      file_out = pth+"/"+file_base+"."+now+".nc"

      ;;------------------------------------------------------------------
      ;; read data file
      ;;------------------------------------------------------------------
      col := readAsciiHead(file_in,1)
      col := str_split(col," ")
      ncol = dimsizes(col)
      col = str_sub_str(col,".","_")

      strdata := asciiread(file_in,-1,"string")
      nrow = dimsizes(strdata)-1
      
      strarrary := new(ncol,string)
      data := new((/nrow,ncol/),double)
      do i = 1,dimsizes(strdata)-1
          strarray := str_split(strdata(i)," ")
          strarray = where(strarray .eq. "NA",fv,strarray)
          strarray = where(strarray .eq. "-888",fv,strarray)
          data(i-1,:) = stringtodouble(strarray)
      end do

      ;;------------------------------------------------------------------
      ;; open output file
      ;;------------------------------------------------------------------
      system("rm -f " + file_out)
      fo = addfile(file_out, "c")
      modified = systemfunc("stat -c '%y' "+file_in)
      
      fAtt = True
      fAtt@source =  "converted from "+file_in+" on "+systemfunc("date")
      fAtt@source_timestamp = "source file was last modified: "+modified
      fileattdef(fo,fAtt)
 
      filedimdef(fo,"time",nrow,True)
      filevardef(fo,"time","double","time")
      do j = 0,ncol-1
          filevardef(fo,col(j),"double","time")
      end do
      
      ;;------------------------------------------------------------------
      ;; write data to file
      ;;------------------------------------------------------------------
      do j = 0,ncol-1
          fo->$col(j)$ = (/data(:,j)/)
          if isatt(units,col(j))
              fo->$col(j)$@units = units@$col(j)$
          else
              print("units@"+col(j)+" = ?")
          end if
          
          if isatt(long_name,col(j))
              fo->$col(j)$@long_name = long_name@$col(j)$
          else
              print("long_name@"+col(j)+" = ?")
          end if

      end do
      time := (fo->DOY-1.)*86400. + fo->UTC
      time@units = "seconds since 2016-01-01 00:00"
      time@calendar = "gregorian"
      time!0 = "time"
      
      time := cd_convert(time,"days since 2016-01-01 00:00")
      
      fo->time = (/time/)
      att = getvaratts(time)
      do i = 0,dimsizes(att)-1
          fo->time@$att(i)$ = time@$att(i)$
      end do


      delete(fo)
      print("done.")
      print("")
  end do

end