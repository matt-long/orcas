loadpath = systemfunc("dirname "+get_script_name())
loadscript(loadpath+"/regrid_var.ncl")
begin
  if task .eq. "gen_latlon_grid_file"
      opt = True
      opt@ForceOverwrite = True
      opt@Debug = True
      opt@LLCorner = (/ -90.0+0.5*dlat,       left_lon_corner+0.5*dlon /)
      opt@URCorner = (/  90.0-0.5*dlat, 360.0+left_lon_corner-0.5*dlon /)
      latlon_to_SCRIP(grid_out_fname, grid_type, opt)

      f = addfile(grid_out_fname,"w")
      nx = 360/dlon
      ny = 180/dlat
      lat = fspan(opt@LLCorner(0),opt@URCorner(0),toint(ny))
      lon = fspan(opt@LLCorner(1),opt@URCorner(1),toint(nx))

      lat!0 = "lat"
      lon!0 = "lon"
      lat&lat = lat
      lon&lon = lon
      f->lat = lat
      f->lon = lon

  else if task .eq. "gen_rectilinear_grid_file"
      opt = True
      opt@ForceOverwrite = True
      opt@Debug = True

      f = addfile(latlon_file,"r")
      lat = f->lat
      lon = f->lon

      lat!0 = "lat"
      lon!0 = "lon"
      lat&lat = lat
      lon&lon = lon

      rectilinear_to_SCRIP(grid_out_fname, lat, lon, opt)
      f = addfile(grid_out_fname,"w")
      f->lat = lat
      f->lon = lon

  else if task .eq. "gen_weight_file"
      opt = True
      opt@ForceOverwrite = True
      opt@InterpMethod = InterpMethod
      opt@Check = True
      opt@RemovePETLog = True
      opt@Debug = True
      ESMF_regrid_gen_weights(srcGridFile, dstGridFile, wgtFile, opt)

  else if task .eq. "regrid_var"
      opt = True
      opt@Debug = True
      opt@src_grid = src_grid
      opt@dst_grid = dst_grid
      opt@postfill_opt = postfill_opt
      opt@prefill_opt = prefill_opt
      opt@outfile_opt = outfile_opt
      regrid_var(wgtFile, fname_in, varname_in, time_coordname, \
                 depth_coordname, vert_grid_file, fname_out, varname_out, opt)
  else
      print("Task unknown: "+task)
      status_exit(1)
  end if
  end if
  end if
  end if

  system("rm -f "+ABNORMAL_EXIT)
  exit
end
status_exit(1)
