get_attribute_info <- function(attribute, level=NULL)
{
  info = list()
  
  # Surface attributes definition: short-name, long-name, unit
  info[["cape"]] <- c("cape_fc", "interpolated cape ensemble forecast", "J kg^-1")
  info[["sp"]]   <- c("sp_fc", "Pa", "interpolated surface pressure ensemble forecast")
  info[["msl"]]  <- c("msl_fc", "Pa", "interpolated mean sea level pressure ensemble forecast")
  info[["u10"]]  <- c("u10_fc", "m s^-1", "interpolated u-wind ensemble forecast")
  info[["v10"]]  <- c("v10_fc", "m s^-1", "interpolated v-wind ensemble forecast")
  info[["t2m"]]  <- c("t2m_fc", "K", "interpolated 2m temperature ensemble forecast")
  info[["skt"]]  <- c("skt_fc", "K", "interpolated skin temperature ensemble forecast")
  info[["d2m"]]  <- c("d2m_fc", "K", "interpolated 2m dew point temperature ensemble forecast")
  info[["sshf"]] <- c("sshf_fc", "W m^-2", "interpolated sensible heat flux ensemble forecast")
  info[["slhf"]] <- c("slhf_fc", "W m^-2", "interpolated latent heat flux ensemble forecast")
  info[["ssr"]]  <- c("ssr_fc", "W m^-2", "interpolated net short wave radiation flux ensemble forecast")
  info[["str"]]  <- c("str_fc", "W m^-2", "interpolated net long wave radiation flux ensemble forecast")
  info[["sm"]]   <- c("sm_fc", "kg m^-3", "interpolated soil moisture ensemble forecast")
  info[["tp"]]   <- c("tp_fc", "kg m^-2", "interpolated total precipitation ensemble forecast")
  info[["tcc"]]  <- c("tcc_fc", "%", "interpolated total cloud cover ensemble forecast")

  # Metadata attributes definition:
  info[["orog"]]  <- c("orog", "m", "interpolated orography information")

  # Pressure level 500hPa attributes definition: 
  if (!is.null(level)){
    info[["t"]]  <- c(paste0("t_pl", level,"_fc"), "K", paste0("interpolated temperature at pressure level ", level,"hPa ensemble forecast"))
    info[["u"]]  <- c(paste0("u_pl", level,"_fc"), "m s^-1", paste0("interpolated u-wind at pressure level ", level,"hPa ensemble forecast"))
    info[["v"]]  <- c(paste0("v_pl", level,"_fc"), "m s^-1", paste0("interpolated v-wind at pressure level ", level,"hPa ensemble forecast"))
    info[["q"]]  <- c(paste0("q_pl", level,"_fc"), "kg kg^-1", paste0("interpolated specific humidity at pressure level ", level,"hPa ensemble forecast"))
    info[["gh"]] <- c(paste0("gh_pl", level,"_fc"), "gpm", paste0("interpolated geopotential height at pressure level ", level,"hPa ensemble forecast"))
  }

  return(info[[attribute]])
}
