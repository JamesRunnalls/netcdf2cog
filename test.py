from netcdf2cog import netcdf_singleband_float_cog, netcdf_singleband_color_cog, netcdf_rgb_cog

#netcdf_singleband_float_cog("input/polymer.nc", "output/polymer_tsm_singleband_float.tiff", "tsm_vantrepotte665", "latitude", "longitude")
#netcdf_singleband_color_cog("input/polymer.nc", "output/polymer_tsm_singleband_color.tiff", "tsm_vantrepotte665", "latitude", "longitude")
#netcdf_rgb_cog("input/polymer.nc", "output/polymer_tsm_rgb.tiff", "Rw400", "Rw560", "Rw779", "latitude", "longitude")


netcdf_singleband_float_cog("/media/jamesrunnalls/JamesSSD/Eawag/DIAS/output_data/datalakes_sui_S2_geneva_2021-09-06_2021-09-06/L2ACOLITE/L2ACOLITE_S2A_MSIL1C_20210906T103021_N0301_R108_T31TGM_20210906T141939.SAFE.nc", "output/s2_tsm_singleband_float.tiff", "TUR_Nechad2016_665", "lat", "lon")
