#!/usr/bin/env python
import os
import numpy as np
from osgeo import osr
from osgeo import gdal
from netCDF4 import Dataset
osr.UseExceptions()


def netcdf_singleband_float_cog(input_file, output_file, band, lat, lon):
    """
    Convert NetCDF data to single-band Cloud Optimised GeoTiff

    Parameters
    ----------

    input_file
        Path to input NetCDF file
    output_file
        Path to output .tiff file
    band
        Variables name in the NetCDF file.
    lat
        Name of latitude variable in NetCDF file
    lon
        Name of longitude variable in NetCDF file
    """
    temp_file = os.path.join(os.path.dirname(output_file), "temp_" + os.path.basename(output_file))
    with Dataset(input_file, "r", format="NETCDF4") as nc:
        lat = np.array(nc.variables[lat][:])
        lon = np.array(nc.variables[lon][:])
        band = np.array(nc.variables[band][:])
        alpha = np.zeros_like(band)
        alpha[~np.isnan(band)] = 255

        image_size = band.shape
        nx = image_size[0]
        ny = image_size[1]
        xmin, ymin, xmax, ymax = [np.nanmin(lon), np.nanmin(lat), np.nanmax(lon), np.nanmax(lat)]
        xres = (xmax - xmin) / float(ny)
        yres = (ymax - ymin) / float(nx)
        geotransform = (xmin, xres, 0, ymax, 0, -yres)

        dst_ds = gdal.GetDriverByName('GTiff').Create(temp_file, ny, nx, 2, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dst_ds.SetProjection(srs.ExportToWkt())
        dst_ds.GetRasterBand(1).WriteArray(band)
        dst_ds.GetRasterBand(2).WriteArray(alpha)
        dst_ds.FlushCache()
        dst_ds = None

        translate_options = gdal.TranslateOptions(gdal.ParseCommandLine(
            '-co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE'))
        ds = gdal.Open(temp_file, gdal.OF_READONLY)
        gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
        ds.BuildOverviews('NEAREST', [2, 4, 8, 16, 32])
        ds = None
        del ds
        ds = gdal.Open(temp_file)
        gdal.Translate(output_file, ds, options=translate_options)
        ds = None
        os.unlink(temp_file)
        os.unlink(temp_file + ".ovr")


def netcdf_singleband_color_cog(input_file, output_file, band, lat, lon, color_palette="white-black", bounds=False):
    """
        Convert NetCDF data to single-band Cloud Optimised GeoTiff

        Parameters
        ----------

        input_file
            Path to input NetCDF file
        output_file
            Path to output .tiff file
        band
            Variables name in the NetCDF file.
        lat
            Name of latitude variable in NetCDF file
        lon
            Name of longitude variable in NetCDF file
        color_palette
            | **Default: False**
            | Color palette to use for single-band outputs. Defaults to a white-black color palette.
        bounds
            | **Default: False**
            | Bounds [lower,upper] to use for single-band color palette transformation, default to min and max of band.
        """
    temp_file = os.path.join(os.path.dirname(output_file), "temp_" + os.path.basename(output_file))
    with Dataset(input_file, "r", format="NETCDF4") as nc:
        lat = np.array(nc.variables[lat][:])
        lon = np.array(nc.variables[lon][:])

        band = np.array(nc.variables[band][:])
        if np.isnan(band).all():
            return False
        if bounds:
            lower = bounds[0]
            upper = bounds[1]
        else:
            lower = np.nanmin(band)
            upper = np.nanmax(band)

        if color_palette == "white-black":
            r_pixels = np.around(((band - lower) / (upper - lower)) * 255)
            b_pixels = r_pixels.copy()
            g_pixels = r_pixels.copy()
        else:
            raise ValueError("Color palette {} not recognised".format(color_palette))

        image_size = r_pixels.shape
        a_pixels = np.zeros_like(r_pixels)
        a_pixels[~np.isnan(r_pixels)] = 255

        nx = image_size[0]
        ny = image_size[1]
        xmin, ymin, xmax, ymax = [np.nanmin(lon), np.nanmin(lat), np.nanmax(lon), np.nanmax(lat)]
        xres = (xmax - xmin) / float(ny)
        yres = (ymax - ymin) / float(nx)
        geotransform = (xmin, xres, 0, ymax, 0, -yres)

        dst_ds = gdal.GetDriverByName('GTiff').Create(temp_file, ny, nx, 4, gdal.GDT_Byte)
        dst_ds.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dst_ds.SetProjection(srs.ExportToWkt())
        dst_ds.GetRasterBand(1).WriteArray(r_pixels)
        dst_ds.GetRasterBand(2).WriteArray(g_pixels)
        dst_ds.GetRasterBand(3).WriteArray(b_pixels)
        dst_ds.GetRasterBand(4).WriteArray(a_pixels)
        dst_ds.FlushCache()
        dst_ds = None

        translate_options = gdal.TranslateOptions(gdal.ParseCommandLine(
            '-co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE'))
        ds = gdal.Open(temp_file, gdal.OF_READONLY)
        gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
        ds.BuildOverviews('NEAREST', [2, 4, 8, 16, 32])
        ds = None
        del ds
        ds = gdal.Open(temp_file)
        gdal.Translate(output_file, ds, options=translate_options)
        ds = None
        os.unlink(temp_file)
        os.unlink(temp_file + ".ovr")


def netcdf_rgb_cog(input_file, output_file, red_band, green_band, blue_band, lat, lon):
    """
        Convert NetCDF data to Cloud Optimised GeoTiff

        Parameters
        ----------

        input_file
            Path to input NetCDF file
        output_file
            Path to output .tiff file
        red_band
            Variables name in the NetCDF file.
        green_band
            Variables name in the NetCDF file.
        blue_band
            Variables name in the NetCDF file.
        lat
            Name of latitude variable in NetCDF file
        lon
            Name of longitude variable in NetCDF file
        """
    temp_file = os.path.join(os.path.dirname(output_file), "temp_" + os.path.basename(output_file))
    with Dataset(input_file, "r", format="NETCDF4") as nc:
        lat = np.array(nc.variables[lat][:])
        lon = np.array(nc.variables[lon][:])

        red = np.array(nc.variables[red_band][:])
        green = np.array(nc.variables[green_band][:])
        blue = np.array(nc.variables[blue_band][:])
        if np.isnan(red).all() or np.isnan(green).all() or np.isnan(blue).all():
            return False
        r_pixels = np.around(((red - np.nanmin(red)) / (np.nanmax(red) - np.nanmin(red))) * 255)
        g_pixels = np.around(((green - np.nanmin(green)) / (np.nanmax(green) - np.nanmin(green))) * 255)
        b_pixels = np.around(((blue - np.nanmin(blue)) / (np.nanmax(blue) - np.nanmin(blue))) * 255)

        image_size = r_pixels.shape
        a_pixels = np.zeros_like(r_pixels)
        a_pixels[~np.isnan(r_pixels)] = 255

        nx = image_size[0]
        ny = image_size[1]
        xmin, ymin, xmax, ymax = [np.nanmin(lon), np.nanmin(lat), np.nanmax(lon), np.nanmax(lat)]
        xres = (xmax - xmin) / float(ny)
        yres = (ymax - ymin) / float(nx)
        geotransform = (xmin, xres, 0, ymax, 0, -yres)

        dst_ds = gdal.GetDriverByName('GTiff').Create(temp_file, ny, nx, 4, gdal.GDT_Byte)
        dst_ds.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dst_ds.SetProjection(srs.ExportToWkt())
        dst_ds.GetRasterBand(1).WriteArray(r_pixels)
        dst_ds.GetRasterBand(2).WriteArray(g_pixels)
        dst_ds.GetRasterBand(3).WriteArray(b_pixels)
        dst_ds.GetRasterBand(4).WriteArray(a_pixels)
        dst_ds.FlushCache()
        dst_ds = None

        translate_options = gdal.TranslateOptions(gdal.ParseCommandLine(
            '-co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE'))
        ds = gdal.Open(temp_file, gdal.OF_READONLY)
        gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
        ds.BuildOverviews('NEAREST', [2, 4, 8, 16, 32])
        ds = None
        del ds
        ds = gdal.Open(temp_file)
        gdal.Translate(output_file, ds, options=translate_options)
        ds = None
        os.unlink(temp_file)
        os.unlink(temp_file + ".ovr")
