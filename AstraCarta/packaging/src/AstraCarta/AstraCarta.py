from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
import argparse
import numpy as np
import os
import math
import csv
from astropy.io import ascii
import matplotlib.pyplot as plt
import sys
import shutil

def astracarta(**kwargs):

    if kwargs.get("ra") is None:
        raise ValueError('-ra must be supplied.')
    ra = kwargs.get("ra")
    
    if kwargs.get("dec") is None:
        raise ValueError('-dec must be supplied.')
    dec = kwargs.get("dec")

    #coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='fk5')

    if kwargs.get("scale") is None:
        raise ValueError('-scale must be supplied.')
    scale = kwargs.get("scale") #arcsec per pixel

    if kwargs.get("pixwidth") is None:
        raise ValueError('-pixwidth must be supplied.')
    pixwidth = kwargs.get("pixwidth")

    if kwargs.get("pixheight") is None:
        raise ValueError('-pixheight must be supplied.')
    pixheight = kwargs.get("pixheight")

    if kwargs.get("maglimit") is None:
        raise ValueError('-maglimit must be supplied.')
    maglimit = kwargs.get("maglimit")

    buffer = 0;
    if kwargs.get("buffer") is not None:
        buffer = kwargs.get("buffer") #arcminutes

    shape = "rectangle"
    shapenum = 1
    if kwargs.get("shape") is not None:
        if kwargs.get("shape") != "circle":
            if kwargs.get("shape") != "rectangle":
                raise ValueError('-shape may only be "circle" or "rectangle".')
        shape = kwargs.get("shape")
        if (shape == "circle"):
            if pixwidth != pixheight:
                raise ValueError('-shape may only be "circle" if -pixwidth and -pixheight are equal.')
            shapenum = 2
            radius = scale / 3600 * pixwidth / 2 + buffer / 60 #degrees
            radiuspix = radius / (scale/3600)

    rotation = 0
    if kwargs.get("rotation") is not None:
        if (shape == "circle") and kwargs.get("rotation") != 0:
            raise ValueError('-rotation doesn''t make sense with a cirlce query.')
        rotation = kwargs.get("rotation") * np.pi / 180 #radians

    cd11 = -scale/3600 * math.cos(rotation)
    cd12 = scale/3600 * math.sin(rotation)
    cd21 = -scale/3600 * math.sin(rotation)
    cd22 = -scale/3600 * math.cos(rotation)
    crval1 = ra
    crval2 = dec
    crpix1 = pixwidth / 2
    crpix2 = pixheight / 2
    det = 1 / ((cd11 * cd22 - cd12 * cd21) * np.pi / 180)
    cdinv11 = det * cd22
    cdinv12 = -det * cd12
    cdinv21 = -det * cd21
    cdinv22 = det * cd11

    if shape == "rectangle":
        #top left
        xpix_topleft = 1 - (buffer / 60 * 3600 / scale)
        ypix_topleft = 1 - (buffer / 60 * 3600 / scale)
        X_intrmdt = cd11 * (xpix_topleft - crpix1) * np.pi / 180 + cd12 * (ypix_topleft - crpix2) * np.pi / 180
        Y_intrmdt = cd21 * (xpix_topleft - crpix1) * np.pi / 180 + cd22 * (ypix_topleft - crpix2) * np.pi / 180
        ra_topleft = (crval1 * np.pi / 180 + math.atan(X_intrmdt / (math.cos(crval2 * math.pi / 180) - Y_intrmdt * math.sin(crval2 * math.pi / 180)))) * 180 / math.pi
        dec_topleft = (math.asin((math.sin(crval2 * math.pi / 180) + Y_intrmdt * math.cos(crval2 * math.pi / 180)) / math.sqrt(1 + X_intrmdt * X_intrmdt + Y_intrmdt * Y_intrmdt))) * 180 / math.pi
        if ra_topleft < 0:
            ra_topleft += 360;

        #top right
        xpix_topright = pixwidth + (buffer / 60 * 3600 / scale)
        ypix_topright = 1 - (buffer / 60 * 3600 / scale)
        X_intrmdt = cd11 * (xpix_topright - crpix1) * np.pi / 180 + cd12 * (ypix_topright - crpix2) * np.pi / 180
        Y_intrmdt = cd21 * (xpix_topright - crpix1) * np.pi / 180 + cd22 * (ypix_topright - crpix2) * np.pi / 180
        ra_topright = (crval1 * np.pi / 180 + math.atan(X_intrmdt / (math.cos(crval2 * math.pi / 180) - Y_intrmdt * math.sin(crval2 * math.pi / 180)))) * 180 / math.pi
        dec_topright = (math.asin((math.sin(crval2 * math.pi / 180) + Y_intrmdt * math.cos(crval2 * math.pi / 180)) / math.sqrt(1 + X_intrmdt * X_intrmdt + Y_intrmdt * Y_intrmdt))) * 180 / math.pi
        if ra_topright < 0:
            ra_topright += 360;

        #bottom right
        xpix_bottomright = pixwidth + (buffer / 60 * 3600 / scale)
        ypix_bottomright = pixheight + (buffer / 60 * 3600 / scale)
        X_intrmdt = cd11 * (xpix_bottomright - crpix1) * np.pi / 180 + cd12 * (ypix_bottomright - crpix2) * np.pi / 180
        Y_intrmdt = cd21 * (xpix_bottomright - crpix1) * np.pi / 180 + cd22 * (ypix_bottomright - crpix2) * np.pi / 180
        ra_bottomright = (crval1 * np.pi / 180 + math.atan(X_intrmdt / (math.cos(crval2 * math.pi / 180) - Y_intrmdt * math.sin(crval2 * math.pi / 180)))) * 180 / math.pi
        dec_bottomright = (math.asin((math.sin(crval2 * math.pi / 180) + Y_intrmdt * math.cos(crval2 * math.pi / 180)) / math.sqrt(1 + X_intrmdt * X_intrmdt + Y_intrmdt * Y_intrmdt))) * 180 / math.pi
        if ra_bottomright < 0:
            ra_bottomright += 360;

        #bottom left
        xpix_bottomleft = 1 - (buffer / 60 * 3600 / scale)
        ypix_bottomleft = pixheight + (buffer / 60 * 3600 / scale)
        X_intrmdt = cd11 * (xpix_bottomleft - crpix1) * np.pi / 180 + cd12 * (ypix_bottomleft - crpix2) * np.pi / 180
        Y_intrmdt = cd21 * (xpix_bottomleft - crpix1) * np.pi / 180 + cd22 * (ypix_bottomleft - crpix2) * np.pi / 180
        ra_bottomleft = (crval1 * np.pi / 180 + math.atan(X_intrmdt / (math.cos(crval2 * math.pi / 180) - Y_intrmdt * math.sin(crval2 * math.pi / 180)))) * 180 / math.pi
        dec_bottomleft = (math.asin((math.sin(crval2 * math.pi / 180) + Y_intrmdt * math.cos(crval2 * math.pi / 180)) / math.sqrt(1 + X_intrmdt * X_intrmdt + Y_intrmdt * Y_intrmdt))) * 180 / math.pi
        if ra_bottomleft < 0:
            ra_bottomleft += 360;

    catalogue = "GaiaDR3"
    cataloguenum = 1;
    if kwargs.get("catalogue") is not None:
        catalogue = kwargs.get("catalogue")
        if catalogue != "GaiaDR3":
            raise ValueError('Catalogue ' + catalogue + " not valid.")        

    filter = "phot_g_mean_mag"
    filternum = 3
    if kwargs.get("filter") is not None:
        filter = kwargs.get("filter")
        if catalogue == "GaiaDR3":
            if filter == 'bp':
                filter = "phot_bp_mean_mag"
                filternum = 1
            elif filter == 'rp':
                filter = "phot_rp_mean_mag"
                filternum = 2
            elif filter == "g":
                filter = "phot_g_mean_mag"
                filternum = 3
            else:
                raise ValueError('Filter ' + filter + " not valid for Catalogue " + catalogue)

    imageout = False
    if kwargs.get("imageout") is not None:
        imageout = True

    outformat = ".csv"
    fitsout = False
    if kwargs.get("fitsout") is not None:
        outformat = ".fits"
        fitsout = True
    
    imageshow = False;
    if kwargs.get("imageshow") is not None:
        imageshow = True

    outdir = os.getcwd();    
    if kwargs.get("outdir") is not None:
        outdir = kwargs.get("outdir")
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    outdir += "\\"
    rawoutdir = os.path.dirname(os.path.abspath(__file__)) + "\\AstraCartaRawQueries\\"
    if not os.path.isdir(rawoutdir):
        os.makedirs(rawoutdir)

    forcenew = False
    if kwargs.get("forcenew") is not None:
        forcenew = True

    nquery = 500
    if kwargs.get("nquery") is not None:
        nquery = kwargs.get("nquery")    

    pmepoch = 0
    if kwargs.get("pmepoch") is not None:
        pmepoch = kwargs.get("pmepoch")

    pmlimit = np.Infinity
    if kwargs.get("pmlimit") is not None:
        if kwargs.get("pmepoch") is None:
            raise ValueError('pmlimit is only valid if pmepoch is specified')
        pmlimit = kwargs.get("pmlimit")

    entries = "ref_epoch, ra, ra_error, dec, dec_error, pmra, pmra_error, pmdec, pmdec_error, pm, phot_bp_mean_mag, phot_g_mean_mag, phot_rp_mean_mag"
    if kwargs.get("entries") is not None:
         if kwargs.get("entries") == "all":
            entries = "*"
         else:
            entries += " " + kwargs.get("entries")

    notableout = False;
    if kwargs.get("notableout") is not None:
        notableout = True

    rmvrawquery = False;
    if kwargs.get("rmvrawquery") is not None:
        rmvrawquery = True

    if (shape == "circle"):
        rawqueryfilenamehash = hash((ra, dec, nquery, cataloguenum, filternum, radius, shapenum))
        fileoutfilenamehash = hash((ra, dec, nquery, cataloguenum, filternum, radius, shapenum, maglimit, pmepoch, pmlimit))
    else:
        rawqueryfilenamehash = hash((ra, dec, nquery, cataloguenum, filternum, rotation, ra_topleft, dec_topleft, ra_topright, dec_topright, ra_bottomright, dec_bottomright, ra_bottomleft, dec_bottomleft, shapenum))
        fileoutfilenamehash = hash((ra, dec, nquery, cataloguenum, filternum, rotation, ra_topleft, dec_topleft, ra_topright, dec_topright, ra_bottomright, dec_bottomright, ra_bottomleft, dec_bottomleft, shapenum, maglimit, pmepoch, pmlimit))
    rawqueryfilenamehash += sys.maxsize + 1
    fileoutfilenamehash += sys.maxsize + 1

    outname = str(fileoutfilenamehash)
    if kwargs.get("outname") is not None:
        outname = kwargs.get("outname")

    silent = False
    if kwargs.get("silent") is not None:
        silent = True

    rawqueryfilename = rawoutdir + str(rawqueryfilenamehash) + ".csv"
    resultsfilename = outdir + outname + outformat
    imagefilename = outdir + outname + ".jpg"
    if os.path.isfile(resultsfilename) or os.path.isfile(imagefilename):
        f = 1
        resultsfilename = outdir + outname + " ({0})".format(f) + outformat
        imagefilename = outdir + outname + " ({0})".format(f) + ".jpg"
        while os.path.isfile(resultsfilename) or os.path.isfile(imagefilename):
            f += 1
            resultsfilename = outdir + outname + " ({0})".format(f) + outformat
            imagefilename = outdir + outname + " ({0})".format(f) + ".jpg"

    if not os.path.isfile(rawqueryfilename) or forcenew:    #query new table download
        if not silent:
            print("Launching " + shape + " query to " + catalogue + "...")
            sys.stdout.flush()
        if  catalogue == "GaiaDR3":
            if nquery == 0:
                jobstr = "SELECT " + entries + " FROM gaiadr3.gaia_source\n"
            else:
                jobstr = "SELECT TOP {0}".format(nquery) + " " + entries + " FROM gaiadr3.gaia_source\n"
            jobstr += "WHERE 1=CONTAINS(POINT('ICRS', gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),"

            if shape == "rectangle":
                jobstr += "POLYGON('ICRS',{0},{1},{2},{3},{4},{5},{6},{7}))\n".format(ra_topleft,dec_topleft,ra_topright,dec_topright,ra_bottomright,dec_bottomright,ra_bottomleft,dec_bottomleft)
            else:
                jobstr += "CIRCLE('ICRS',{0},{1},{2}))\n".format(ra,dec,radius)

            jobstr += "ORDER by gaiadr3.gaia_source." + filter + " ASC"

            #print(jobstr)
            job = Gaia.launch_job_async(jobstr, dump_to_file=False)
            #print(job)

        results = job.get_results()
    
        if fitsout:
            # Strip object columns from FITS table
            removelist = []
            for col in results.columns:
                if results[col].dtype == 'object':
                    removelist += [col]
            results.remove_columns(removelist)

        results.write(rawqueryfilename, overwrite=True, format="csv")


    rowswritten = 0
    with open(rawqueryfilename, 'r') as rawquery:
        with open(resultsfilename, 'w', newline='') as results:
            reader = csv.DictReader(rawquery)
            writer = csv.DictWriter(results, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                if float(row[filter]) <= maglimit: #valid to write
                    if pmepoch != 0: #update ra and dec
                        if row['pm'] != '' and float(row['pm']) < pmlimit:
                            dt = pmepoch - float(row['ref_epoch'])
                            row['ref_epoch'] = float(row['ref_epoch']) + dt
                            row['ra'] = float(row['ra']) + dt * float(row['pmra']) / 3600 / 1000
                            row['dec'] = float(row['dec']) + dt * float(row['pmdec']) / 3600 / 1000
                            writer.writerow(row)
                            rowswritten += 1;
                    else:
                        writer.writerow(row)
                        rowswritten += 1;

    if rmvrawquery:
        shutil.rmtree(rawoutdir, ignore_errors=True)
    
    if rowswritten == 0:
        os.remove(resultsfilename)
        if not silent:
            print("No sources found.")
        return ""
    else:
        if not silent:
            print("Found {0} sources.".format(rowswritten))
            sys.stdout.flush()

        if imageshow or imageout:
            #plot RA, DEC locations of sources on sky-coordinate image...
            #this needs to plot on a proper sky grid...ok when near equator but no good near poles
            #ras = np.array([])
            #dec = np.array([])
            #with open(resultsfilename, 'r') as datafile:
            #    reader = csv.DictReader(datafile)
            #    for row in reader:
            #        ras = np.append(ras, float(row['ra']))
            #        dec = np.append(dec, float(row['dec']))
            #fig0 = plt.figure(0,figsize=(6,6))
            #axes0 = fig0.add_subplot(111)
            #axes0.scatter(ras,dec)
            #plt.savefig(imagefilename)
            #if imageshow:
                #plt.show()
                #plt.close()

            #plot x,y pixel locations of sources on detector image...
            x = np.array([])
            y = np.array([])
            a0 = crval1 * np.pi / 180
            d0 = crval2 * np.pi / 180
            with open(resultsfilename, 'r') as datafile:
                reader = csv.DictReader(datafile)
                for row in reader:
                    ra = float(row['ra'])
                    de = float(row['dec'])
                    a = ra * np.pi / 180
                    d = de * np.pi / 180
                    X_intrmdt = math.cos(d) * math.sin(a - a0) / (math.cos(d0) * math.cos(d) * math.cos(a - a0) + math.sin(d0) * math.sin(d))
                    Y_intrmdt = (math.cos(d0) * math.sin(d) - math.cos(d) * math.sin(d0) * math.cos(a - a0)) / (math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(a - a0))
                    X_pix = cdinv11 * X_intrmdt + cdinv12 * Y_intrmdt + crpix1
                    Y_pix = cdinv21 * X_intrmdt + cdinv22 * Y_intrmdt + crpix2
                    x = np.append(x, X_pix)
                    y = np.append(y, Y_pix)
            fig0 = plt.figure(0,figsize=(6,6))
            axes0 = fig0.add_subplot(111)
            axes0.scatter(x,y)
            xlimmin = 1
            xlimmax = pixwidth
            ylimmin = 1
            ylimmax = pixheight
            if shape == "rectangle":
                if buffer > 0:
                    axes0.plot([1, 1, pixwidth, pixwidth, 1], [1, pixheight, pixheight, 1, 1], 'k')
                    xlimmin = xpix_topleft
                    xlimmax = xpix_topright
                    ylimmin = ypix_topleft
                    ylimmax = ypix_bottomleft
                elif buffer < 0:
                    axes0.plot([xpix_topleft, xpix_topright, xpix_bottomright, xpix_bottomleft, xpix_topleft], [ypix_topleft, ypix_topright, ypix_bottomright, ypix_bottomleft, ypix_topleft], 'k')
            elif shape == "circle":
                theta = np.linspace(0,360,361)
                if buffer >= 0:
                    rimage = crpix1                    
                    if buffer > 0:
                        xlimmin -= (radiuspix - rimage)
                        xlimmax += (radiuspix - rimage)
                        ylimmin -= (radiuspix - rimage)
                        ylimmax += (radiuspix - rimage)
                elif buffer < 0:
                    rimage = radiuspix
                x = rimage*np.sin(np.radians(theta)) + crpix1
                y = rimage*np.cos(np.radians(theta)) + crpix2
                axes0.plot(x, y, 'k')

            plt.title('Image Plot of Sources Found with MagLimit of ' + str(maglimit))
            plt.ylabel('Vertical Image Axis (Pixels)')
            plt.xlabel('Horizontal Image Axis (Pixels)')
            axes0.set_xlim([xlimmin, xlimmax]);
            axes0.set_ylim([ylimmax, ylimmin]);

            if imageout:
                plt.savefig(imagefilename)
                if not silent:
                    print("Wrote output image to: " + imagefilename)
                    sys.stdout.flush()
        
        if notableout:
            os.remove(resultsfilename)
        elif fitsout:
            t = ascii.read(resultsfilename, format='csv')
            os.remove(resultsfilename)
            base = os.path.splitext(resultsfilename)[0]
            resultsfilename = base + ".fits"
            t.write(resultsfilename, overwrite=True, format='fits')
        if not notableout:
            if not silent:
                print("Wrote output table to: " + resultsfilename)
                sys.stdout.flush()

        if imageshow:
            plt.show()
            plt.close()   
           
        if notableout:
            return ""

        return resultsfilename

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Search for sources in a region which satisfy a magnitude limit, and return outputs in a table and optionally an image. See https://github.com/user29A/AstraCarta/wiki for more info.')
    parser.add_argument('-ra', dest='ra', action='store', const=None, default=None, type=float, help='MANDATORY: RA of field center (decimal degrees, J2000).')
    parser.add_argument('-dec', dest='dec', action='store', const=None, default=None, type=float, help='MANDATORY: Dec of field center (decimal degrees, J2000).')
    parser.add_argument('-scale', dest='scale', action='store', const=None, default=None, type=float, help='MANDATORY: Plate scale of image. Identical to CDELT of WCS (arcseconds per pixel).')
    parser.add_argument('-pixwidth', dest='pixwidth', action='store', const=None, default=None, type=int, help='MANDATORY: Width of image in pixels.')
    parser.add_argument('-pixheight', dest='pixheight', action='store', const=None, default=None, type=int, help='MANDATORY: Height of image in pixels.')
    parser.add_argument('-maglimit', dest='maglimit', action='store', const=None, default=None, type=float, help='MANDATORY: Magnitude limit below which to flag bright sources and save to output table. Pass an unreasonably low magnitude (30, say) to pass all.')
    parser.add_argument('-buffer', dest='buffer', action='store', const=None, default=None, type=float, help='OPTIONAL: Tolerance buffer around image field, in arcminutes. This field can be negative, if one wishes to mitigate image padding in the query.')
    parser.add_argument('-shape', dest='shape', action='store', const=None, default=None, type=str, help='OPTIONAL: Shape of field to query: "rectangle" (default) or "circle". Circle may only be used if pixwidth and pixheight are equal. Rectangle query uses a polygon query with corners defined by an ad-hoc WCS given the supplied field parameters, whereas circle uses a radius.')
    parser.add_argument('-rotation', dest='rotation', action='store', const=None, default=None, type=float, help='OPTIONAL: Field rotation, applicable to a rectangle query. Raises an exception if used for a circle query.')
    parser.add_argument('-catalogue', dest='catalogue', action='store', const=None, default=None, type=str, help='OPTIONAL: Catalogue or Service to query. Valid options are currently: "GaiaDR3"')
    parser.add_argument('-filter', dest='filter', action='store', const=None, default=None, type=str, help='OPTIONAL: Filter of the catalogue to sort on. Options are: for GaiaDR3: ''rp'', ''bp'', ''g'' (default)')
    parser.add_argument('-forcenew', dest='forcenew', action='store_true', default=None, help='OPTIONAL: Force new astroquery. The raw query is saved with a filename based on a hash of the astroquery parameters, and therefore should be unique for unique queries, and the same for the same queries. The exception is for the "entries" option which cannot be hashed non-randomly. Therefore if everything else stays the same except for "entries", one would need to force a new query.')
    parser.add_argument('-imageout', dest='imageout', action='store_true', default=None, help='OPTIONAL: Output an image field plot with maglimit sources marked.')
    parser.add_argument('-imageshow', dest='imageshow', action='store_true', default=None, help='OPTIONAL: Show the image field plot.')
    parser.add_argument('-outdir', dest='outdir', action='store', const=None, default=None, type=str, help='OPTIONAL: Directory to save files. By default files are saved in the current working directory.')
    parser.add_argument('-outname', dest='outname', action='store', const=None, default=None, type=str, help='OPTIONAL: The name to use for output files. If not supplied a settings-consistent but random hash will be used for output file names. Existing filenames will not be overwritten.')
    parser.add_argument('-fitsout', dest='fitsout', action='store_true', default=None, help='OPTIONAL: Output results table format as FITS binary table (instead of csv).')
    parser.add_argument('-rmvrawquery', dest='rmvrawquery', action='store_true', default=None, help='OPTIONAL: Remove the raw query folder and its contents after running. This will force future astroqueries.')
    parser.add_argument('-nquery', dest='nquery', action='store', const=None, default=None, type=int, help='OPTIONAL: Number of brightest sources in the filter to retreive from the query service. Pass 0 to retreive all sources. Default 500.')
    parser.add_argument('-pmepoch', dest='pmepoch', action='store', const=None, default=None, type=float, help='OPTIONAL: Pass the year.year value of the observation to update the RA and Dec entries of the table with their proper motion adjustments, given the catalogue reference epoch. Only entries in the query which have valid proper motion entries will be saved to output.')
    parser.add_argument('-pmlimit', dest='pmlimit', action='store', const=None, default=None, type=float, help='OPTIONAL: Limit the output to proper motions whose absolute value are less than pmlimit. Milliarcseconds per year.')
    parser.add_argument('-entries', dest='entries', action='store', const=None, default=None, type=str, help='OPTIONAL: A commaspace ", " separated list of source columns to request from the query. Pass entries="all" to get everything from the query source. Default is for GaiaDR3, entries="ref_epoch, ra, ra_err, dec, dec_err, pmra, pmra_err, pmdec, pmdec_err, pm, phot_bp_mean_mag, phot_g_mean_mag, phot_rp_mean_mag". Thus, if entries is supplied, it appends additional entries to the default. For example if you additionally wanted the absolute value of proper motion errors then passing entries="pm_error" would append " pm_error" to the string.')
    parser.add_argument('-notableout', dest='notableout', action='store_true', default=None, help='OPTIONAL: Do not write an output file even when sources have been found. Useful if only wanting to view images but don''t want to fill up a directory with table results.')
    parser.add_argument('-silent', dest='verbose', action='store_false', default=None, help='OPTIONAL: Do not output process milestones to command window. Default false.')
    args = parser.parse_args()
    astracarta(**vars(args))

