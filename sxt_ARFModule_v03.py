#/bin/python

import numpy as np
import os
from astropy.io import fits 
import matplotlib.pyplot as plt


'''
#+++++++++++++++++++++++++++++++++++++++++++++++++++ Readme +++++++++++++++++++++++++++++++++++++++++++++

Plese use this script for estimating the useful fraction to correct the normalization of your spectral fitting or overall ARF (using recommneded tools only)

===========================================================
=========================================================================
	 Running Task: /home/chandra/Desktop/SXTEEFModule/sxt_ARFModule_v03.py 
 	 	 Version: 0.03 Release Date: 2019-07-03 Updated on: 2022-01-05 
	 	 Developed By : Sunil Chandra, 
             SAAO, Cape Town, South Africa
 	         CSR, North-West University,
                 Last Updated on 05, January, 2022 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The estimations provided by this script is based on actual SXT observations. The profiles used here are corrected for background counts. The values estimated from this script is deemed to have 2-3 % uncertainities. This script is best suited for point sources. We may expect ~5 % uncertainity in extracting ARF for sources with fluxes ~ 1 Crab or brighter... 
 The parametric way of estimating the factor uses  kings's profile as functional form + best fit archived parameters
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
=========================================================================

Usage: sxt_ARFModule_v02.py [options] 

Options:
  -h, --help            show this help message and exit
  -m MODE, --mode=MODE  Mode for extraction;   Two options a) using table and
                        b) using fit parameters   Default is b)   Input
                        options are ['a', 'tbl', 'tab', 't', 'table'] for a)
                        and one of ['b', 'par', 'param', 'fit', 'fitpar'] for
                        b)
  -r RADIUS, --radius=RADIUS
                        Source Region Size in arcmin.. Always use for cicular
                        regions...  Multiple radii (r1, r2, r3) can be entered
                        in format or 'r1nr2nr3' or 'r1ANDr2ANDr3' or
                        'r1ORr2ORr3'  Default is 15 arcm
  --sxtarf=SXTARF       The Full Path of Base ARF provided by SXT Team...
                        For Example 'sxt_pc_excl00_v04.arf'
  -e EVTFILE, --evtfile=EVTFILE
                        The Events file in case you are entering source
                        postion in sky pixels instead of detector
                        coordinates.....Default = None
  --sxtpha=SXTPHA       The Name of the input SXT spectrum [needed to fetch
                        source region related information]...   For Example
                        'sxt_pc_mytarget_src.pha'  Default is None..[This
                        means source coordinates should be entered manually]
  -x XPIX, --xpix=XPIX   X coordinate (RA) of the source position [sky/RAW
                        pixels] ...  Default is None [If None you should
                        provide spectrum as input file
  -y YPIX, --ypix=YPIX   Y coordinate (DEC) of the source position [sky/RAW
                        pixels] ...  Default is None [If None you should
                        provide spectrum as input file
  -o OUTSTEM, --outstem=OUTSTEM
                         Stem for output ARF File...
  --coordtyp=COORDTYP   Input Coordinate Type, .i.e., Detector Coordinates or
                        Sky Coordinates...Default is Sky   The input options
                        are ['Sky', 'sky', 'SKY', 's', 'S', 1, True] and/or
                        ['RAW', 'raw', 'Raw', 'R', 'r', 0]   Default is 'RAW'
  --vigcorrflag=VIGCORRFLAG
                        The flag for Vignetting Correction for off-axis SXT
                        observations, if needed...  Accepted options are [1.,
                        True, 'yes', 'YES', 'y', 'Y']   Default is 'no'
  --pltflag=PLTFLAG     The flag used for making the ARF diagnostics plot to
                        display various versions...  Accepted options are [1.,
                        True, 'yes', 'YES', 'y', 'Y']   Default is 'no'
  --onlyceoff=ONLYCEOFF
                        The flag to print only the coefficients for the input
                        radii..No ARF files will be generated.....  Accepted
                        options are [1., True, 'yes', 'YES', 'y', 'Y']
                        Default is 'no'

'''

def print_preamble(scriptname = ''):
    preambtxt = "=========================================================================\n"
    preambtxt += "\t Running Task: {} \n \t \t Version: 0.03 Release Date: 2019-07-03 \n \t \t Updated on: 2022-01-05\n".format(scriptname)
    preambtxt += "\t \t Developed By : Sunil Chandra, \n \t \t South African Astronomical Observatory, Cape Town \n \t \t CSR, North-West University, South Africa \n \t \t Last Updated on 05, January, 2022 \n"
    preambtxt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    preambtxt += "The estimations provided by this script is based on actual SXT observations. The profiles used here are corrected for background counts. This script is best suited for point sources. You may expect ~5 % uncertainity in flux estimations while extracting ARFs for sources with fluxes ~ 1 Crab or brighter... \n The parametric way of estimating the factor uses  kings's profile as functional form [with best fit archived parameters].\n"
    
    preambtxt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    preambtxt += "=========================================================================\n"
    print (preambtxt)



def GetIntegKingsProf(p, x) :
    #----------------------------------------------------------------------
    try :
        from scipy.integrate import quad, dblquad
    except :
        print ("The module \'kapteyn\' not found .. Installing!!")
        cmd = "sudo pip install kapteyn"
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal", -retcode,  sys.stderr)
        else:
            print("Child returned", retcode,  sys.stderr)
        from scipy.integrate import quad, dblquad
    #----------------------------------------------------------------------
    A, x0, r = p

    def KingsProf(x, A, x0, r):
        return A * (1 + (x/x0)**2)**-r
    xlim0 = 0
    val2ret = ([])
    for xlim in x :
        val2ret.append([xlim, quad(KingsProf, xlim0, xlim, args = (A, x0, r) )[0]] )
    val2ret = np.array(val2ret)    
    return val2ret[:,1]


def GetCoeff4romProfile(radius = 15, mode = 'A', datamatrix = None, parammatrix = None, strout = "") :

    #Using Table tempelate...
    if mode.upper() in ['A', 'TBL', 'TAB', 'T', 'TABLE'] :
        #print (radius, type(radius), isinstance(radius, np.ndarray))
        tbl2use = datamatrix

        xp, fp = np.array(datamatrix[0], float), np.array(datamatrix[1], float)
        #print (xp, np.interp(6.544, xp, fp))

        if isinstance(radius, np.ndarray) :
            finalcoeff = []
            for rad2use in radius :
                 coeff2prnt = np.interp(rad2use, xp, fp)
                 strout += "The Norm. EEF Coeff. for radius {}' = {} \n".format(rad2use, np.round(coeff2prnt, 3))
                 finalcoeff.append(coeff2prnt)
            finalcoeff = np.array(finalcoeff, float)
        else :
            coeff2prnt = np.interp(radius, xp, fp, left = None, right = None)
            strout += "The Norm. EEF Coeff. for radius {}' = {} \n".format(radius, np.round(coeff2prnt, 3))
            finalcoeff = float(Coeff2prnt[0])
        tbl4confreg = tbl2use

        strout += " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        strout += " The radius for 60% containment of photons  : {} arcmin \n".format(np.round(np.interp(0.60, tbl4confreg[1], tbl4confreg[0]),2))
        strout += " The radius for 90% containment of photons  : {} arcmin \n".format(np.round(np.interp(0.90, tbl4confreg[1], tbl4confreg[0]),2))
        strout += " The radius for 95% containment of photons  : {} arcmin \n".format(np.round(np.interp(0.95, tbl4confreg[1], tbl4confreg[0]),2))
        strout += " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"


    elif mode.upper() in ['B', 'PAR', 'PARAM', 'FIT', 'FITPAR'] :

        #Using Best fit parameters       
        param2fit = parammatrix[0, ]
        if isinstance(radius, np.ndarray) :
            finalcoeff = []
            for rad2use in radius :
                Coeff2prnt = GetIntegKingsProf(param2fit, [rad2use])
                strout += "The Norm. EEF Coeff. for radius {}' = {} \n".format(rad2use, np.round(Coeff2prnt[0], 3))
                finalcoeff.append(Coeff2prnt[0])
            finalcoeff = np.array(finalcoeff, float)
        else :
            Coeff2prnt = GetIntegKingsProf(param2fit, [radius])
            strout += "The Norm. EEF Coeff. for radius {}' = {} \n".format(radius, np.round(Coeff2prnt[0], 3))
            finalcoeff = float(Coeff2prnt[0])
        
        xconfreg = np.arange(0, 21, 0.25)
        tbl4confreg = [xconfreg, GetIntegKingsProf(param2fit, xconfreg)]
        #print tbl4confreg[1]
        strout += " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        strout += " The radius for 60% containment of photons  : {} arcmin \n".format(np.round(np.interp(0.60, tbl4confreg[1], tbl4confreg[0]),2))
        strout += " The radius for 90% containment of photons  : {} arcmin \n".format(np.round(np.interp(0.90, tbl4confreg[1], tbl4confreg[0]),2))
        strout += " The radius for 95% containment of photons  : {} arcmin \n".format(np.round(np.interp(0.95, tbl4confreg[1], tbl4confreg[0]),2))
        strout += " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        
    else :
        print ('No idea about the mode you wish to use for parameterization...quitting!!!!!')
        os._exit(0)
        strout += "Failed to Get Your results!!! \n Review your inputs and rerun this script\n Better Luck Next Time!!!"

    return strout, finalcoeff


def Str4romRADvec(Radvec=[15], shape='C', outstem = None, delimiter = "to"):
    #a string based on source radii to be appended in output stem..
    trp="{}".format(shape)
    for k in range(len(Radvec)):
        Rad = 4.122*Radvec[k]/60.
        #trp += "{}to".format(Rad)
        if k != len(Radvec)-1:
            trp += "{}{}".format(str(round(Rad, 1)).replace('.','p'), delimiter)
        else:
            trp += "{}".format(str(round(Rad, 1)).replace('.','p'))            
    radii4text = trp         
    #Output ARF file name..
    if outstem == None:
        outstem = "Null"           
    finarfflname = "{}_{}.arf".format(outstem, radii4text)
    return finarfflname


def MainRunOverShapes(phafile = None, outstem = None, str4logGetCoeff = "", refarf = None, evtfile = None, mode = 'b', 
                          datamatrix = None, parammatrix=None, vigcorrflag = True, coordtyp = 1):

    from astropy.io import fits

    fpha = fits.open(phafile)

    if len(fpha) >= 4 :

        try :
            regextn = fpha['REG00101'].data
        except :
            regextn = fpha[3].data

        pixscale = 4.122/60.
        
        shape, coordinate, RAD = regextn['SHAPE'][0], (regextn['X'][0][0], regextn['Y'][0][0]), regextn['R'][0]
        rotang, component = regextn['ROTANG'][0], regextn['COMPONENT'][0]
        
        txt2prnt = "++++++++ Region Info ++++++++++\n"
        txt2prnt += "SHAPE : {}\n".format(shape)
        txt2prnt += "Centered at SKY :{}\n".format(coordinate)
        txt2prnt += "Of Radius (in pix): {}\n".format(RAD)
        
        #print (len(RAD))
        if len(rotang) > 0:
            txt2prnt += "Rotation : {}\n".format(rotang)
        if len(component) > 0:
            txt2prnt += "Component : {}\n".format(component)
        txt2prnt += "++++++++++++++++++++++++++++++\n"
        print(txt2prnt)
        convfacpix2arcm = 4.122/60
        
        print ("outstem : {}".format(outstem))
        if shape == 'CIRCLE' :
            
            shape2use = 'C'
            
            radinarcm = convfacpix2arcm * RAD[0]
            
            #print(radinarcm)
            str4logGetCoeff, Coeff4corr = GetCoeff4romProfile(radius = radinarcm, mode = mode, datamatrix = datamatrix, 
                                                                  parammatrix = parammatrix, strout = str4logGetCoeff)
            
            finStr4prod = Str4romRADvec(Radvec = RAD, shape = shape2use, delimiter = "to", outstem = outstem)
            
            file2ret, str4logGetCoeff = ARFMakerUsingREFSXTARFv2(refarf = refarf, coord = coordinate, outstem = finStr4prod, 
                                                                     vigcorrflag = vigcorrflag, mode = mode, datamatrix = datamatrix,
                                                                     parammatrix = parammatrix, evtfile = evtfile, coordtyp = coordtyp, 
                                                                     Coeff4corr = Coeff4corr, str4logGetCoeff = str4logGetCoeff)


        elif shape == 'ANNULUS' :

            shape2use = 'A'
            
            radinarcm1 = convfacpix2arcm * RAD[0]
            radinarcm2 = convfacpix2arcm * RAD[-1]
            str4logGetCoeff1, Coeff4corr1 = GetCoeff4romProfile(radius = radinarcm1, mode = mode, datamatrix = datamatrix, 
                                                                    parammatrix = parammatrix, strout = str4logGetCoeff)
            
            str4logGetCoeff2, Coeff4corr2 = GetCoeff4romProfile(radius = radinarcm2, mode = mode, datamatrix = datamatrix, 
                                                                    parammatrix = parammatrix, strout = str4logGetCoeff)
            fincoeff = Coeff4corr2 - Coeff4corr1
            #print ([radinarcm1, radinarcm2, Coeff4corr2, Coeff4corr1, fincoeff])
            finStr4prod = Str4romRADvec(Radvec = RAD, shape = shape2use, delimiter = "to", outstem = outstem)
            
            file2ret, str4logGetCoeff = ARFMakerUsingREFSXTARFv2(refarf = refarf, coord = coordinate, outstem = finStr4prod, 
                                                                     vigcorrflag = vigcorrflag, mode = mode, datamatrix = datamatrix,
                                                                     parammatrix = parammatrix, evtfile = evtfile, coordtyp = coordtyp,
                                                                     Coeff4corr = fincoeff, str4logGetCoeff = str4logGetCoeff)


        elif shape == 'BOX' :

            shape2use = 'B'
            
            eqarea = RAD[0] * RAD[-1]
            
            eqrad = np.sqrt(eqarea/np.pi)
            
            eqrad = eqrad * convfacpix2arcm
            
            str4logGetCoeff, Coeff4corr = GetCoeff4romProfile(radius = eqrad, mode = mode, datamatrix = datamatrix, 
                                                                  parammatrix = parammatrix, strout = str4logGetCoeff)
            
            finStr4prod = Str4romRADvec(Radvec = RAD, shape = shape2use, delimiter = "N", outstem = outstem)
            
            file2ret, str4logGetCoeff = ARFMakerUsingREFSXTARFv2(refarf = refarf, coord = coordinate, outstem = finStr4prod, 
                                                                     vigcorrflag = vigcorrflag, mode = mode, datamatrix = datamatrix,
                                                                     parammatrix = parammatrix, evtfile = evtfile, coordtyp = coordtyp,
                                                                     Coeff4corr = Coeff4corr, str4logGetCoeff = str4logGetCoeff)


            
        elif shape == 'ELLIPSE':

            shape2use = 'E'
            
            eqarea = RAD[0] * RAD[-1]
            
            eqrad = convfacpix2arcm * np.sqrt(eqarea)
            
            str4logGetCoeff, Coeff4corr = GetCoeff4romProfile(radius = eqrad, mode = mode, datamatrix = datamatrix, parammatrix = parammatrix, strout = str4logGetCoeff)
            
            finStr4prod = Str4romRADvec(Radvec = RAD, shape = shape2use, delimiter = "N", outstem = outstem)
            
            file2ret, str4logGetCoeff = ARFMakerUsingREFSXTARFv2(refarf = refarf, coord = coordinate, outstem = finStr4prod, 
                                                                     vigcorrflag = vigcorrflag, mode = mode, datamatrix = datamatrix,
                                                                     parammatrix = parammatrix, evtfile = evtfile, coordtyp = coordtyp,
                                                                     Coeff4corr = Coeff4corr, str4logGetCoeff = str4logGetCoeff)
            
        elif shape == 'BOXANNULUS':
            
            shape2use = 'BA'
            
            eqarea1 = RAD[0] * RAD[1]
            eqrad1 = convfacpix2arcm * np.sqrt(eqarea1/np.pi)
            
            eqarea2 = RAD[2] * RAD[3]
            eqrad2 = convfacpix2arcm * np.sqrt(eqarea2/np.pi)


            str4logGetCoeff1, Coeff4corr1 = GetCoeff4romProfile(radius = eqrad1, mode = mode, datamatrix = datamatrix, 
                                                                    parammatrix = parammatrix, strout = str4logGetCoeff)
            
            str4logGetCoeff2, Coeff4corr2 = GetCoeff4romProfile(radius = eqrad2, mode = mode, datamatrix = datamatrix, 
                                                                    parammatrix = parammatrix, strout = str4logGetCoeff)
            fincoeff = Coeff4corr2 - Coeff4corr1
            
            finStr4prod = Str4romRADvec(Radvec = RAD, shape = shape2use, delimiter = "to", outstem = outstem)
            
            file2ret, str4logGetCoeff = ARFMakerUsingREFSXTARFv2(refarf = refarf, coord = coordinate,  outstem = finStr4prod,
                                                                     vigcorrflag = vigcorrflag, mode = mode, datamatrix = datamatrix,
                                                                     parammatrix = parammatrix, evtfile = evtfile, coordtyp = coordtyp,
                                                                     Coeff4corr = fincoeff, str4logGetCoeff = str4logGetCoeff)


        else :

            print ("The unsupported shape for the region found in PHA file....")
            file2ret, str4logGetCoeff = None, None
    else :

        print ("Input pha file has no REGION info so can't fetch details...Enter manually !!!")
        file2ret, str4logGetCoeff = None, None

    return  file2ret, str4logGetCoeff   



def ARFMakerUsingREFSXTARFv2(refarf = None, coord = (0, 0), outstem = 'OutARF_v00_auto', vigcorrflag = True, mode = 'a', 
                                 datamatrix = 'datamatrix', parammatrix = 'parammatrix', evtfile = None, coordtyp = 0, 
                                 Coeff4corr = 0.5, str4logGetCoeff = 'log') :

    #print ("outStem: {} ".format(outstem))
    def SrcRadiiCorr(refarf = None, coeff2prnt = 1.0, outstem = outstem) :
        
        outstem = os.path.splitext(outstem)[0]
        #Reading Reference ARF
        refsxtarf = fits.open(refarf)
        #Copy ref in a proxy ARF
        outarf2wrt = refsxtarf
        #Modify ARF based on estimated coeff. for a particular source region..
        outarf2wrt[1].data['SPECRESP'] =  coeff2prnt * outarf2wrt[1].data['SPECRESP']
        finarfflname = "{}.arf".format(outstem)   
        #Writing the modified ARF 
        outarf2wrt.writeto(finarfflname, overwrite=True)
        outarf2wrt.close()
        refsxtarf.close()
        #print ("outStem: {} \n outfl = {}".format(outstem, finarfflname))

        if os.path.exists(finarfflname) :
            outfl = finarfflname
        else :
            outfl = None
        return outfl


    def ARFOffAngleCorr(inparf = None, offangle = None, outarf=None) :

        if ( inparf != None ) and (offangle != None) :
            print ("ARF Correction for Off-Axis Angle : \"yes\" \n Input ARF File Name : {} \n Input Off-Axis Angle : {}".format(inparf, offangle))
            if not os.path.exists(inparf) :

                inparf = input("Enter The Full Path of ARF File distributed by SXT Team as you had either provided it wrongly or somehow deleted from the path you enterd : ")

                if not inparf :
                    os._exit(0)

            arffl = fits.open(inparf)
            #Important Vignetting Parameters... 
            P0 = -0.00052850083
            P1 = 0.99962130
            P2 = 0.0020385522

            E = (arffl[1].data['ENERG_LO'] +   arffl[1].data['ENERG_HI'])/2.
            C = P0 * (P1**E) + P2
            V = 1 - (C*(offangle**2))

            arffl[1].data['SPECRESP'] = arffl[1].data['SPECRESP']*V
            if outarf == None :
                outarfstem = os.path.splitext(os.path.basename(inparf))[0]
                outarf = "{}_ArfVigCorr.arf".format(outarfstem)
            arffl.writeto(outarf, overwrite=True)
            #arffl.close()
            return  outarf

    
    #Write a new ARF based on coeff...
    outarf = SrcRadiiCorr(refarf = refarf, coeff2prnt = Coeff4corr, outstem = outstem)
    print(str4logGetCoeff)
    if outarf != None :
        str4logGetCoeff += "ARF File {} written successfully...\n".format(outarf)
        file2ret = outarf
        if vigcorrflag :
            #Correcting Newly generated ARF for vignetting effect...
            #Bore-sight Information..
            X0, Y0 = 301.8, 287.1
            X, Y = coord[0], coord[1]
            if coordtyp == 0 :
                offangle = np.sqrt( (X - X0)**2 + (Y - Y0)**2 )*(4.122/60.)
            else :
                if evtfile == None :
                    evtfile = input("Enter Events file to go ahead : ")
                # Extracting 
                X0, Y0 = GetCentralcoordLookup(evtfile, rawx=X0, rawy = Y0, coordtol=2e0, stepcoordtol=5e-1)
                offangle = np.sqrt( (X - X0)**2 + (Y - Y0)**2 )*(4.122/60.)
            #ARF after Vig. Correction
            outfilevigcorr = "{}_VigCorr.arf".format(os.path.splitext(outarf)[0])
            outarfvigcorr = ARFOffAngleCorr(inparf = outarf, offangle = offangle, outarf = outfilevigcorr)
            if os.path.exists(outarfvigcorr) :
                str4logGetCoeff += "Resulting ARF File {} is corrected for Vig. Effect....\n".format(outarfvigcorr)
                file2ret = outarfvigcorr
    else :
        str4logGetCoeff += "ARF File Generation failed....Please recheck your inputs! \n"
        file2ret = None
    return file2ret, str4logGetCoeff



def main() :

    import optparse

    usage = "usage: %prog [options] "
    parser = optparse.OptionParser(usage)

    parser.add_option("-m", "--mode", dest = "mode", help = "Mode for extraction; \n Two options a) using table and b) using fit parameters \n Default is b) \n Input options are [\'a\', \'tbl\', \'tab\', \'t\', \'table\'] for a) and one of [\'b\', \'par\', \'param\', \'fit\', \'fitpar\'] for b) ", default = 'b')

    parser.add_option("-r", "--radius", dest = "radius", help = "Source Region Size in arcmin.. Always use for cicular regions...\n Multiple radii (r1, r2, r3) can be entered in format or \'r1nr2nr3\' or \'r1ANDr2ANDr3\' or \'r1ORr2ORr3\'\n Default is 15 arcm", default = 15)

    parser.add_option("", "--sxtarf", dest = "sxtarf", help = "The Full Path of Base ARF provided by SXT Team...\n  For Example \'sxt_pc_excl00_v04.arf\'", default = None)

    parser.add_option("-e", "--evtfile", dest = "evtfile", help = "The Events file in case you are entering source postion in sky pixels instead of detector coordinates.....Default = None", default = None)

    parser.add_option("", "--sxtpha", dest = "sxtpha", help = "The Name of the input SXT spectrum [needed to fetch source region related information]...\n  For Example \'sxt_pc_mytarget_src.pha\'\n Default is None..[This means source coordinates should be entered manually]", default = None)

    parser.add_option("-x", "--xpix", dest = "xpix", help = " X coordinate (RA) of the source position [sky/RAW pixels] ...\n Default is None [If None you should provide spectrum as input file", default = None)

    parser.add_option("-y", "--ypix", dest = "ypix", help = " Y coordinate (DEC) of the source position [sky/RAW pixels] ...\n Default is None [If None you should provide spectrum as input file", default = None)

    parser.add_option("-o", "--outstem", dest = "outstem", help = " Stem for output ARF File...", default = "OutARF_v04")

    parser.add_option("", "--coordtyp", dest = "coordtyp", help = "Input Coordinate Type, .i.e., Detector Coordinates or Sky Coordinates...Default is Sky \n The input options are ['Sky', 'sky', 'SKY', 's', 'S', 1, True] and/or ['RAW', 'raw', 'Raw', 'R', 'r', 0] \n Default is 'RAW' ", default = 'sky')

    parser.add_option("", "--vigcorrflag", dest = "vigcorrflag", help = "The flag for Vignetting Correction for off-axis SXT observations, if needed...\n Accepted options are [1., True, 'yes', 'YES', 'y', 'Y'] \n Default is 'no'", default = 'no')

    parser.add_option("", "--pltflag", dest = "pltflag", help = "The flag used for making the ARF diagnostics plot to display various versions...\n Accepted options are [1., True, 'yes', 'YES', 'y', 'Y'] \n Default is 'no'", default = 1)


    parser.add_option("", "--coeffonly", dest = "coeffonly", help = "The flag to print only the coefficients for the input radii..No ARF files will be generated.....\n Accepted options are [1., True, 'yes', 'YES', 'y', 'Y'] \n Default is 'no'", default = 'no')


    (options, args) = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    #global mode, evtfile, refarf, vigcorrflag, sxtpha, coordtyp, pltflag, coeffonly, str4logGetCoeff
     
    mode = options.mode 
    evtfile = options.evtfile
    refarf = options.sxtarf
    vigcorrflag  = options.vigcorrflag
    sxtpha = options.sxtpha
    coordtyp = options.coordtyp
    pltflag = options.pltflag
    radius = options.radius
    srcx, srcy = options.xpix, options.ypix 
    coeffonly, outstem = options.coeffonly, options.outstem
    
    print (vigcorrflag)
    #Vignetting Correction flag input ..
    if vigcorrflag in ['YES', 'Yes', 'Y', 'y', 'yes']:
        vigcorrflag = True 
    else :
        vigcorrflag = False

    #Plot flag input ..
    if (pltflag in ['YES', 'Yes', 'Y', 'y', 'yes']) or (pltflag == 1) or (pltflag == True):
        pltflag = True 
    else :
        pltflag = False

    #mode flag input ..
    if mode.lower() in ['b', 'par', 'param', 'fit', 'fitpar']:
        mode = 'b' 
    else :
        mode = 'a'

    #Plot flag input ..
    if (pltflag in ['YES', 'Yes', 'Y', 'y', 'yes']) or (pltflag == 1) or (pltflag == True):
        pltflag = True 
    else :
        pltflag = False

    #coordtyp flag input ..
    if (coordtyp in ['RAW', 'raw', 'Raw', 'R', 'r']) or (coordtyp == 0):
        coordtyp = 0 
    else :
        coordtyp = 1

    #Coeffonly flag input ..
    if (coeffonly in ['NO', 'No', 'N', 'n', 'no']) or (coeffonly == 0) or (coeffonly == False):
           
        #Extracting Source region related information...
        if sxtpha != None :
            str4logGetCoeff = "Input SXT-Spectrum : {} \t Used for source region info..\n".format(sxtpha)
            if os.path.exists(sxtpha):
                file2ret, str4logGetCoeff = MainRunOverShapes(phafile = sxtpha, outstem = outstem, str4logGetCoeff = "",
                                                              refarf = refarf, evtfile = evtfile, mode = mode, 
                                                              datamatrix = datamatrix, parammatrix=parammatrix,
                                                              vigcorrflag = vigcorrflag, coordtyp= coordtyp)
                #MainRunOverShapes(phafile = sxtpha, outstem = outstem)
                if pltflag :
                    try :
                        from datetime import datetime
                        datetimenow = datetime.now().strftime('%Y-%m-%d %H:%M:%S').replace('-','').replace(' ','T').replace(":","h").split('h')[0]+'h'
                    except :
                        datetimenow = 'temp_now'
                    outdigpltname = "{}_ARFCompPlt_{}".format(outstem, datetimenow)
                    arflist2plt = [[refarf,'orig'], [file2ret, 'new']]
                    ARFListPlotter(arflist = arflist2plt, figsize = [12,9], outdiagplt = outdigpltname)

            else:
                sxtpha = input("The pha file parsed from command line is not found ... Please enter here with full path : ")
                if not sxtpha:
                    os._exit(0)
                file2ret, str4logGetCoeff = MainRunOverShapes(phafile = sxtpha, outstem = outstem)
                
            print(str4logGetCoeff)
        else :
            print ("coeffonly = False with Insufficient information!! \n QUITING!! ")
            os._exit(999)
    else : #when coeffonly is True 
        radiivec = clarifyradius(radius)
        str4logGetCoeff, finalcoeff = GetCoeff4romProfile(radius = radiivec, mode = mode, datamatrix = datamatrix, parammatrix = parammatrix, strout = "")
        print(str4logGetCoeff)


def clarifyradius(radius) :

    def radius2radiivec(radius) :
        # Using Input Radius/Radii
        #--------------------------------------------------------
        if ( len( str(radius).split('n') ) > 1) :
            radiivec = np.array(str(radius).split('n'), float)
        elif (len( str(radius).split('AND') ) > 1) :
            radiivec = np.array(str(radius).split('AND'), float)
        elif (len( str(radius).split('OR') ) > 1) :
            radiivec = np.array(str(radius).split('OR'), float)
        else :
            radiivec = np.array([radius], float)
        #--------------------------------------------------------
        return radiivec

    #Use it as per need...
    if radius != None :
        radiivec = radius2radiivec(radius)
    else :
        radius = input("You had not provided input spectrum and also no information about source radius was entered...\n Please Enter radius of circular source region or equivalent radii [in arcm] : ")
        if not radius :
            os._exit(0)
        radiivec = radius2radiivec(radius)
    return radiivec


def clarifysrccoordinate(srcx, srcy) :

    if (srcx != None) and (srcy != None) :
        srccoord = (srcx, srcy)
    else :     
        srcx = input("The input spectrum is not found and also the source position = None... Please Enter X  coordinate of useful source position : ")
        srcx = input(" Please Enter Y  coordinate of useful source position : ")
        srccoord = (srcx, srcy)
    return srccoord


def GetCentralcoordLookup(evtfile, rawx=300.0, rawy = 300.0, coordtol=2e0, stepcoordtol=5e-1) :

    import numpy as np
    from astropy.io import fits
    from astropy.table import Table

    evtfl = fits.open(evtfile)
    evthdu = evtfl[1].data
    evtfl.close()

    evttbl = Table(evthdu)
    N = 0

    while N <= 3 :
        xrawlowlim, xrawuplim, yrawlowlim, yrawuplim = rawx - coordtol, rawx + coordtol, rawy - coordtol, rawy + coordtol

        fltrdhdu = evttbl[(evttbl['RAWX']>=xrawlowlim) & (evttbl['RAWX']<=xrawuplim) & (evttbl['RAWY']>=yrawlowlim) & (evttbl['RAWY']<=yrawuplim)]

        coordtol += stepcoordtol

        N = len(fltrdhdu)

    fltrdcoordtbl = Table()
    fltrdcoordtbl.add_column(fltrdhdu['RAWX'], index = 0)
    fltrdcoordtbl.add_column(fltrdhdu['RAWY'], index = 1)
    fltrdcoordtbl.add_column(fltrdhdu['X'], index = 2)
    fltrdcoordtbl.add_column(fltrdhdu['Y'], index = 3)
    #print ( fltrdcoordtbl )

    #print ( (np.mean(fltrdhdu['RAWX']), np.mean(fltrdhdu['RAWY']) ), (np.mean(fltrdhdu['X']), np.mean(fltrdhdu['Y'])) )

    return [np.mean(fltrdhdu['X']), np.mean(fltrdhdu['Y'])]



def ARFListPlotter(arflist=[], figsize=[12,9], outdiagplt = 'outARfComp.png') :

    linevec = ['-', '-.', ':', '--', '-', '-.', ':', '--', '-', '-.', ':', '--']
    colorvec = ['k', 'k', 'k', 'k', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'g', 'g', 'g', 'g', 'c', 'c', 'c', 'c', 'm', 'm', 'm', 'm', 'y' ]
    
    if len(arflist) > 0 :

        plt.figure(figsize=figsize)
        for i in range(len(arflist)) :
            print (arflist[i][0])
            arffile, arfflag = arflist[i][0], arflist[i][1]

            arfhdu = fits.open(arffile)

            plt.plot((arfhdu[1].data['ENERG_LO']+ arfhdu[1].data['ENERG_HI'])/2., arfhdu[1].data['SPECRESP'], lw=1.0, linestyle = linevec[i], color = colorvec[i], label = arfflag)

        plt.legend(loc='best', fontsize=13)
        plt.xscale('log')
        #plt.yscale('log')
        plt.minorticks_on()
        plt.xlabel("Energy [keV]", fontsize = 19)
        plt.ylabel("SPECRESP [cm$^2$]", fontsize = 19)

        plt.title('Comparison of Input and Resulting ARFs ', fontsize = 17)
        plt.xlim([0.25, 10.05])
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.tick_params(axis='both', which='minor', labelsize=13)

        plt.savefig(outdiagplt)
        #plt.show()



#Important data for profiles ...based on experiments with various observations from different classes... Please do not change this part
#------------------------------------------------------------------------------------------------------------------------------------------------------

datamatrix = [[ 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3.0, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4.0, 4.125, 4.25, 4.375, 4.5, 4.625, 4.75, 4.875, 5.0, 5.125, 5.25, 5.375, 5.5, 5.625, 5.75, 5.875, 6.0, 6.125, 6.25, 6.375, 6.5, 6.625, 6.75, 6.875, 7.0, 7.125, 7.25, 7.375, 7.5, 7.625, 7.75, 7.875, 8.0, 8.125, 8.25, 8.375, 8.5, 8.625, 8.75, 8.875, 9.0, 9.125, 9.25, 9.375, 9.5, 9.625, 9.75, 9.875, 10.0, 10.125, 10.25, 10.375, 10.5, 10.625, 10.75, 10.875, 11.0, 11.125, 11.25, 11.375, 11.5, 11.625, 11.75, 11.875, 12.0, 12.125, 12.25, 12.375, 12.5, 12.625, 12.75, 12.875, 13.0, 13.125, 13.25, 13.375, 13.5, 13.625, 13.75, 13.875, 14.0, 14.125, 14.25, 14.375, 14.5, 14.625, 14.75, 14.875, 15.0, 15.125, 15.25, 15.375, 15.5, 15.625, 15.75, 15.875, 16.0, 16.125, 16.25, 16.375, 16.5, 16.625, 16.75, 16.875, 17.0, 17.125, 17.25, 17.375, 17.5, 17.625, 17.75, 17.875, 18.0, 18.125, 18.25, 18.375, 18.5, 18.625, 18.75, 18.875, 19.0, 19.125, 19.25, 19.375, 19.5, 19.625, 19.75, 19.875], 
[ -0.0278 , -0.0139, -0.0, 0.0138, 0.0277, 0.0416, 0.0555, 0.0693, 0.0832, 0.0972, 0.1118, 0.1263, 0.1409, 0.1555, 0.1698, 0.184, 0.1983, 0.2125, 0.2266, 0.2405, 0.2545, 0.2684, 0.2823, 0.2955, 0.3088, 0.322, 0.3353, 0.3483, 0.3612, 0.374, 0.3869, 0.3997, 0.4119, 0.4242, 0.4365, 0.4488, 0.4604, 0.472, 0.4835, 0.495, 0.5063, 0.5174, 0.5284, 0.5394, 0.5505, 0.5608, 0.5711, 0.5814, 0.5917, 0.6015, 0.6112, 0.6208, 0.6305, 0.6399, 0.6488, 0.6576, 0.6665, 0.6754, 0.6836, 0.6917, 0.6999, 0.708, 0.7158, 0.7233, 0.7308, 0.7383, 0.7457, 0.7525, 0.7593, 0.7661, 0.7729, 0.7793, 0.7855, 0.7917, 0.7979, 0.8039, 0.8095, 0.8151, 0.8207, 0.8263, 0.8315, 0.8366, 0.8417, 0.8468, 0.8516, 0.856, 0.8605, 0.8649, 0.8693, 0.8733, 0.8773, 0.8813, 0.8852, 0.8889, 0.8925, 0.8961, 0.8997, 0.9032, 0.9063, 0.9094, 0.9126, 0.9157, 0.9185, 0.9212, 0.9239, 0.9266, 0.9292, 0.9316, 0.934, 0.9364, 0.9388, 0.941, 0.9432, 0.9453, 0.9475, 0.9495, 0.9513, 0.9532, 0.955, 0.9568, 0.9586, 0.9603, 0.9621, 0.9638, 0.9654, 0.9669, 0.9684, 0.9699, 0.9713, 0.9726, 0.9739, 0.9752, 0.9765, 0.9777, 0.9789, 0.9801, 0.9812, 0.9823, 0.9834, 0.9845, 0.9855, 0.9865, 0.9875, 0.9884, 0.9893, 0.9903, 0.9911, 0.9919, 0.9928, 0.9936, 0.9943, 0.995, 0.9957, 0.9964, 0.997, 0.9977, 0.9983, 0.9989, 0.9995, 1.0002 ] ]
 

parammatrix = np.array([[0.10973051034777614, 33.81373398681025, 11.300603962874558], [ 0.11137994, 26.28498828,  7.27950504]], float)

#------------------------------------------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    import sys

    sys.settrace

    scriptname = sys.argv[0]

    print_preamble(scriptname = scriptname)

    sys.setrecursionlimit(500000)

    main()

