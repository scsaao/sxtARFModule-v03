#!/usr/bin/env python

#------------------------------------ Readme Content ----------------------------------------------------------------------
'''
chandra@blazars:~$ python3 ~/auxpyscrpt/make_sxt_EEF_v02_py2.py -h
=========================================================================
	  Running 'SXT EEF Maker' Tool 
 	  Task: /home/chandra/auxpyscrpt/make_sxt_EEF_v02_py2.py 
 	  Version: 02 Release Date: 2018-11-14 
	  Originally Developed By : 

 	 	 Sunil Chandra, 
 	 TIFR, Mumbai, 15 January, 2017 
 	 Updated on 08 July, 2019 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Description: This tool makes Encirled Energy Fraction (EEF) profile for sources by using merged clean events files as input. The profile is shown to have radii containing 90% and 60% of total incoming photons. This version provides an option to correct profiles for background [Default].
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
=========================================================================

Usage: make_sxt_EEF_v02_py2.py [options] 

Options:
  -h, --help            show this help message and exit
  -e EVTFILE, --evtfile=EVTFILE
                        Input Events file Name;   The most important input of
                        script give with full path if not in same directory
  -r RA, --ra=RA        Input RA in Sky pixels   Note that open your events
                        file with ds9 and note the position of centriod for
                        the source in pixels
  -d DEC, --dec=DEC     Input DEC in Sky pixels
  -o OUTFILE, --outfile=OUTFILE
                        Stem for output files
  -x MINR, --min=MINR   Radius of inner circle for annulii (in Sky pixel);
                        Default value is 8 pix (~0.5 arcm)
  -y MAXR, --max=MAXR   Radius of outer circle for annulii (in Sky pixel);
                        Default value is 328 pix (~22.5 arcm)
  -s STEPS, --steps=STEPS
                        Steps size to decide number of data points in outpur
                        profile [N = int((r_out - r_in)/stepsize)];  default
                        is 8 pixel
  -p PROCESS, --process=PROCESS
                        Plot type: 'plot' or 'both'; Whether you wish to make
                        plot only or make tables as well; You may use any of
                        ["BOTH", "AUTO", "B", "A", "PROCESS"] for default way
                        [making profile and plot] and any of ["PLOT", "PL",
                        "P"] for plot-only option if you already have tables
                        with you
  --estring=ESTRING     String for energy selection; the acceptable input
                        format is 0p3to7p0 (for 0.3-7.0 keV)   Default is
                        0.3-7.0 keV
  -g GRADE, --grade=GRADE
                        String for grade selections in the format 0-12
                        Default is all grades 0-12
  -t SIMEEF, --simeef=SIMEEF
                        Filename with full path for EEF from simulation;
                        Needed only for POC releted bussiness, you make ignore
                        it, Default is None, which mean no
  -l IPACFILE, --ipacfile=IPACFILE
                        IPAC file name with data; needed for 'plot' option for
                        process, Default = None
  --coordtyp=COORDTYP   Input Coordinate Type, .i.e., Detector Coordinates or
                        Sky Coordinates...Default is Sky   The input options
                        are ['Sky', 'sky', 'SKY', 's', 'S', 1, True] and/or
                        ['RAW', 'raw', 'Raw', 'R', 'r', 0]   Default is 'RAW'
  --bkgcorr=BKGCORR     The flag for back-ground correction, Default options
                        is 'yes'
  --NrefStyl=NREFSTYL   The flag for using method to normalize the profile,
                        options are either use total counts from largest safe
                        radius of source extraction or the counts from a
                        source region of 20 arcm [extra-polated, in cases
                        20arcm is not possible] to use.   Input options =
                        ['extrapolate', 'ff', 'full', '20arcm', 'auto', 1] for
                        using 20arcm anything else except ['both', 'BOTH',
                        'debug', 2] will result in using largest possible safe
                        radius for normalization of profiles; You may also
                        chose for any of ['kings', 'fit', 'model', 0]
  --normradius=NORMRADIUS
                        The radius [in arcm] what script should use in case
                        user is confident with source PSF
  --logfile=LOGFILE     The Name of logfile ..Default is
                        'LogFile_DateTimeFormat.log'
  --usepoly3flag=USEPOLY3FLAG
                        The Flag whether polynomial of order 3 should be used
                        for fitting or not.. Default = 'no'
'''
#-------------------------------------------------------------------------------------------------------------------


def print_preamble(inpstr = 'scriptname'):

    version = (scriptname.split('_')[-2]).split('v')[-1]

    updatereleasedate = '08 July, 2019'

    creationdate = '15 January, 2017'

    preambtxt = "=========================================================================\n"
    preambtxt += "\t  Running 'SXT EEF Maker' Tool \n \t  Task: {} \n \t  Version: {} Release Date: 2018-11-14 \n".format(scriptname, version)
    preambtxt += "\t  Originally Developed By : \n\n \t \t Sunil Chandra, \n \t TIFR, Mumbai, {} \n \t Updated on {} \n".format(creationdate, updatereleasedate)
    preambtxt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n"
    preambtxt += "Description: This tool makes Encirled Energy Fraction (EEF) profile for sources by using merged clean events files as input. The profile is shown to have radii containing 90% and 60% of total incoming photons. This version provides an option to correct profiles for background [Default].\n"
    preambtxt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    preambtxt += "=========================================================================\n"
    print (preambtxt)


def FrameProductMaker(evtfile,regfile=None,CHANMIN=30,CHANMAX=701,productstem="SrcNameOut", grade_flag="0-12", curvebinsize = 120, sky_coord = 0):
    
    import os, shutil, subprocess

    if regfile == None :
        print ("The product is contaminated with corner sources")
        status = "The product is contaminated with corner sources"
        regfile = ''
    else :
        status = "The product is made after excluding the corner sources"

    xco_outstr="{}_scrpt.xco".format(productstem)
    
    #regfInp=regfile

    if len(evtfile.split("/")) > 1 :

        evtfname = os.path.split(evtfile)[-1]

        if os.path.dirname(evtfile) == os.getcwd() :
            extfile = False

        else :

            shutil.copyfile(evtfile, os.path.join("./", evtfname))
            extfile = True

        evtdir = "./"
        evtfile = evtfname
        
    else :
        evtdir="./"
        evtfile = os.path.split(evtfile)[-1]
        extfile = False

    if os.path.exists(xco_outstr) :
        os.remove(xco_outstr)

    fl_xco=open(str(xco_outstr),'a')
    
    fl_xco.write("xsel\n")
    
    fl_xco.write("read events\n")
    
    fl_xco.write(("%s\n")%(evtdir))
    
    fl_xco.write(("%s\n")%(evtfile))
    
    fl_xco.write("yes\n")

    if sky_coord == 0 :
        fl_xco.write("set xyname RAWX RAWY\n")
    
    fl_xco.write("extract all\n")
    
    fl_xco.write(("filter pha_cutoff %d %d \n")%(CHANMIN,CHANMAX))
    
    fl_xco.write(("filter grade %s\n")%(grade_flag))
    
    fl_xco.write(("filter region %s \n")%(regfile))
    
    fl_xco.write("extract image\n")
    
    fl_xco.write("save image {}.img \n".format(productstem))

    #fl_xco.write(("set binsize %f \n")%(curvebinsize))
    
    #fl_xco.write("extract curve\n")
    
    #fl_xco.write("save curve {}.lc \n".format(productstem))
    
    #fl_xco.write("clear pha_cutoff\n")
    
    #fl_xco.write("extract spectrum\n")
    
    #fl_xco.write("save spectrum {}.pha \n".format(productstem))
    
    fl_xco.write("quit\n")
    
    fl_xco.write("no\n")
    
    fl_xco.close()

    cmd = "xselect @{}".format(xco_outstr)

    try:
        retcode = subprocess.call(cmd, shell=True)

        if retcode < 0:
            print("Child was terminated by signal", -retcode,  sys.stderr)
        else:
            print("Child returned", retcode,  sys.stderr)
    except OSError as e:
        print("Execution failed:", e,  sys.stderr)

    os.remove("{}".format(xco_outstr))

    if extfile == True :
        os.remove(evtfile)
    return status 


def main() :

    import numpy as np
    import os, glob, shutil, optparse, glob 
    from astropy.io import fits 
    from astropy.table import Table, Column
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt

    usage = "usage: %prog [options] "

    parser = optparse.OptionParser(usage)

    parser.add_option("-e", "--evtfile", dest = "evtfile", help = "Input Events file Name; \n The most important input of script give with full path if not in same directory", default = None)

    parser.add_option("-r", "--ra", dest = "ra", help = "Input RA in Sky pixels \n Note that open your events file with ds9 and note the position of centriod for the source in pixels", default = None)

    parser.add_option("-d", "--dec", dest = "dec", help = "Input DEC in Sky pixels", default = None)

    parser.add_option("-o", "--outfile", dest = "outfile", help = "Stem for output files", default = "outfile.png")

    parser.add_option("-x", "--min", dest = "minr", help = "Radius of inner circle for annulii (in Sky pixel); \n Default value is 8 pix (~0.5 arcm) ", default = 8)

    parser.add_option("-y", "--max", dest = "maxr", help = "Radius of outer circle for annulii (in Sky pixel); \n Default value is 328 pix (~22.5 arcm)", default = 328)

    parser.add_option("-s", "--steps", dest = "steps", help = "Steps size to decide number of data points in outpur profile [N = int((r_out - r_in)/stepsize)];\n default is 8 pixel", default = 8)

    parser.add_option("-p", "--process", dest = "process", help = "Plot type: 'plot' or 'both'; Whether you wish to make plot only or make tables as well; You may use any of  [\"BOTH\", \"AUTO\", \"B\", \"A\", \"PROCESS\"] for default way [making profile and plot] and any of [\"PLOT\", \"PL\", \"P\"] for plot-only option if you already have tables with you", default = 'both')

    parser.add_option("", "--estring", dest = "estring", help = "String for energy selection; the acceptable input format is 0p3to7p0 (for 0.3-7.0 keV) \n Default is 0.3-7.0 keV", default = "all")

    parser.add_option("-g", "--grade", dest = "grade", help = "String for grade selections in the format 0-12 \n Default is all grades 0-12", default = "0-12")

    parser.add_option("-t", "--simeef", dest = "simeef", help = "Filename with full path for EEF from simulation; \n Needed only for POC releted bussiness, you make ignore it, Default is None, which mean no", default = None)

    parser.add_option("-l", "--ipacfile", dest = "ipacfile", help = "IPAC file name with data; needed for 'plot' option for process, Default = None", default = None)

    parser.add_option("", "--coordtyp", dest = "coordtyp", help = "Input Coordinate Type, .i.e., Detector Coordinates or Sky Coordinates...Default is Sky \n The input options are ['Sky', 'sky', 'SKY', 's', 'S', 1, True] and/or ['RAW', 'raw', 'Raw', 'R', 'r', 0] \n Default is 'RAW' ", default = 'sky')

    parser.add_option("", "--bkgcorr", dest = "bkgcorr", help = "The flag for back-ground correction, Default options is 'yes' ", default = 'yes')

    parser.add_option("", "--NrefStyl", dest ="NrefStyl", help= "The flag for using method to normalize the profile, options are either use total counts from largest safe radius of source extraction or the counts from a source region of 20 arcm [extra-polated, in cases 20arcm is not possible] to use. \n Input options = ['extrapolate', 'ff', 'full', '20arcm', 'auto', 1] for using 20arcm anything else except ['both', 'BOTH', 'debug', 2] will result in using largest possible safe radius for normalization of profiles; You may also chose for any of ['kings', 'fit', 'model', 0]", default = 0)

    parser.add_option("", "--normradius", dest = "normradius", help = "The radius [in arcm] what script should use in case user is confident with source PSF ", default = 20)

    parser.add_option("", "--logfile", dest = "logfile", help = "The Name of logfile ..Default is 'LogFile_DateTimeFormat.log'", default = "LogFile_DateTimeFormat.log")

    parser.add_option("", "--usepoly3flag", dest = "usepoly3flag", help = "The Flag whether polynomial of order 3 should be used for fitting or not.. Default = 'no'", default = "LogFile_DateTimeFormat.log")

    (options, args) = parser.parse_args()

    
    if len(sys.argv[1:]) == 0 :
        parser.print_help()
        parser.exit()
    
    evtfilename, process = options.evtfile, options.process
    srcra, srcdec, outfile = options.ra, options.dec, options.outfile
    minsizesky = int(options.minr)
    maxsizesky = int(options.maxr)
    setpsizepix = float(options.steps)
    outradivec = np.arange(minsizesky, maxsizesky, setpsizepix)
    estring, grade = options.estring, options.grade
    simeeffile = options.simeef
    ipacfilename = options.ipacfile
    coordtyp = options.coordtyp
    bkgcorrflag = options.bkgcorr 
    NrefStyl = options.NrefStyl
    normradii = options.normradius
    logfile = options.logfile
    usepoly3 = options.usepoly3flag

    logfl = open(logfile, 'a')

    str4log = "+++++++++++++++ Important Inputs +++++++++++++++++++++++\n"
    str4log += "Input Events File Name : {}\n".format(evtfilename)
    str4log += "Input Coordinate Type [RAW/SKY] : {}\n".format(coordtyp)
    str4log += "Input source coordinates (X,Y)/(RAWX, RAWY) : ({}, {})\n".format(srcra, srcdec)
    str4log += "String for Stem of output files and product directory : {}\n".format(outfile)
    str4log += "Minimum and Maximum radii of circles used to generate EEF profile (arcm) : ({}, {})\n".format(np.round(minsizesky*4.122/60.,2),np.round(maxsizesky*4.122/60.,2))
    str4log += "Stepsize in arcm used for making profile : {}\n".format(np.round(setpsizepix*4.122/60,2))
    str4log += "Input energy String for energy selection for EEF profile : {}\n".format(estring)
    str4log += "Input grade selection : {}\n".format(grade)
    str4log += "Simulated Profile [Only for POC purpose] \n".format(simeeffile)
    str4log += "Input IPAC File for profile in case you want to use plot-only option : {} \n".format(ipacfilename)
    str4log += "Input Flag for Background correction [Default = yes] : {} \n".format(bkgcorrflag)
    str4log += "Input Flag for type of normalization selection...[Default is through fitting the specific profile...] : {} \n".format(NrefStyl)
    str4log += "Input for radius of circle used for normalization of the EEF profile [Default = 20 arcm] : {} \n".format(normradii)
    str4log += "Input Flag for using polynomial of order 3  : {} \n".format(usepoly3)

    #Method of normalization of profiles..
    if (NrefStyl in ['extrapolate', 'ff', 'full', '20arcm', 'auto']) or (int(NrefStyl) == 1) :

        NrefStyl = 1

    elif (NrefStyl in ['both', 'b', 'BOTH', 'debug']) or (int(NrefStyl) == 2) :

        NrefStyl = 2

    elif (NrefStyl in ['fix', 'custom','useradii']) or (int(NrefStyl) == 3) :

        NrefStyl = 3

    else :

        NrefStyl = 0

    #Whther to use background corrections or not...
    if bkgcorrflag in [1, True, 'yes', 'Yes', 'YES', 'y', 'Y'] :
        bkgcorrflag = True

    else :
        bkgcorrflag = False

    # Whther to use sky coordinated or detector coordinates for entered position..
    if (coordtyp in ['Sky', 'sky', 'SKY', 's', 'S']) or (coordtyp == 1) or (coordtyp == True) :

        sky_coord = 1

    elif ( coordtyp in ['RAW', 'raw', 'Raw', 'R', 'r'] ) or (coordtyp == 0) or (coordtyp == False) :

        sky_coord = 0

    else :

        print ("Input Options for Coordinate Type is not understandable ... I am assuming it is SKY coordinates... !!")

        sky_coord = 1

    print (NrefStyl, sky_coord)
    #os._exit(0)
    color2plt = ['k', 'r', 'b', 'g', 'y', 'c', 'm', 'k', 'r', 'b', 'g', 'y', 'c', 'm']
    getval  =True
    rad2ret = [1.0, 2.0, 3.0]

    outputdir2prod = os.path.splitext(os.path.basename(outfile))[0]

    outfile = outputdir2prod

    outpltfl = "{}_Outplt.png".format(outfile)

    if not os.path.exists(outputdir2prod) :
        os.makedirs(outputdir2prod)


    #Making the selections for energy range to look for....

    if estring in ["ALL", "all", "All", True, "Default", "default"] :

        envec = ["0p3to7p0", "0p3to0p7", "0p7to2p0", "2p0to7p0"] 

    elif len(estring.split("_and_")) :

        envec = estring.split("_and_")

    else :

        envec = [estring]

    if evtfilename == None :

        print ("You have not parsed any events file!!! Exiting !!!")
        os._exit()

    try :
        evtf = fits.open(evtfilename)
        evthdr = evtf[0].header
        evtf.close()
        srcname = (evthdr['OBJECT'].replace(' ', ''))
        obsid = evthdr['OBSID']
        exposure = float(evthdr['EXPOSURE'])

    except :
        srcname = os.popen('gethead {} OBJECT'.format(evtfilename)).read()
        exposure = os.popen('gethead {} EXPOSURE'.format(evtfilename)).read()
        obsid = os.popen('gethead {} OBS_ID'.format(evtfilename)).read()
        #srcname = raw_input("Enter the object name mannually : ")
        #exposure = raw_input("Enter the exposure in second : ")
        #obsid = raw_input("Enter the Observation ID : ")
        dictexmp = {'object' : '{}'.format(srcname), 'exposure' : '{}'.format(exposure), 'OBS_ID' : '{}'.format(obsid)}
        evthdr = dictexmp

    str4log += "Important Header Information about observation : \n \t\t\t OBSID : {} \n \t\t\t EXPOSURE : {} ks \n \t\t\t TARGET NAME : {} \n".format(obsid, np.round(float(exposure)/1e3, 2), srcname)
    str4log += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"


    if str(process).upper() in ["BOTH", "AUTO", "B", "A", "PROCESS"] :

        if (srcra == None) and (srcdec == None) :
            srcra = float(raw_input("Enter RA in SKY X : "))
            srcdec = float(raw_input("Enter DEC in SKY Y: "))

        #plt.figure(figsize=(10,8))

        for j in range(len(envec)) :

            enstring2use = envec[j]

            chanvec = fromestring2chrange(enstring2use)

            CHANMIN, CHANMAX = chanvec[0], chanvec[1]
            figsize = [12,9]
            outfileStem = "{}_{}".format(outfile, enstring2use)

            #print ("CHANNEL MIN = {} \n CHANNEL MAX = {} \n GRADE SELECTION = {} \n SOURCE Coordinates (SKY X, SKY Y) = ({}, {})".format(CHANMIN, CHANMAX, grade, srcra, srcdec))

            enstringtmp = "{}-{}".format(round(CHANMIN/1e2,1), round(CHANMAX/1e2,1))

            #enstring = r"$\Delta$E : {} keV".format(enstringtmp) 

            eef, eeftbl = extract_cnts(evtfile = evtfilename, raddivec= outradivec, CHANMIN = CHANMIN, CHANMAX = CHANMAX, grade_flag = grade, sky_coord = sky_coord, coord = [srcra, srcdec], outfile = outfileStem)

            #eeftbl.write("test_before_bkgcrr_{}.tab".format(enstring2use), format='ascii.ipac', overwrite = True)
            if len(eeftbl) >= 1 :

                str4log += "Congratulations!!! The EEF table with counts for various radii is generated successfully : \n \t \t The Energy Range = {} \n \t \t Length of the output table = {} \n".format(enstringtmp, len(eeftbl))

            else :

                str4log += "Sad News !!! The EEF table with counts for various radii has failed \n"


            # Background Correction...of EEF profile
            if bkgcorrflag :
                str4log += "Background Correction Flag was 'yes', so performing it... \n"
                print ("Exposure : {}".format(float(evthdr['exposure'])))
                rad4corr = (eef[:,0] * 4.122)/60.0
                eef[:,1] = eef[:,1] - (float(evthdr['exposure']) * bkg_profile(enstring = enstring2use, x = rad4corr))
                eef[:,2] = np.sqrt(eef[:,1])
                eeftbl['Counts'] = eeftbl['Counts'] - (float(evthdr['exposure']) * bkg_profile(enstring = enstring2use, x = eeftbl['Radius2']))
                eeftbl['Error'] = np.sqrt(eeftbl['Counts'])

            #eeftbl.write("test_after_bkgcrr_{}.tab".format(enstring2use), format='ascii.ipac', overwrite = True)

            if j < 1 : 
                fig2, ax2 = plt.subplots(figsize = [12,9])

                normeef, tbl2ret, ax2 = eef_profilePlooter(eef, outfile = "{}.png".format(outfileStem), evtheader = evthdr, simeeffile = simeeffile, ipacfilename = ipacfilename, enstring = enstringtmp, plttype = 'lead4multi', color2plt = color2plt[j], getval  = getval, rad2ret = rad2ret, tbl2ret = None, NrefStyl = NrefStyl, plttypflg4romreg = j, ax2 = ax2, usepoly3 = usepoly3)

            else :

                normeef, tbl2ret, ax2 = eef_profilePlooter(eef, outfile = "{}.png".format(outfileStem), evtheader = evthdr, simeeffile = simeeffile, ipacfilename = ipacfilename, enstring = enstringtmp, plttype = None, color2plt = color2plt[j], getval = True, rad2ret = rad2ret, tbl2ret = tbl2ret, NrefStyl = NrefStyl, plttypflg4romreg = j, ax2 = ax2, usepoly3 = usepoly3)

        if len(tbl2ret) >= 1 : 
            str4log += "The EEF plotter has successfully generated important plots...!!\n"
        else :
            str4log += "The EEF plotter has failed to generate important plots...!!\n"

        outtblname = "{}_finlTbl.ipac".format(outfileStem)

        tbl2ret.write(outtblname, format="ascii.ipac", overwrite = True)

        ax2.legend(loc='best', fontsize=12)

        fig2.savefig(outpltfl)
        
        movefileswithstr(str2search =  "{}*ipac".format(os.path.splitext(outfile)[0]), SrcDir = './', OutDir = outputdir2prod)
        movefileswithstr(str2search =  "{}*png".format(os.path.splitext(outfile)[0]), SrcDir = './', OutDir = outputdir2prod)


        #plt.show()
  

    elif str(process).upper() in ["PLOT", "PL", "P"] :

        enstring = r"$\Delta$E : {}-{} keV".format(round(CHANMIN/1e2,1), round(CHANMAX/1e2,1))

        if ipacfilename != None :

            filetoplot = ipacfilename

        else :

            filetoplot = raw_input("Enter the name of ipac file to be plotted : ")

        inptbl = Table.read(filetoplot, format="ascii.ipac")

        eef = np.column_stack((inptbl['Radius2'], inptbl['Counts'], inptbl['Error']))

        normeef, tab2ret = eef_profilePlooter(eef, outfile = "{}.png".format(outfile), evtheader = evthdr, simeeffile = simeeffile, ipacfilename = ipacfilename, enstring = enstring)


        outtblname = "{}_finlTbl.ipac".format(outfile)

        tbl2ret.write(outtblname, format="ascii.ipac", overwrite = True)

        plt.legend(loc='best', fontsize=12)

        plt.savefig(outpltfl)

        movefileswithstr(str2search =  "{}*ipac".format(os.path.splitext(outfile)[0]), SrcDir = './', OutDir = outputdir2prod)
        movefileswithstr(str2search =  "{}*png".format(os.path.splitext(outfile)[0]), SrcDir = './', OutDir = outputdir2prod)
        #plt.show()

    
    else :

        tbl2ret = None


    print (str4log)

    logfl.write(str4log)

    logfl.close()

    return tbl2ret


def movefileswithstr(str2search =  None, SrcDir = None, OutDir = None) :

    import glob, os

    if SrcDir == None :
        SrcDir = os.getcwd()

    fileslist = glob.glob(os.path.join(SrcDir, str2search))

    if len(fileslist) > 0 :

        for file2cpy in fileslist :

            basefile2cpy = os.path.basename(file2cpy)

            os.rename(file2cpy, os.path.join(SrcDir, OutDir, basefile2cpy))


def fromestring2chrange(estring) :
    import os

    lowchan = float( estring.split("to")[0].split("p")[0]+"."+ estring.split("to")[0].split("p")[1] ) * 100 + 1
    highchan = float( estring.split("to")[1].split("p")[0]+"."+ estring.split("to")[1].split("p")[1] ) * 100 + 1
    chanvec = [lowchan, highchan]

    return chanvec
   



def extract_cnts(evtfile = "evtfile.evt", raddivec= None, CHANMIN=30, CHANMAX=801, grade_flag="0-12", sky_coord = 1, coord = [500, 500], outfile = "outfile", curvebinsize = 118.873645 ) :

    import numpy as np
    import matplotlib.pyplot as plt
    import glob, os
    from astropy.io import fits as pyfits
    from astropy.table import Table, Column
    from scipy.interpolate import interp1d

    srcra = coord[0]; srcdec = coord[1]
    #print ("Coordinates : {} & {}".format(srcra, srcdec))
    cwdir = os.getcwd()
    loopcnt=0
    eef = ([])
    eeftbl = Table()

    for i in raddivec :

        #print ("S. No. of loop : {}".format(i))
        radius = i
        regfilename = "regfile_{}.reg".format(loopcnt)
        if not os.path.exists(os.path.join(cwdir,regfilename)) :
            regf = open("{}".format(regfilename),'w')
            regf.write('# Region file format: DS9 version 4.1\n')
            regf.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
            regf.write("image\n")
            regf.write("circle({},{},{})\n".format(srcra, srcdec, radius))
            regf.close()

        outimgname = "tmpimg_{}".format(loopcnt)

        imgStts = FrameProductMaker(evtfile,regfile = regfilename, CHANMIN = CHANMIN, CHANMAX = CHANMAX, productstem = outimgname, grade_flag = grade_flag, curvebinsize = 120, sky_coord = sky_coord)

        #imgfname = glob.glob("tmpimg*.img")
    
        if os.path.exists(os.path.join("./","{}.img".format(outimgname))) :

            outimgname = "{}.img".format(outimgname)

            imgf = pyfits.open(outimgname)
    
            counts = np.sum(imgf[0].data)
            errorcnts = np.sqrt(counts)

            imgf.close()
            eef.append([i, counts, errorcnts])
            #dtout.append()
    
            os.remove("{}".format(outimgname))
            os.remove("{}".format(regfilename))
        loopcnt = loopcnt + 1

    eef = np.array(eef)
    
    rad = eef[:,0]*(4.122/60.)
    eeftbl.add_column(Column(eef[:,0], name = "Radius"), index=0)
    eeftbl.add_column(Column(eef[:,1], name = "Counts"), index=1)
    eeftbl.add_column(Column(eef[:,2], name = "Error"), index=2)
    eeftbl.add_column(Column(rad, name = "Radius2"), index=3)
    eeftbl.write("{}_out.ipac".format(outfile), format="ascii.ipac", overwrite=True)
    return eef, eeftbl 


def eef_profilePlooter(eef, outfile = "outfile.png", evtheader = 'evthdr', simeeffile = 'simeeffile', ipacfilename = 'ipacfilename', enstring = "0.3-7.0 keV", getval=True, rad2ret = [1.0, 2.0, 3.0], plttype = 'lead4multi', color2plt = 'r', tbl2ret = None, NrefStyl = 0, stepsize = 8, dtmarker = 'o', markersize = 8, fillstyle = 'none', plttypflg4romreg = 0, normradii = 20.0, nstd =2., ax2 = '', usepoly3 = 'no' ) :

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.table import Table, Column
    from scipy.interpolate import interp1d
    from scipy import optimize
    import os


    srcname = evtheader['object']
    exposure = float(evtheader['exposure'])
    obsid = evtheader['OBS_ID']
    
    enstringtmp = enstring.split("-")[0].split(".")[0] + "p" + enstring.split("-")[0].split(".")[1] + "to" + enstring.split("-")[1].split(".")[0] + "p" + enstring.split("-")[1].split(".")[1] 
    
    srcname = (srcname.replace(" ", "")).upper()

    exposure = round(exposure/1e3, 2)
    
    x = eef[:,0]; y = eef[:,1]; ey = eef[:,2]

    x = x*(4.122/60.); ex = stepsize * (4.122/60.)

    ex =np.repeat(ex, len(x))
    #-----------------------------------

    y4fit = y/np.mean(y)
    ey4fit = ey/np.mean(y)

    #     Model EEF using a suits of models....... Edit it in Morning.....!!!!!!!!!!!
    #------------------------------------------------------------------
    enstring4fl = ((enstring.replace(".", 'p')).replace('-', 'to')).replace(' ', '')

    testpltfname = "{}_{}_chk.png".format(os.path.splitext(outfile)[0], enstring4fl)

    rfit = np.arange(0, 22.0, 0.125)

    #tframe = ''
    
    bestfitredchisq, dof4bestfit, bestfitwerr, bestfitmdl, bestfitmdlstr = run_eef_fitsuits(x, y4fit, ex, ey4fit, pltindx = 0, testpltfname = testpltfname,  figsize = [12, 9], debug = True, usepoly3 = usepoly3)

    print (bestfitwerr[:,0])

    Norm4romfit = bestfitmdl(bestfitwerr[:,0], [normradii])[0]

    print ("-------------------------------{}\n".format(Norm4romfit))

    fit_up = bestfitmdl( bestfitwerr[:,0] + nstd * bestfitwerr[:,1], rfit)/Norm4romfit
    fit_dw = bestfitmdl( bestfitwerr[:,0] - nstd * bestfitwerr[:,1], rfit)/Norm4romfit

    fittab = np.column_stack((rfit, bestfitmdl(bestfitwerr[:,0], rfit)/Norm4romfit))
    #plot the bestfit mdl...

    if plttypflg4romreg < 1:

        #fig2, ax2 = plt.subplots(figsize = [12, 9])

        ax2.plot(rfit, bestfitmdl(bestfitwerr[:,0], rfit)/Norm4romfit, linestyle = '--', color = color2plt, lw = 2.5, label = bestfitmdlstr)

    else :

        ax2.plot(rfit, bestfitmdl(bestfitwerr[:,0], rfit)/Norm4romfit, linestyle = '--', color = color2plt, lw = 2.5, label = '')

    ax2.fill_between(rfit, fit_up, fit_dw, color = 'y', alpha=.25, label='')


    #For linear interpoolation of profile...
    profinterpol = interp1d(x, y4fit, kind='linear', fill_value = 'extrapolate')

    print ("The Mode of Profile Normalization : {}".format(NrefStyl))

    if NrefStyl == 1 :

        yref = profinterpol(20.0)
        print ("NormRefStyle : {}".format("Fixed 20 arcm from extrapolated data in case we do not cover that energy..."))
        y = y4fit/yref  
        ey = ey4fit/yref
        ax2.errorbar(x, y, yerr=ey, xerr=0, marker=dtmarker, color = color2plt, fillstyle = fillstyle, markersize=markersize, linestyle ='', label = "data : {}".format(enstring))

        func2use4levels = interp1d(y[::-1], x[::-1], kind='linear', fill_value = 'extrapolate')

    elif NrefStyl == 2 :

        yref1 = profinterpol(20.0)
        yref2 = y4fit[-1]
        ytemp = y4fit/yref1
        eytemp = ey4fit/yref1
        ax2.errorbar(x, ytemp, yerr=eytemp, xerr=0, marker=dtmarker, color = color2plt, fillstyle = fillstyle, markersize=markersize, linestyle ='', label = "data-ExtraPol : {}".format(enstring))

        ytemp = y4fit/yref2
        eytemp = ey4fit/yref2
        ax2.errorbar(x, ytemp, yerr=eytemp, xerr=0, marker=dtmarker, color = color2plt, fillstyle = fillstyle, markersize=markersize, linestyle ='', label = "data-LatRad : {}".format(enstring))

        ytemp = y4fit/Norm4romfit
        eytemp = ey4fit/Norm4romfit
        ax2.errorbar(x, ytemp, yerr=eytemp, xerr=0, marker=dtmarker, color = color2plt, fillstyle = fillstyle, markersize=markersize, linestyle ='', label = "data-BestFit : {}".format(enstring))

        yref = y4fit[-1]
        y = y4fit/yref  
        ey = ey4fit/yref

        func2use4levels = interp1d(y[::-1], x[::-1], kind='linear', fill_value = 'extrapolate')
  
    elif NrefStyl == 3 :

        print ("NormRefStyle : {}".format("Largest possible radius.."))

        yref = profinterpol(normradii)
        y = y4fit/yref  
        ey = ey4fit/yref
        ax2.errorbar(x, y, yerr=ey, xerr=0, marker=dtmarker, color = color2plt, fillstyle = fillstyle, markersize=markersize, linestyle ='', label = "data : {}".format(enstring))

        func2use4levels = interp1d(y[::-1], x[::-1], kind='linear', fill_value = 'extrapolate')

    else :

        print ("NormRefStyle : {}".format("20 arcm at Fitted Profile; Kings' Profile.."))

        yref = Norm4romfit
        y = y4fit/yref  
        ey = ey4fit/yref
        ax2.errorbar(x, y, yerr=ey, xerr=0, marker=dtmarker, color = color2plt, fillstyle = fillstyle, markersize=markersize, linestyle ='', label = "data-BestFit : {}".format(enstring))

        #Function to use for extract info about 60 % and 90 % containments
        
        #fitsintblform = bestfitproffunc(fitpar, np.arange(0, 23, 0.125))
        normfit = fittab[:,1]

        func2use4levels = interp1d(normfit[::-1], fittab[:,0][::-1], kind='linear', fill_value = 'extrapolate')

    normtbl = Table()
 
    # Interpolating the normlized curve

    if getval == True :

        f1 = interp1d(x, y, fill_value="extrapolate")
        xnew = np.linspace(0, 22, 100, endpoint = True)

        ax2.plot(xnew, f1(xnew),'-', color = color2plt, lw=2.5, label = '')

    if simeeffile != None :

        simeefdt = np.loadtxt(simeeffile, skiprows = 8, dtype="str")
        simpsf = np.array(simeefdt[:,1],float)
        simrad = np.array(simeefdt[:,0],float)*(4.13/60) 
        simeef = np.array(simeefdt[:,2],float)
        ax2.plot(simrad, simeef, color = 'k', linestyle = "-", lw = 2.5, label = "Simulated EEF")


    if plttype == 'lead4multi' :
        ylim = [min(y) - 0.25*(max(y) - min(y)), max(y) + 0.25*(max(y) - min(y)) ]

        if ylim[0] < 0:
            ylim[0]=0

        ax2.set_ylim(ylim)
        xlim = [0, 24]
        ax2.set_xlim(xlim)
        ax2.set_xlabel(r"Radius; r [$^{'}$]", fontsize=20)
        ax2.set_ylabel("Norm. Counts", fontsize=20)
        ax2.set_title("Encircled Energy Fraction; Source : {} & Exp. : {} ks".format(srcname, exposure), fontsize=17)
        ax2.text(15., 0.05, obsid, color='black', fontsize=15)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.tick_params(axis='both', which='minor', labelsize=13)

    
    # Interpolating the values to get at 90% position; not fitting a proper function
    xcomp90 = func2use4levels(0.90)
    xcomp60 = func2use4levels(0.60) 

    if plttype == 'lead4multi' :
        #print (xcomp90)
        ax2.axvline(x = xcomp90, linestyle = '--', color=color2plt, lw = 2.5, label = "90 % line")
        ax2.axvline(x = xcomp60, linestyle = '-.', color=color2plt, lw = 2.5, label = "60 % line")
        ax2.axhline(y = 1.0, linestyle = '-.', color=color2plt, lw = 2.0)

        #x4refplt = np.arange(0, 22, 0.25)
        #print (x4refplt)
        #print (profinterpol(x4refplt))
        #plt.plot(x4refplt, profinterpol(x4refplt)/profinterpol(20.), color = 'k', lw =5.5, linestyle = '--', label ="ExtrpolLine" )

    else :
        ax2.axvline(x = xcomp90, linestyle = '--', color=color2plt, lw = 2.5, label = "")
        ax2.axvline(x = xcomp60, linestyle = '-.', color=color2plt, lw = 2.5, label = "")
        

    ax2.text((int(xcomp90) + 1 + plttypflg4romreg), 0.87-(0.06)*plttypflg4romreg, str(round(xcomp90,2)), color=color2plt, fontsize=15)
    ax2.text((int(xcomp60) + 1 + plttypflg4romreg), 0.57-(0.06)*plttypflg4romreg, str(round(xcomp60,2)), color=color2plt, fontsize=15)
    ax2.minorticks_on()
    
    normtbl.add_column(Column(x, name = "X"), index=0)
    normtbl.add_column(Column(y, name = "Y"), index=1)
    normtbl.add_column(Column(ey, name = "EY"), index=2)
    tbloutname = os.path.splitext(outfile)[0] + "_tab.ipac"
    normtbl.write(tbloutname, format="ascii.ipac", overwrite = True)
    
    #rad2ret = [1.0, 2.0, 3.0] 
    if getval == True :
        if tbl2ret == None :
            tbl2ret = Table()

        tbl2ret.add_column(Column(rad2ret, name = "R_{}".format(enstringtmp)), index = 0)
        tbl2ret.add_column(Column(f1(rad2ret), name = "lin_{}".format(enstringtmp)), index = 1)
        try :
            tbl2ret.add_column(Column(f2(rad2ret), name = "cub_{}".format(enstringtmp)), index = 2)
        except :
            print ("No output for cubic spline fit...")
    else :
        tbl2ret = None
    
    return eef, tbl2ret, ax2


def bkg_profile(enstring = 'e0p3to7p0', x = 'x') :

    import numpy as np

    fitpardict = {}
    fitpardict['0p3to7p0'] = {'a' : 3.76517e-4, 'b' : 6.92770e-4, 'c' : -7.50153e-4}
    fitpardict['0p3to0p7'] = {'a' : 9.328271e-05, 'b' : 9.890137e-05, 'c' : -1.227563e-4}
    fitpardict['0p7to2p0'] = {'a' : 1.7934009e-4, 'b' : 4.361011e-04, 'c' : -4.640430e-4}
    fitpardict['2p0to7p0'] = {'a' : 1.0650766e-4, 'b' : 1.645912e-04, 'c' : -1.675936e-4}
    fitpardict['0p3to5p0'] = {'a' : 3.4914665e-4, 'b' : 7.061383e-04, 'c' : -7.915850e-4}
    fitpardict['0p5to2p0'] = {'a' : 2.1830701e-4, 'b' : 5.0626015e-04, 'c' : -5.338131e-4}
    fitpardict['0p5to5p0'] = {'a' : 2.9667819e-4, 'b' : 6.85923397e-04, 'c' : -7.459751e-4}
    fitpardict['0p5to7p0'] = {'a' : 3.24109126e-4, 'b' : 6.71443674e-04, 'c' : -7.01485e-4}
    fitpardict['0p7to5p0'] = {'a' : 2.57705021e-4, 'b' : 6.15889145e-04, 'c' : -6.76593e-4}
    fitpardict['0p7to7p0'] = {'a' : 2.84025679e-4, 'b' : 6.16493538e-04, 'c' : -6163716e-4}
    fitpardict['1p0to7p0'] = {'a' : 2.29400128e-4, 'b' : 4.304600974e-04, 'c' : -4368962e-4}
    fitpardict['3p0to7p0'] = {'a' : 6.69394462e-05, 'b' : 8.97810539e-05, 'c' : -7.348099e-05}
    fitpardict['5p0to7p0'] = {'a' : 2.6519467e-05, 'b' : 3.7466461e-06, 'c' : -4.48550894e-06}

    def func(p, x):

        return p['a'] * x*x + p['b']*x + p['c']

    try :
        p = fitpardict[enstring]
        #print (p)
        return func(p, x)
        stts = 0

    except :
        stts = 999
        return np.repeat(0.0, len(x))


def get_centralcoord_lookup(evtfile, rawx=300.0, rawy = 300.0, coordtol=2e0, stepcoordtol=5e-1) :

    import numpy as np

    from astropy.io import fits
    from astropy.table import Table, Column, vstack

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
    print ( fltrdcoordtbl )

    print ( (np.mean(fltrdhdu['RAWX']), np.mean(fltrdhdu['RAWY']) ), (np.mean(fltrdhdu['X']), np.mean(fltrdhdu['Y'])) )

    return [np.mean(fltrdhdu['X']), np.mean(fltrdhdu['Y'])]



def run_eef_fitsuits(x, y, errx, erry, pltindx = 0, testpltfname = 'ModelFit_DebugFile.png', figsize = [12, 9], debug = True, usepoly3 = 'no') :

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.table import Table, Column
    import os


    # Important Definition Modules...
    def runodr(model, data, beta0 = [-0.00736605,  0.22761323, -0.18054454], modelstr = "Quadratic Model Fitting using ODR") :

        from scipy.odr import Data, Model, ODR, RealData, odr_stop

        x, y, errx, erry = data

        model2use = Model(model)

        data2use = RealData(x, y, sx = errx, sy = erry)

        odr2use = ODR(data2use, model2use, beta0 = beta0, maxit = 5000)

        out = odr2use.run()

        txtreport = "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        txtreport += "\n====== {} =======\n".format(modelstr)
        txtreport += "Fitted Parameters : \t {} \n".format(out.beta)
        txtreport += "Covariance errors: \t {} \n".format(out.cov_beta.diagonal())  
        txtreport += "Standard errors : \t {} \n".format(out.sd_beta)
        txtreport += "Minimum chi^2 : \t {} \n".format(out.sum_square)
        txtreport += "Minimum chi^2 [ reduced ] : \t {} \n".format(out.res_var)
        txtreport += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"

        print (txtreport)

        return [out.beta, out.sd_beta] , out.sum_square, (len(y) - len(out.beta))


    def runkmfitter(model, data,  beta0, modelstr = "Quadratic Model Fitting using Fitter") :

        import subprocess
        import numpy as np
        #----------------------------------------------------------------------
        try :
            from kapteyn import kmpfit

        except :

            print ("The module \'kapteyn\' not found .. Installing!!")
            cmd = "sudo pip install kapteyn"

            retcode = subprocess.call(cmd, shell=True)

            if retcode < 0:
                print("Child was terminated by signal", -retcode,  sys.stderr)

            else:
                print("Child returned", retcode,  sys.stderr)

            import kapteyn

        #----------------------------------------------------------------------
        # Residulal used for fitting data

        def residuals(p, data):

            # Fitting Quadratic for data with Errors in Y only
            #a, b, c = p

            x, y, ex, ey = data
            w = ey*ey
            wi = np.sqrt(np.where(w==0.0, 0.0, 1.0/(w)))
            d = wi*(y-model(p,x))

            return d


        # Prepare fit routine
        fitobj = kmpfit.Fitter(residuals = residuals, data = (x, y, errx, erry), xtol=1e-12, gtol=1e-12)
    
        fitobj.fit(params0 = beta0)

        txtreport = "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        txtreport += "\n====== {} =======\n".format(modelstr)
        txtreport += "Fitted parameters: \t {} \n ".format(fitobj.params)
        txtreport += "Covariance errors: \t {} \n ".format(fitobj.xerror)
        txtreport += "Standard errors: \t {} \n ".format(fitobj.stderr)
        txtreport += "Chi^2 min :  \t {} \n ".format(fitobj.chi2_min)
        txtreport += "Reduced Chi^2 min \t {} \n ".format(fitobj.rchi2_min)
        txtreport += "Status Message \t : {} \n".format(fitobj.message)
        txtreport += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"

        print (txtreport)

        return [fitobj.params, fitobj.stderr], fitobj.chi2_min, (len(y) - len(fitobj.params))


    def quadmodel(p, x):

       # Model: Y = a + b*x + c*x*x
       a,b,c = p
       x = np.array(x)
       return a + b*x + c*x*x


    def poly3model(p, x):
       # Model: Y = a + b*x + c*x*x
       a,b,c,d = p
       x = np.array(x)
       return d + c*x + b*x*x + a*x*x*x 


    def GetIntegKingsProf(p, x) :

        from scipy.integrate import quad, dblquad

        A, x0, r = p

        def KingsProf(x, A, x0, r) :

            return A * (1 + (x/x0)**2)**-r

        xlim0 = 0

        val2ret = ([])

        for xlim in x :

            val2ret.append([xlim, quad(KingsProf, xlim0, xlim, args = (A, x0, r) )[0]] )

        val2ret = np.array(val2ret)
    
        return val2ret[:,1]

    #Running ODR  using scipy.odr & KMFITTER from Kapteyn

    #Defining the data list
    data = [x, y, errx, erry]

    # Tested Initial Values...and Models to fit..
    if usepoly3.upper() in ["YES", "Y"] : 
        invallist =  [ [ [-0.00736605,  0.22761323, -0.18054454], [1.72445564e-04, -1.04219887e-02,  2.18305129e-01, -1.30875794e-01], [0.35, 8., 3.] ], [quadmodel, poly3model, GetIntegKingsProf], ["Quad. ", "Poly. O3", "Integ. Kings"] ]
    else :
        invallist =  [ [ [-0.00736605,  0.22761323, -0.18054454],  [0.35, 8., 3.] ], [quadmodel, GetIntegKingsProf], ["Quad. ", "Integ. Kings"] ]

    #Model to use for tests
    models2use = invallist[1]

    #Initial Values to start with
    invals2use = invallist[0]

    print ("Length of invallist : {} \n Models to use : {} \n Initial Values to use : {} \n".format( len(invallist), models2use, invals2use)   )

    #Main def files which calls two mode of fitting algorithms
    modelinstance = [[runodr, runkmfitter],["ODR", "KMFIT"]]

    #initalize dictionary to store fit output [free-style dict.]
    fitoutdict = {}

    
    #Define (sub)plot object for storing the Fit Quality Cross-Check figure... 
    if testpltfname == None :
        testpltfname = "Cross_Check_plot_4random_energyband.png"
 
    fig1, tframe = plt.subplots(figsize = figsize)

    #Setting up Axis Labels for Cross-check plot
    tframe.set_xlabel("Radius [arcm]", fontsize =17)
    tframe.set_ylabel("Norm. EEF Prof.", fontsize =17)
    tframe.set_title(r"ODR/kmpfit with weighted fit. Model: $y=a+bx+cx^2$ and  A [ 1 + ($\frac{x}{x0})^2]^{-r}$ ", fontsize = 17)

    tframe.errorbar(x, y, xerr=errx, yerr=erry,  marker='o', color ='k', ecolor = 'gray', markersize = 9, linestyle = '', label = '-data', fillstyle = 'none')
    

    linestylevec = ['-', '--', ':', '-.']

    colorvec = ['r', 'b', 'g', 'c', 'm', 'y']

    prodoutlist = ([])

    for i in range(len(modelinstance)) :

        mdlinstance = modelinstance[0][i]

        mainfitcallerstr = modelinstance[1][i]

        #print (modelinstance[0][i].func_name)

        for j in range(len(invallist[0])) :

            inval = invallist[0][j]
            mdl2fit = invallist[1][j]
            #print (mdl2fit)
            modstr4plt = invallist[2][j]

            str2prntOnplt = "{}-{}".format(modstr4plt, mainfitcallerstr)

            try :

                params, minchisq, dof = mdlinstance(mdl2fit, data, beta0 = inval, modelstr = str2prntOnplt)

            except :

                continue 

            fitoutdict['PARAMS{}{}'.format(i, j)] = np.array(np.column_stack((params[0], params[1])), float)

            fitoutdict['MINCHISQR{}{}'.format(i, j)] = minchisq

            prodoutlist.append([minchisq, params, dof, mdl2fit, str2prntOnplt ])

            #Adding specific plot to plot objects...main and cross-check, both...
            tframe.plot(x, mdl2fit(params[0], x),  linestyle = linestylevec[j], color = colorvec[i], lw = 2.5, label = str2prntOnplt)

    prodoutlist = np.array(prodoutlist)

    #print prodoutlist

    #index of minimum chi^2 value for best fit model..
    index2use = np.where(prodoutlist[:,0]==min(prodoutlist[:,0]))

    #print (np.array(prodoutlist[index2use, 1][0][0],float)[1])
    text2debug = "____________________________________________________________________________________________\n"
    text2debug += "\tMin. \chi^2 : {} \n ".format( prodoutlist[index2use, 0][0][0]) 

    for par, perr in zip(prodoutlist[index2use, 1][0][0][0], prodoutlist[index2use, 1][0][0][1]) :

        text2debug += "Best Fit Parameters with Errors : {} +- {} \n ".format(par, perr)

    text2debug += "Degree of Freedom for Best Fit Model : {} \n ".format( prodoutlist[index2use, 2][0][0]  )
    
    text2debug += "The Name of the Best Fit Model : {} \n ".format( (prodoutlist[index2use, 3][0][0]).func_name  )

    text2debug += "The String Specifier for Best Fit Model : {}\n ".format( prodoutlist[index2use, 4][0][0]  )

    text2debug += "____________________________________________________________________________________________\n"

    if debug :
        print (text2debug)
    #print(np.array(prodoutlist)[:,1][2][0][0])

    #Adding plot to main plotobj...
    bestfitredchisq = prodoutlist[index2use, 0][0][0]
    bestfitparams = prodoutlist[index2use, 1][0][0][0]
    bestfitparerr = prodoutlist[index2use, 1][0][0][1]
    dof4bestfit = prodoutlist[index2use, 2][0][0]
    bestfitmdl = prodoutlist[index2use, 3][0][0]
    bestfitmdlstr = prodoutlist[index2use, 4][0][0]

    bestfitwerr = np.array(np.column_stack((bestfitparams, bestfitparerr)), float)
    
    #Final tweaks to cross-check plot
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tick_params(axis='both', which='minor', labelsize=13)
    plt.minorticks_on()
    tframe.grid(True)
    tframe.legend(loc=2)
    fig1.savefig(testpltfname)
    #plt.close()
    
    return bestfitredchisq, dof4bestfit, bestfitwerr, bestfitmdl, bestfitmdlstr


def GetIntegKingsProf2run(A, x0, r, x) :

    from scipy.integrate import quad, dblquad
    import numpy as np

    #A, x0, r = p

    def KingsProf(x, A, x0, r) :

        return A * (1 + (x/x0)**2)**-r

    xlim0 = 0

    val2ret = ([])

    for xlim in x :

        val2ret.append([xlim, quad(KingsProf, xlim0, xlim, args = (A, x0, r) )[0]] )

    val2ret = np.array(val2ret)
    
    return val2ret[:,1]


if __name__ == "__main__":
    import sys
    sys.settrace

    scriptname = sys.argv[0]

    print_preamble()

    sys.setrecursionlimit(500000)

    main()

 




