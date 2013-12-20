#authors Vlasios Vasileiou with help from Giacomo Vianello, Nicola Omodei

import ROOT,os,sys

try: 
	print("Loading BKGE...");
	ROOT.gSystem.Load("libBKGE.so")
except:
	print("Could not load BKGE library (file libBKGE.so). Did you compile it?")
	sys.exit()
print "Success!"	

from ROOT import TOOLS, BKGE_NS
ROOT.gStyle.SetOptStat(0)
TOOLS.Set("GRB_NAME","") #this is to create the datadir directory of the bkg estimator (kluge)
if(os.environ.get("BKGE_DATADIR")==None):
  TOOLS.Set("BKGE_DATADIR",os.path.join(os.path.dirname(__file__),"BKGE_Datafiles")+os.path.sep)
else:
  TOOLS.Set("BKGE_DATADIR",(os.environ.get("BKGE_DATADIR")+os.path.sep).replace("%s%s" % (os.path.sep,os.path.sep),os.path.sep))
pass
ResponseFunction="P7REP_TRANSIENT_V15"     #only P7TRANSIENT reprocessed class is currently publicly supported -- this is the same P7REP_TRANSIENT
TOOLS.Set("BKG_ESTIMATE_ERROR",0.15)  #see associated publication for where this number can from

def CalculateBackground(start, stop , grb_trigger_time, RA, DEC, FT1, FT2, OUTPUT_DIR="output/", emin=-1, emax=-1, ebins=-1, chatter=1, overwrite=False, EvaluateMaps=True, CalcResiduals=True,
    ROI_Calculate=1, ROI_Containment=0.95, ROI_Localization_Error=0, ROI_Radius=12, ROI_Max_Radius=12, GRB_NAME="",ROI_RadiusFile="" ):
    '''
    Calculate the Background for a single observation.
    '''
    
    if not os.path.isfile(FT1):
       print "FT1 %s file is not found" %FT1
       sys.exit()

    if not os.path.isfile(FT2):
       print "FT2 %s file is not found" %FT2
       sys.exit()
       

    #These will be used for naming folders and files  
    if (GRB_NAME==""):
       GRB_NAME="GRB_MET_%.2f" %grb_trigger_time
    TOOLS.Set("GRB_NAME",GRB_NAME)
    TOOLS.Set("OUTPUT_DIR",OUTPUT_DIR)

    
    #Setup the ROI radius
    
    TOOLS.Set("CALCULATE_ROI",ROI_Calculate)
    if (ROI_Calculate==1):
        TOOLS.Set("ROI_MAX_RADIUS",ROI_Max_Radius)
        TOOLS.Set("ROI_CONTAINMENT",ROI_Containment)
    elif (ROI_Calculate==0) :
        TOOLS.Set("ROI_MAX_RADIUS",ROI_Radius) #the C++ code picks up the desired radius from the MAX_RADIUS variable if ROI_Calculate==2
    elif (ROI_Calculate==2):
        if (ROI_RadiusFile=="") :
           print("Please provide a ROI_RadiusFile")
           sys.exit()
        else:
           TOOLS.Set("ROI_RADIUSFILE",ROI_RadiusFile)
    TOOLS.Set("ROI_LOCALIZATION_ERROR",ROI_Localization_Error)
    ###############################
  
    TOOLS.Set("GRB_RA",RA); TOOLS.Set("GRB_DEC",DEC)
        
    duration = stop - start
    if chatter>1: print ' ====> CalculateBackground(start=%.2f,stop=%.2f) (duration =%.2f)' %(start,stop,duration)
    myname=sys._getframe().f_code.co_name

    if chatter>2: TOOLS.PrintConfig()

    BKGE_NS.CalculateBackground("%.2f_%.2f/" %(start,stop), grb_trigger_time+start, duration, FT1,FT2,ResponseFunction,emin,emax,ebins,chatter,CalcResiduals)
    ResultsFilename=""
    if EvaluateMaps:
        if chatter: print "%s: Counting actually detected events and plotting estimated background...\n" %myname
        #if you ask for bkg on an interval that starts before the trigger time, calculate the ROI for the pointing exactly at trigger time
        if start<0: 
              MET_FOR_THETA=grb_trigger_time
              if chatter>1: print "Calculating ROI Radius for an off-axis angle corresponding to the beginning of the observation."
        else :
              MET_FOR_THETA=grb_trigger_time+start
              if chatter>1: print "Calculating ROI Radius for an off-axis angle corresponding to the middle of the observation."
             
        ####
        ResultsFilename = BKGE_NS.PlotBackground("%.2f_%.2f" %(start,stop), grb_trigger_time+start, duration, FT1,FT2,ResponseFunction,emin,emax,ebins,overwrite,chatter,MET_FOR_THETA)
    
        ROOTFile=ROOT.TFile(ResultsFilename,"open")
        BKGE_NDET = float(ROOTFile.BKGE_NDET.GetTitle())
        BKGE_NEXP = float(ROOTFile.BKGE_NEXP.GetTitle())
        BKGE_SIGNIF = float(ROOTFile.BKGE_SIGNIF.GetTitle())
        BKGE_SIGNIF_WITH_UNCERTAINTY = float(ROOTFile.BKGE_SIGNIF_WITH_UNCERTAINTY.GetTitle())
        return BKGE_NDET,BKGE_NEXP,BKGE_SIGNIF,BKGE_SIGNIF_WITH_UNCERTAINTY
    else :
        return -1,-1,-1,-1



def MakeGtLikeTemplate(start, stop, grb_trigger_time, RA, DEC, FT1, FT2, OUTPUT_DIR="output/", chatter=1, ROI_Radius=15, GRB_NAME="" ):
    if (GRB_NAME==""):
       GRB_NAME="GRB_MET_%.2f" %grb_trigger_time
       
    if chatter: print("Using a constant ROI radius of %d deg\n" %ROI_Radius)
    CalculateBackground(start,stop,grb_trigger_time, RA,DEC, FT1, FT2, OUTPUT_DIR, EvaluateMaps=False, ROI_Calculate=0, ROI_Radius=ROI_Radius,  GRB_NAME=GRB_NAME, chatter=chatter )
    BKGE_NS.MakeGtLikeTemplate(ROI_Radius, OUTPUT_DIR+'/'+'Bkg_Estimates'+'/%.2f_%.2f/' %(start,stop), ResponseFunction)





def Make_BKG_PHA(start, stop, grb_trigger_time, RA, DEC, FT1, FT2, emin, emax, ebins, ROI_Calculate, OUTPUT_DIR="output/", chatter=1, overwrite=False, ROI_Containment=0.95, 
    ROI_Localization_Error=0, ROI_Max_Radius=12, ROI_Radius=12, GRB_NAME="" , ROI_RadiusFile=""):
    ''' Produces a PHA I file containing the expected background of a single observation'''
    print "Warning: Code is in BETA version!"    
    import pyfits,commands
    
    if   ROI_Calculate==0:   suffix='flat_ROI'
    elif ROI_Calculate==1:   suffix='energy_dependent_ROI'
    elif ROI_Calculate==2:   suffix='user_defined_ROI'
    else :
       print("Invalid value for ROI_Calculate %d" %ROI_Calculate)
       sys.exit()

    duration=stop-start
    if (GRB_NAME==""):
       GRB_NAME="GRB_MET_%.2f" %grb_trigger_time

    pha_filename = "%s_%.2f_%.2f_%s.pha" %(ResponseFunction,start,duration,suffix)
    output_path  = "%s/%s/%.2f_%.2f" %(OUTPUT_DIR,'Bkg_Estimates',start,stop)
    if chatter: print "Results will be saved in file %s/%s" %(output_path,pha_filename)
        
    CalculateBackground(start, stop , grb_trigger_time, RA, DEC, FT1, FT2, OUTPUT_DIR, emin, emax, ebins, chatter, overwrite,True, True,\
       ROI_Calculate, ROI_Containment, ROI_Localization_Error, ROI_Radius, ROI_Max_Radius, GRB_NAME, ROI_RadiusFile)
     
    dfile = ROOT.TFile("%s/%s_bkg_%.0f_%.0f.root" %(output_path,ResponseFunction,emin,emax),"r")
    bkg = dfile.hCtsvsEnergy_Est
    roi = dfile.hROI
    tmp_filename = "%s/%s_bkg_for_PHA_%.0f_%.0f.txt" %(output_path,ResponseFunction,emin,emax)
    tmp_file = open(tmp_filename,"w")
    try: 
       for i in range(1,bkg.GetNbinsX()+1):
           tmp_file.write("%d %e\n" %(i,bkg.GetBinContent(i)/duration))
       tmp_file.close()
    except:
       print "BKGE did not finish ok!"
       sys.exit()
       
    #I am doing this cd to dir thing and then execute ascii2pha locally because (if I remember correctly) ascii2pha had some maximum arguments length that was being broken by a too long command.   
    tmp_filename = "%s_bkg_for_PHA_%.0f_%.0f.txt" %(ResponseFunction,emin,emax)
    cmd = "cd %s/%s/%.2f_%.2f ;ascii2pha infile=%s chantype=PI outfile=%s chanpres=yes dtype=2 qerror=no rows=- fchan=1 detchans=%d pois=no telescope=GLAST instrume=LAT detnam=LAT filter=NONE tlmin=1 exposure=%f clobber=yes" \
        %(OUTPUT_DIR,'Bkg_Estimates',start,stop,tmp_filename,pha_filename,ebins,duration)
    status,output=commands.getstatusoutput(cmd)
    if chatter:
        print cmd
        print output
    pass
    
    pha_file = pyfits.open("%s/%s" %(output_path,pha_filename),"update")
    pha_file[1].header.add_comment("Created by the BKGE - vlasisva@gmail.com")
    del pha_file[1].header['SYS_ERR']
    del pha_file[1].header['STAT_ERR']

    #add statistical error (1%)
    stat_err=[];sys_err=[]
    for ibkg in range(0,len(pha_file[1].data)):
        #print pha_file[1].data[ibkg][1]
        stat_err.append(0.01*pha_file[1].data[ibkg][1]) #STAT_ERR is in units of events
        sys_err.append(TOOLS.Get("BKG_ESTIMATE_ERROR")) #SYST_ERR is in units of percent (fractional)
    #stat_err_col=pyfits.Column(name="STAT_ERR",format="E",array=numpy.ones(len(pha_file[1].data)))
    stat_err_col=pyfits.Column(name="STAT_ERR",format="E",array=stat_err)
    sys_err_col=pyfits.Column(name="SYS_ERR",format="E",array=sys_err)
    new_table = pyfits.new_table(pha_file[1].columns+stat_err_col+sys_err_col)
    for i in pha_file[1].header.items():
        if not (i[0] in new_table.header):
            new_table.header.update(i[0],i[1])
    pha_file[1]=new_table

    #add energies
    mine=[];maxe=[]
    for i in range(1,bkg.GetNbinsX()+1):
        mine.append(pow(10,3+bkg.GetXaxis().GetBinLowEdge(i)))
        maxe.append(pow(10,3+bkg.GetXaxis().GetBinUpEdge(i)))
    chan_col=pyfits.Column(name="CHANNEL",format="I",array=range(1,bkg.GetNbinsX()+1))
    mine_col=pyfits.Column(name="E_MIN",format="E",array=mine)
    maxe_col=pyfits.Column(name="E_MAX",format="E",array=maxe)
    ebounds_tab=pyfits.new_table(pyfits.ColDefs([chan_col,mine_col,maxe_col]))
    ebounds_tab.header.update('EXTNAME','EBOUNDS')
    ebounds_tab.header.update('CHANTYPE','PI')
    ebounds_tab.header.update('DETCHANS',ebins)
    ebounds_tab.header.update('TLMIN1',1)
    ebounds_tab.header.update('TLMAX1',ebins)
    pha_file.append(ebounds_tab)

    #add GTI
    start_col=pyfits.Column(name="START",format="D",array=[grb_trigger_time+start])
    stop_col=pyfits.Column(name="STOP",format="D",array=[grb_trigger_time+stop])
    gti_tab=pyfits.new_table(pyfits.ColDefs([start_col,stop_col]))
    gti_tab.header.update('EXTNAME','GTI')
    gti_tab.header.update('TIMESYS','TT')
    gti_tab.header.update('TIMEREF','LOCAL')
    gti_tab.header.update('TIMEZERO','0')
    gti_tab.header.update('TIMEUNIT','s')
    pha_file.append(gti_tab)


    #add ROI Radius
    roi_radius=[];
    for i in range(1,roi.GetNbinsX()+1):
        roi_radius.append(roi.GetBinContent(i))
    roi_col=pyfits.Column(name="ROI_RADIUS",format="E",array=roi_radius)
    roi_tab=pyfits.new_table(pyfits.ColDefs([chan_col,roi_col]))
    roi_tab.header.update('EXTNAME','ROI_RADIUS')
    roi_tab.header.update('DETCHANS',ebins)
    roi_tab.header.update('TLMIN1',1)
    roi_tab.header.update('TLMAX1',ebins)
    pha_file.append(roi_tab)


    pha_file.close()
    return ("%s/%s" %(output_path,pha_filename))



def Make_BKG_PHA2(grb_trigger_time, RA, DEC, FT1, FT2, emin, emax, ebins,
    ROI_Calculate, ROI_Localization_Error=0, ROI_Max_Radius=12, ROI_Radius=12, ROI_Containment=0.95, ROI_RadiusFile="", 
    OUTPUT_DIR="output/",GRB_NAME="" , 
    chatter=1, overwrite=False,     
    Time_bins_def_file="", start=-1, stop=-1, dt=-1  #these are used to define the time intervals
    ):
    ''' Produces a PHA II file containing the expected background of a series of observations'''

    print "Warning: Code is in BETA version!"
    if   ROI_Calculate==0:   suffix='flat_ROI'
    elif ROI_Calculate==1:   suffix='energy_dependent_ROI'
    elif ROI_Calculate==2:   suffix='user_defined_ROI'
    else :
       print("Invalid value for ROI_Calculate %d" %ROI_Calculate)
       sys.exit()


    import commands,pyfits,math
    import numpy as num

    
    if (GRB_NAME==""):
      GRB_NAME="GRB_MET_%.2f" %grb_trigger_time

    output_path  = "%s/%s/%.2f_%.2f" %(OUTPUT_DIR,'Bkg_Estimates',start,stop)
    
    #time bins setup
    at0=[]; at1=[]; adt=[]
    if Time_bins_def_file!="":
	if ".fit" in Time_bins_def_file:
	    dfile=pyfits.open(Time_bins_def_file)
	    at0=dfile[1].data.field('START')
	    at1=dfile[1].data.field('STOP')
	else :
	    f=open(Time_bins_def_file,"r")
	    while (1):
		string=f.readline()
		if string=="": break
		at0.append(float(string.split(" ")[0]))
		at1.append(float(string.split(" ")[1]))
	    pass
	pass
	Time_bins_def_file=Time_bins_def_file.split("/")[-1]
	pha_filename = "%s_bkg_%s" %(ResponseFunction,Time_bins_def_file)
	print "Using time definition file %s and ignoring any passed start,stop,dt" %Time_bins_def_file
    else: 
	if start==stop:
    	    print "start==stop ?\n"
    	    sys.exit()
	if dt==0:
	    print "dt=0?"
	    sys.exit()
	pass
	
	nbins = int(math.floor((stop-start+1e-1)/dt))
	print "Using %d time bins between %.2f and %.2f s post-trigger" %(nbins,start,stop)
	att0=start
	dt=float(dt) #keep this line, otherwise at0 and at1 might end up being integers and saving them in the FITS file will not work
	for i in range(0,nbins):
	    at0.append(att0)
	    at1.append(att0+dt)
	    att0+=dt
	pass
	if(at1[-1]-stop >= 1E-2):
	  at0.append(at1[-1])
	  at1.append(stop)
	pass
	pha_filename = "%s_bkg_%.2f_%.2f_%.2f" %(ResponseFunction,start,stop,dt)
    #########################

    pha_filename="%s_%s.pha" %(pha_filename,suffix)
    if chatter: print "Results will be saved in file %s/%s" %(output_path,pha_filename)

    adt=num.array(at1)-num.array(at0)
    if chatter: print "Using %d bins between %f-%f MeV" %(ebins,float(emin),float(emax))

    bkg_data=[]
    roi_data=[]
    stat_err_data=[]
    sys_err_data=[]; 
    for i in range(0,len(at0)): #1e-5 is for rounding errors
        CalculateBackground(at0[i],at1[i],grb_trigger_time, RA, DEC, FT1, FT2, OUTPUT_DIR, emin, emax, ebins, chatter, overwrite, True, True, ROI_Calculate, \
        ROI_Containment, ROI_Localization_Error, ROI_Radius, ROI_Max_Radius, GRB_NAME,ROI_RadiusFile)

        an_output_path  = "%s/%s/%.2f_%.2f" %(OUTPUT_DIR,'Bkg_Estimates',at0[i],at1[i])
        
        results_filename = "%s/%s_bkg_%.0f_%.0f.root" %(an_output_path,ResponseFunction,emin,emax)
        dfile = ROOT.TFile(results_filename,"r")
	bkg = dfile.hCtsvsEnergy_Est
	roi = dfile.hROI
	abkg_data=[]
	aroi_data=[]
	for ibin in range(0,ebins):
	    abkg_data.append(bkg.GetBinContent(ibin+1)/adt[i]) #convert from counts in dt to rate
	    aroi_data.append(roi.GetBinContent(ibin+1)) #convert from counts in dt to rate
	pass
	bkg_data.append(abkg_data)
	roi_data.append(aroi_data)
	stat_err_data.append(num.array(abkg_data)*0.01)
	sys_err_data.append(TOOLS.Get("BKG_ESTIMATE_ERROR"))
    pass

    try:
	os.remove(pha_filename)
    except: 
	pass
    hdu=pyfits.PrimaryHDU()
    hdu.header.update('FILETYPE','PHAII')
    hdu.header.update('FILE-VER','1.0.0')
    
    #add energies
    mine=[];maxe=[]
    for i in range(1,bkg.GetNbinsX()+1):
        mine.append(pow(10,3+bkg.GetXaxis().GetBinLowEdge(i)))
        maxe.append(pow(10,3+bkg.GetXaxis().GetBinUpEdge(i)))
    chan_col=pyfits.Column(name="CHANNEL",format="I",array=range(1,bkg.GetNbinsX()+1))
    mine_col=pyfits.Column(name="E_MIN",unit="keV",format="1E",array=mine)
    maxe_col=pyfits.Column(name="E_MAX",unit="keV",format="1E",array=maxe)
    ebounds_tab=pyfits.new_table(pyfits.ColDefs([chan_col,mine_col,maxe_col]))
    ebounds_tab.header.update('TLMIN1',1)
    ebounds_tab.header.update('TLMAX1',ebins)
    ebounds_tab.header.update('EXTNAME','EBOUNDS')
    
    #add rates  
    
    counts_col  =pyfits.Column(name="RATE",format="%dD" %ebins,unit="count/s",array=bkg_data)
    roi_col     =pyfits.Column(name="ROI_RADIUS",format="%dD" %ebins,unit="deg",array=roi_data)
    exposure_col=pyfits.Column(name="EXPOSURE",unit="s",format="1E",array=adt)
    quality_col =pyfits.Column(name="QUALITY",format="1I",array=[0,0])   
    time_col    =pyfits.Column(name="TIME",unit="s",format="1D",array=at0)
    etime_col   =pyfits.Column(name="ENDTIME",unit="s",format="1D",array=at1)
    stat_err_col=pyfits.Column(name="STAT_ERR",unit="",format="%dD" %ebins,array=stat_err_data)
    sys_err_col =pyfits.Column(name="SYS_ERR",unit="",format="1D",array=sys_err_data)
    
    data_tab=pyfits.new_table(pyfits.ColDefs([counts_col,exposure_col,quality_col,time_col,etime_col,stat_err_col,sys_err_col,roi_col]))
    data_tab.header.update('EXTNAME','SPECTRUM')
    data_tab.header.update('TELESCOP','GLAST')
    data_tab.header.update('INSTRUME','LAT')
    #data_tab.header.update('EXPOSURE','')
    data_tab.header.update('BACKFILE','none')
    data_tab.header.update('CORRFILE','none')
    data_tab.header.update('CORRSCL',1)
    data_tab.header.update('RESPFILE','none')
    data_tab.header.update('ANCRFILE','none')
    data_tab.header.update('POISERR',False)
    #data_tab.header.update('SYS_ERR',9990.15)
    data_tab.header.update('GROUPING',0)
    '''
    The CHANTYPE keyword specifies whether the channels referred to in the file 
    are as assigned by the detector electronics at the time of data collection
     (in which case CHANTYPE = PHA), or whether any corrections have been 
     applied. An example of the latter case is when the data has been re-mapped 
     onto a standard "pulse-Invariant" (PI) channel grid, in which case 
     CHANTYPE = PI (see also CAL/GEN/92-002, George et al. 1992, Section 7). 
    '''
    data_tab.header.update('CHANTYPE',"PI")
    data_tab.header.update('TIMESYS','TT')
    data_tab.header.update('TIMEREF','LOCAL')
    data_tab.header.update('TIMEZERO','0')
    data_tab.header.update('TIMEUNIT','s')
    data_tab.header.update('TSTART',"%f" %(start+grb_trigger_time))
    data_tab.header.update('TSTOP',"%f" %(stop+grb_trigger_time))
    data_tab.header.update('RA_OBJ',"%f" %RA)
    data_tab.header.update('DEC_OBJ',"%f" %DEC)

    data_tab.header.update('HDUCLASS','OGIP')
    data_tab.header.update('HDUCLAS1','SPECTRUM')
    data_tab.header.update('HDUVERS','1.2.1')    
    data_tab.header.update('HDUCLAS2','BKG')
    data_tab.header.update('HDUCLAS3','RATE')
    data_tab.header.update('HDUCLAS4','TYPE:II')
    data_tab.header.update('DETCHANS',ebins)
    data_tab.header.update('TLMIN1',1)
    data_tab.header.update('TLMAX1',ebins)
    data_tab.header.update('TRIGTIME',grb_trigger_time)
    
    #add GTI
    start_col=pyfits.Column(name="START",format="D",array=[grb_trigger_time+start])
    stop_col=pyfits.Column(name="STOP",format="D",array=[grb_trigger_time+stop])
    gti_tab=pyfits.new_table(pyfits.ColDefs([start_col,stop_col]))
    gti_tab.header.update('EXTNAME','GTI')
    gti_tab.header.update('TIMESYS','TT')
    gti_tab.header.update('TIMEREF','LOCAL')
    gti_tab.header.update('TIMEZERO','0')
    gti_tab.header.update('TIMEUNIT','s')


    #write!
    hdulist=pyfits.HDUList([hdu,data_tab,ebounds_tab,gti_tab])    
    hdulist.writeto("%s/%s" %(output_path,pha_filename))
    
    return ("%s/%s" %(output_path,pha_filename))
