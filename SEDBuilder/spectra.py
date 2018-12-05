from __future__ import division, print_function
import numpy as np
import astroquery.vizier
import astroquery.simbad
import astroquery.gaia
import astroquery.irsa
import astropy
import astropy.units as u

from pathlib import Path
import datetime

astroquery.vizier.Vizier.ROW_LIMIT = -1

# IRSA gives annoying warnings, this will disable all warnings (so should be commented out for debugging)
#import warnings
#warnings.filterwarnings("ignore")

h=6.63e-34
k=1.38e-23
c=3e8

def getInfo(ra=-1,dec=-1,name=-1):
    simbad=astroquery.simbad.Simbad()
    simbad.add_votable_fields('sp')
    simbad.add_votable_fields('sp_bibcode')
    simbad.add_votable_fields('ids')
    simbad.add_votable_fields('otype')
    simbad.add_votable_fields('plx')
    simbad.add_votable_fields('plx_bibcode')
    if (name!=-1):
        star=simbad.query_object(name)[0]
    elif (ra!=-1) & (dec!=-1):
        coord=astropy.coordinates.SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle, u.deg))
        star=simbad.query_region(coord,radius=100*u.arcsecond)
        #finds nearest object with a 2MASS name (in case something spurious is closer to given coords)
        whichEntry=0
        while 1-('2MASS' in star[whichEntry]['IDS'].decode('utf-8')):
            ids=star[whichEntry]['IDS'].decode('utf-8')
            whichEntry+=1
        star=star[whichEntry]
    else:
        print('\n Must enter name or coords (ra, dec) of star')

    
    ra=star['RA']
    dec=star['DEC']
    names=star['IDS'].decode('utf-8')
    startPoint = names.find('2MASS ')+len('2MASS ')
    startName = names[startPoint:]
    endPoint = startName.find('|')
    if endPoint==-1:
        tmName=startName
    else:
        tmName = startName[:endPoint]
    if (tmName[0]!='J') & (name!=-1): #if this named star doesn't have a 2mass name, looks for the nearest one that does
        print('\n No 2MASS entry for star: ',mainName)
        print('Looking for nearest object with a 2MASS name')
        return getInfo(ra=ra,dec=dec)
    
    starDict={}
    
    starDict['SimbadName']=star['MAIN_ID'].decode('utf-8') #getting rid of small b...
    starDict['2MASSID']=tmName
    starDict['ObjectType']=star['OTYPE'].decode('utf-8')
    starDict['StellarType']=star['SP_TYPE'].decode('utf-8')
    starDict['StellarTypeSource']=getRef(star['SP_BIBCODE'].decode('utf-8'))
    starDict['RA']=ra
    starDict['DEC']=dec
    starDict['CoordSource']=getRef(star['COO_BIBCODE'].decode('utf-8'))
    if np.ma.is_masked(star['PLX_VALUE']): # No parralax data for this star
        starDict['Distance']=-1
        starDict['DistanceSource']=''
    else:
        starDict['Distance']=1000/star['PLX_VALUE'] #1000 as Simbad gives value in mas
        starDict['DistanceSource']=getRef(star['PLX_BIBCODE'].decode('utf-8'))
    
    Teff,rad,lum,Ag=gaiaData(ra,dec)
    gaiaSource="Gaia Collaboration 2018"
    starDict['Teff']=Teff
    starDict['TeffSource']=gaiaSource
    if Teff<=0:
        starDict['TeffSource']=""
    starDict['Radius']=rad
    starDict['RadiusSource']=gaiaSource
    if rad<=0:
        starDict['RadiusSource']=""
    starDict['Luminosity']=lum
    starDict['LuminositySource']=gaiaSource
    if lum<=0:
        starDict['LuminositySource']=""
    starDict['Ag']=Ag
    starDict['AgSource']=gaiaSource
    if Ag<=0:
        starDict['AgSource']=""
    
    return starDict

def gaiaData(ra,dec):
    coord=astropy.coordinates.SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle, u.degree))
    width = u.Quantity(5, u.arcsecond)
    data=astroquery.gaia.Gaia.query_object(coordinate=coord, width=width, height=width)
    if len(data)==0:
        return -1,-1,-1,-1
    Teff=data['teff_val'][0]
    if np.ma.is_masked(Teff):
        Teff=-1
    rad=data['radius_val'][0]
    if np.ma.is_masked(rad):
        rad=-1
    lum=data['lum_val'][0]
    if np.ma.is_masked(lum):
        lum=-1
    Ag=data['a_g_val'][0]
    if np.ma.is_masked(Ag):
        Ag=-1
    return Teff,rad,lum,Ag

def getCoords(name):
    simbad=astroquery.simbad.Simbad()
    star=simbad.query_object(name)[0]
    return star['RA'],star['DEC']

def getRef(ref): #finds the reference (in the form "Herczeg+ 2014") for a given bibcode (in form "2014ApJ...786...97H" or "J/ApJ/786/97")
    catalogs=astroquery.vizier.Vizier.find_catalogs(ref)
    if len(catalogs)==0: # Seems some stars have data in Simbad but not Vizier - should only appear in metadata
        return ''
    description=catalogs[list(catalogs.items())[0][0]].description
    reference=description[description.rfind('(')+1:description.rfind(')')]
    return reference.replace(",","") #removes any commas
    
def getVizierSED(ra,dec,windowSize=2):
    #star=astroquery.simbad.Simbad.query_object(mainName)
    #starDict=getInfo(name=name,ra=ra,dec=dec)
    raString=str(ra).replace(' ','+').strip()
    decString=str(dec).replace(' ','+').strip()
    
    url="http://vizier.u-strasbg.fr/viz-bin/sed?-c="+raString+decString+"&-c.r="+str(windowSize)+"&-c.u=arcsec"
    data=astropy.table.Table.read(url)
    
    lambdas=c/np.array(1e9*data['sed_freq']) # converting to wavelength in m
    fluxs=1e-26*c*np.array(data['sed_flux'])/lambdas**2 # converting to F_lambda in W m
    sigmas=np.nan_to_num(1e-26*c*np.array(data['sed_eflux'])/lambdas**2) # converting to error in F_lambda in W m

    source=data['sed_filter']
    tables=data['_tabname']

    order=np.argsort(lambdas)

    lambdas=lambdas[order]
    fluxs=fluxs[order]
    sigmas=sigmas[order]
    source=source[order]
    tables=tables[order]
    
    sources=np.empty_like(lambdas,dtype=object)
    telescopes=np.empty_like(lambdas,dtype=object)
    
    for i,table in enumerate(tables):
        tableName = table[:table.rfind('/')]
        catalogs=astroquery.vizier.Vizier.find_catalogs(tableName)
        if len(catalogs)==0:
            continue #for some reason can't find table of data...
        description=catalogs[tableName].description
        reference=description[description.rfind('(')+1:description.rfind(')')].replace(',','')
        #print('ref: ',reference)
        #sources[i]='"'+reference.replace(',','')+'"'
        sources[i]=reference
        
        # Manually maps each data point to the telescope it comes from (probably should tabulate this elsewhere and work from that...)
        if ('Herschel' in source[i]) | (reference == 'Marsh+ 2016'):
            telescopes[i]='Herschel'
        elif 'Spitzer' in source[i]:
            telescopes[i]='Spitzer'
        elif 'WISE' in source[i]:
            telescopes[i]='WISE'
        elif ('2MASS' in source[i]) | ('Johnson' in source[i]):
            telescopes[i]='2MASS'
        elif ('GAIA' in source[i]) | ('Gaia' in source[i]):
            telescopes[i]='Gaia'
        elif 'PAN' in source[i]:
            telescopes[i]='Pan-STARRS'
        elif ('IRAS' in source[i]) | (reference == 'Saunders+ 2000'):
            telescopes[i]='IRAS'
        elif 'POSS' in source[i]:
            telescopes[i]='POSS'
        elif 'AKARI' in source[i]:
            telescopes[i]='AKARI'
        elif 'SDSS' in source[i]:
            telescopes[i]='SDSS'
        elif ('SMA' in source[i]) | (reference == 'Andrews+ 2013') | (reference == 'Harris+ 2012'):
            telescopes[i]='SMA'
        elif 'Cousins' in source[i]:
            telescopes[i]='Cousins'
        elif ('ALMA' in source[i]) | (reference == 'Pascucci+ 2016'):
            telescopes[i]='ALMA'
        elif 'ISO' in source[i]:
            telescopes[i]='ISO'
        elif 'DENIS' in source[i]:
            telescopes[i]='DENIS'
        elif 'HIP' in source[i]:
            telescopes[i]='Hipparcos'
        elif 'GALEX' in source[i]:
            telescopes[i]='GALEX'
        elif 'XMM' in source[i]:
            telescopes[i]='XMM'
        elif ('SCUBA' in source[i]) | (reference == 'Mohanty+ 2013'):
            telescopes[i]='SCUBA'
        elif 'UKIDSS' in source[i]:
            telescopes[i]='UKIDSS'
        elif 'MKO' in source[i]:
            telescopes[i]='MKO'
        elif 'MSX' in source[i]:
            telescopes[i]='MSX'
        elif (reference == 'Dzib+ 2013') | (reference == 'Dzib+ 2015'):
            telescopes[i]='VLA'
        elif (reference == 'Belloche+ 2011'):
            telescopes[i]='APEX'
        else:
            telescopes[i]='Unspecified'
            if source[i][0]!=':':
                print('___did you know about the ',source[i],' telescope?')
                
        if sources[i]=='Meng+ 2017': #points from this paper seem completely wrong? (https://ui.adsabs.harvard.edu/#abs/2017ApJ...836...34M/abstract)
            telescopes[i]='Unspecified'
        #telescopes[i]='"'+telescopes[i]+'"'
    return lambdas,fluxs,sigmas,sources,telescopes

def queryIrsa(ra,dec,windowSize=2): #at the moment only queries the Herchel point source catalogs, but open to suggestions
    wavelengths=[70,100,160,250,350,500]
    coord=astropy.coordinates.SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle, u.deg))
    ls=[]
    fs=[]
    es=[]
    for wl in wavelengths:
        if wl<200:
            catalog='ppsc_'+str(wl)
        else:
            catalog='spsc'+str(wl)
        table=astroquery.irsa.Irsa.query_region(coord,catalog=catalog,radius=windowSize*u.arcsecond,verbose=False)
        if len(table)==0:
            continue
        else:
            flux=table['flux'][0]
            error=flux/table['snr'][0]
            if 3*error > flux: # definition of upper limit
                flux=3*error
                error=0
            ls.append(wl*1.0e-6)
            fs.append(flux*c*1e-17/wl**2) #converting to F_lambda in SI
            es.append(error*c*1e-17/wl**2)
    return np.array(ls),np.array(fs),np.array(es)

#def getIrsaSED(tmName):
#    mainName,tmName,ra,dec,spType,objectType=getInfo(name='2MASS '+tmName)
#    irsaLs,irsaFs,irsaEs=queryIrsa(ra,dec)
#    if irsaLs[irsaLs<2e-4].size>0:
#        addToSED(tmName,irsaLs[irsaLs<2e-4],irsaFs[irsaLs<2e-4],irsaEs[irsaLs<2e-4],2017,'Marton+ 2017','Herschel')
#    if irsaLs[irsaLs>2e-4].size>0:
#        addToSED(tmName,irsaLs[irsaLs>2e-4],irsaFs[irsaLs>2e-4],irsaEs[irsaLs>2e-4],2017,'Schulz+ 2017','Herschel')

# Finds and returns all SED data readable online (from Vizier & IRAS) for a given object, can also save it to file        
def getSED(name=-1,ra=-1,dec=-1,saveDir=-1,windowSize=2,overwrite=0):
    starDict=getInfo(name=name,ra=ra,dec=dec)
    ra=starDict['RA']
    dec=starDict['DEC']
    
    if (saveDir!=-1) & (overwrite!=1):
        if saveDir[-1]!='/':
            saveDir=saveDir+'/'
        fName=saveDir+starDict['2MASSID']+'.ecsv'
        if Path(fName).is_file():
            if overwrite==0:
                print('\n Data file for this star already exists')
                print('You can find it at:')
                print(fName)
                print('Remember, all stars must have a 2MASS ID to be read, otherwise will default to nearest star in 2MASS catalog')
                print('The 2MASS ID of this star is: ',starDict['2MASSID'])
                print('Also known as ',starDict['SimbadName'],' at coordinates ',starDict['RA'],', ',starDict['DEC'])
                print('Returning the data from that file')
                print('To overwrite the file function use argument overwrite = 1 (default = 0)')
                print('Or to supress this warning set overwrite = -1')
            #table=astropy.io.ascii.read(fName)
            return getSEDFromFile(starDict['2MASSID'],saveDir)
    
    columnNames=('lambda', 'flux', 'error','source','telescope')
    #HOW CAN I SAVE AS SCIENTIFIC NOTATION FOR A FIXED NUMBER OF SIGNIFICANT FIGURES???
    dataTypes=('f4','f4','f4','object','object') # note: strings must be treated as objects not strings to allow variable length 
    vLs,vFs,vEs,vSs,vTs = getVizierSED(ra,dec,windowSize=windowSize) #SED data from vizier
    table=astropy.table.Table([vLs,vFs,vEs,vSs,vTs],names=columnNames,dtype=dataTypes)
    
    iLs,iFs,iEs = queryIrsa(ra,dec,windowSize=windowSize)
    for i in range(iLs.size):
        telescope='Herschel'
        if iLs[i]>2e-4:
            source='Schulz+ 2017'
        else:
            source='Marton+ 2017'
        table.add_row([iLs[i],iFs[i],iEs[i],source,telescope])
    table.sort('lambda')
    if saveDir!=-1:
        starDict=addPropertyToFile(starDict)
    
    for key in starDict.keys():
        table.meta[key]=str(starDict[key])
    comments=['All units for flux+wavelength are SI, I leave it to the user to convert/add astropy units',
    'Meta data about the star is stored under indivdual fields',
    'e.g. if data stored in variable called "dataTable" the R.A. of the star can be found via "dataTable.meta["RA"]".',
    'Everything intended to be read into astropy tables - either directly or via the getSEDFromFile() function.',
    'See https://github.com/zpenoyre/SEDBuilder for more details.',
    'Please cite Penoyre+2018 if you make use of this tool or data.']
    table.meta['comments']=comments
    
    if saveDir!=-1:
        table.meta['FileCreated']=str(datetime.datetime.today()).split()[0]
        if saveDir[-1]!='/':
            saveDir=saveDir+'/'
        fName=saveDir+starDict['2MASSID']+'.ecsv'
        astropy.io.ascii.write(table,output=fName, format='ecsv',overwrite=True)
    return table
    
# There are some supplementary properties we'd like to keep track of but must add by hand (e.g. Mdot)
# This function just helps reserve a space for them in the data file
def addPropertyToFile(starDict):
    properties=['Av','Age','Mass','Mdot','DiskMass','DiskRadius','BinaryFlag'] #would love to add more binary details but not immediately obvious how to
    starDict['Region']='' # May want to add this by hand (though obviously not very well defined)
    starDict['RegionSource']='' # Can be used to keep track of which works have led you to a particular star
    for field in properties:
        starDict[field]=-1
        starDict[field+'Source']=''
    starDict['ExtraField']='' #adding a spare field for whatever data a user might want to put in
    return starDict
    
# Retrieves an SED from file
def getSEDFromFile(twoMassID,saveDir):
    if saveDir[-1]!='/':
        saveDir=saveDir+'/'
    fName=saveDir+twoMassID+'.ecsv'
    if Path(fName).is_file():
        table=astropy.io.ascii.read(fName)
        #unfortunately all the meta data is converted to strings on saving, translating appropriate data back here
        table.meta['Distance']=float(table.meta['Distance'])
        table.meta['Teff']=float(table.meta['Teff'])
        table.meta['Radius']=float(table.meta['Radius'])
        table.meta['Luminosity']=float(table.meta['Luminosity'])
        table.meta['Ag']=float(table.meta['Ag']) # note - if value is known we recommend you use Av over Ag, but Ag comes from Gaia so is readily available for many stars
        table.meta['Av']=float(table.meta['Av'])
        table.meta['Age']=float(table.meta['Age'])
        table.meta['Mass']=float(table.meta['Mass'])
        table.meta['Mdot']=float(table.meta['Mdot'])
        table.meta['DiskMass']=float(table.meta['DiskMass'])
        table.meta['DiskRadius']=float(table.meta['DiskRadius'])
        table.meta['BinaryFlag']=int(table.meta['BinaryFlag'])
        return table
    else:
        print('\n No file found at: ')
        print(fName)
        return -1
        
# Adds new data points to an existing SED, tries not to duplicate anything
def addToSED(twoMassID,saveDir,ls,fs,es,source,telescope,overwrite=0):
    table=getSEDFromFile(twoMassID,saveDir)
    if (source in np.unique(table['source'])) & (overwrite!=1):
        print('\n Already data in table from this source (',source,')')
        print("We assume you don't want to duplicate data points so just returning the table saved to file")
        print('If you really do want to do this you can run the function with overwrite=1')
        return table
    for i in range(ls.size):
        table.add_row([ls[i],fs[i],es[i],source,telescope])
    table.sort('lambda')
    if overwrite!=-1: # expect to save to file (unless the user REALLY doesn't want to)
        if saveDir[-1]!='/':
            saveDir=saveDir+'/'
        fName=saveDir+twoMassID+'.ecsv'
        astropy.io.ascii.write(table,output=fName, format='ecsv')
    return table
    
def addToMetadata(twoMassID,saveDir,field,value,source,overwrite=0):
    if type(field)==list: # If given a list of properties
        nProperty=len(field)
        nValue=len(value)
        nSource=1
        if type(source)==list:
            nSource=len(source)
        if nProperty!=nValue:
            print('\n If supplying a list of properties then must give a matching length list of values')
        if (nSource!=1) & (nSource!=nProperty):
            print('\n Must either supply a single source or one for each property')
        for i in range(len(field)):
            if nSource!=1:
                thisSource=source[i]
            else:
                thisSource=source
            table=addToMetadata(twoMassID,saveDir,field[i],value[i],thisSource,overwrite=overwrite)
        return table
    if np.ma.is_masked(value):
        print('\n Trying to record a masked (i.e. not valid) value of ',field)
        print('Returning original table instead')
        return getSEDFromFile(twoMassID,saveDir)
    table=getSEDFromFile(twoMassID,saveDir)
    if (field not in table.meta.keys()) | ('ource' in field):
        print('\n No property called ',field,' in table')
        print('Valid properties stored in metadata are:')
        for key in table.meta.keys():
            if 'ource' not in key:
                print(key)
    if type(source)!='string':
        if len(source)<5:
            print('\n Must give a source of the form "Herczeg+ 2014", specifically must end with 4 digit year')
    year=int(source[-4:])
    oldSource=table.meta[field+'Source']
    if oldSource=='':
        oldYear=1632
    else:
        oldYear=int(oldSource[-4:])
    if (year<oldYear) & (overwrite!=1):
        if overwrite==0:
            print('\n Trying to add data about property: '+field+' from an older source')
            print('Existing value of ',table.meta[field],' comes from ',table.meta[field+'Source'])
            print('(Note - if exisiting data is from Gaia and not astrometric you may want to replace it, they are fitted to very few data points)')
            print('You can forcefully write to file by setting overwrite=1')
            print('For now returning original table')
            print('You can suppress this warning by setting overwrite=-1')
        return table
    table.meta[field]=value
    table.meta[field+'Source']=source
    
    if saveDir[-1]!='/':
        saveDir=saveDir+'/'
    fName=saveDir+twoMassID+'.ecsv'
    astropy.io.ascii.write(table,output=fName, format='ecsv') 
    
    return table
    
# Very similar to above, except now we may want multiple sources for one star (specifying that all these papers touch on this star)
def addRegionToMetadata(twoMassID,saveDir,region,source):
    table=getSEDFromFile(twoMassID,saveDir)
    if table.meta['Region']=='': # no record of this stars region at the moment
        table.meta['Region']=region
        table.meta['RegionSource']=source
    else:
        if source in table.meta['RegionSource']:
            print('\n Already recorded this source in addRegionToMetadata')
        else:
            table.meta['RegionSource']=table.meta['RegionSource']+', '+source
    if saveDir[-1]!='/':
        saveDir=saveDir+'/'
    fName=saveDir+twoMassID+'.ecsv'
    astropy.io.ascii.write(table,output=fName, format='ecsv') 
    return table
    
# Have an extra field! Write whatever you want! Go crazy!
# Note - I've sometimes used this to record the parameters of fitting models
# You can put lists or dictionaries in here, as long as you decide how to convert them to strings and back
def rewriteExtraField(twoMassID,saveDir,entry):
    table=getSEDFromFile(twoMassID,saveDir)
    table.meta['ExtraField']
    if saveDir[-1]!='/':
        saveDir=saveDir+'/'
    fName=saveDir+twoMassID+'.ecsv'
    astropy.io.ascii.write(table,output=fName, format='ecsv')
    return table
    
def plotSED(thisPlot,table,
colours=['#6699CC','#FFD23F','#FF8C42','#FF3C38','#A23E48'],
telescopes=['Gaia','2MASS','WISE','Spitzer','Herschel']):
    ls=table['lambda']
    fs=table['flux']
    es=table['error']
    ss=table['source']
    ts=table['telescope']
    for i in range(len(telescopes)):
        tel=telescopes[i]
        col=colours[i]
        points=np.flatnonzero((ts==tel) & (es>0))
        thisPlot.errorbar(ls[points],ls[points]*fs[points],yerr=ls[points]*es[points],color=col,fmt='o',lw=1,zorder=4,label=tel)
        uppers=np.flatnonzero((ts==tel) & (es<=0))
        thisPlot.scatter(ls[uppers],ls[uppers]*fs[uppers],facecolors='w',edgecolors=col,marker='v',zorder=3,label='')
    noTel=np.flatnonzero((~np.isin(ts,telescopes)) & (es>0))
    thisPlot.errorbar(ls[noTel],ls[noTel]*fs[noTel],yerr=ls[noTel]*es[noTel],color='grey',fmt='o',lw=1,zorder=2,label='Other')
    noTelUppers=np.flatnonzero((~np.isin(ts,telescopes)) & (es<=0))
    thisPlot.scatter(ls[noTelUppers],ls[noTelUppers]*fs[noTelUppers],facecolors='w',edgecolors='grey',marker='v',zorder=1,label='')
    
    thisPlot.set_yscale('log')
    thisPlot.set_ylim(0.1*np.min(ls*fs),10*np.max(ls*fs))
    thisPlot.set_xscale('log')
    thisPlot.set_xlim(0.5*np.min(ls),2*np.max(ls))
    
    thisPlot.set_ylabel(r'$\lambda F_\lambda$ (W m)')
    thisPlot.set_xlabel(r'$\lambda$ (m)')
    
    thisPlot.legend(title=table.meta['SimbadName'],frameon=False)
        
    