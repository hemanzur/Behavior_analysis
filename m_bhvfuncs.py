#bhvfuncs.py

import scipy.io as spio
import numpy as np
import os, sys, re, string, pdb
from glob import glob
from PyQt4 import QtGui

import pdb

import guidata.dataset.datatypes as dt
import guidata.dataset.dataitems as di
import guidata
app=guidata.qapplication()

'''
import __main__

if hasattr(__main__, 'BhvDir') or\
   hasattr(__main__, 'MPCDir') or\
   hasattr(__main__, 'ImgDir') or\
   hasattr(__main__, 'HomeDir'):
    BhvDir  = __main__.BhvDir
    MPCDir  = __main__.MPCDir
    ImgDir  = __main__.ImgDir
    HomeDir = __main__.HomeDir
else:
    class SelDir(dt.DataSet):
        """Select the directories to work with"""
        HomeDir = di.DirectoryItem(label='HomeDir', default = os.environ['HOME']).set_pos(col=0)
        MPCDir  = di.DirectoryItem(label='MPC Files Dir', default = os.environ['HOME']).set_pos(col=1)
        BhvDir  = di.DirectoryItem(label='Beh Files Dir', default = os.environ['HOME']).set_pos(col=0)
        ImgDir  = di.DirectoryItem(label='ImgDir', default = os.environ['HOME']).set_pos(col=0)

    SelDirs = SelDir()
    if SelDirs.edit()==1:
        HomeDir = SelDirs.HomeDir
        ImgDir  = SelDirs.ImgDir
        BhvDir  = SelDirs.BhvDir
        MPCDir  = SelDirs.MPCDir
'''
    
MSN1=['NP0','1T','ODD','Odd','1L']
MSN2=['LEFT','9k','6k','9K','6K','2T','UNBIAS','LeftSip','MAP2',
      'CentWhtR','12K','1K','15K','2K','13K']

########################################################################################################################

def mpc2beh(filesdir = '', savedir = '', overwrite = False, calcParams = True):

    '''Routine to transform the MPC files into matlab files\n
    Input:
            *filesdir=can be a direcotyr or a filename
            *savedir= must be a directory
            *overwrite= whether to overwrite the existing file or not
    Output: data.
    '''

    if not filesdir or not os.path.isdir(filesdir):
        filesdir = str(QtGui.QFileDialog.getExistingDirectory(caption = 'MPC Files Dir'))
        if not filesdir: return
        
    if not savedir or not os.path.isdir(savedir):
        savedir = str(QtGui.QFileDialog.getExistingDirectory(caption = 'Save dir'))
        if not savedir: return
        
    if os.path.isdir(filesdir):
        files = glob(os.path.join(filesdir, '!*.Subject*'))

        if files: files.sort()
        else:
            print 'Ther were no files in here !'
            return
        
        for j,k in enumerate(files):
            filesdir, files[j] = os.path.split(k)
            
    elif os.path.isfile(filesdir):
        if not re.search('!.*.Subject .*', filesdir):
            raise SystemExit('It seems that this is not an MPC file !')
        filesdir, files = os.path.split(filesdir)
        files = [files]

    for k in files:

        YY  = k[3:5]
        MM  = k[6:8]
        DD  = k[9:11]
        HH  = k[12:14]
        MI  = k[15:17]
        RAT = k[k.find('Subject')+8:]

        filename = '%s_%s%s%s_%s%s*.beh' % (RAT,YY,MM,DD,HH,MI)

        if glob(os.path.join(savedir,filename))!=[]\
           and overwrite==0 \
           or not os.path.isfile(os.path.join(filesdir,k)):
            continue
        else:
            print 'Processing... %s' % k,

        fid = open(os.path.join(filesdir, k), 'rU')
        w   = fid.readlines()
        fid.close()

        Data={}

        for m, n in enumerate(w):
            if n == '\n': continue

            st    = n.split(': ')
            st[0] = st[0].replace(' ','')
            st[1] = st[1].replace('\n','')

            if st[0] in ['File','Subject','MSN']:
                Data[st[0]] = str(st[1].strip())
            elif st[0].find('Date')!=-1 or st[0].find('Time')!=-1:
                Data[st[0]] = [int(x) for x in re.split('[:/]',st[1].strip())]
            elif st[0] in ['Box','Experiment','Group']:
                if len(st[1].strip()) > 0:
                    Data[st[0]] = int(st[1].strip())
                else:
                    Data[st[0]] = 0

            if st[0] == 'MSN': break

        for l,n in enumerate(w[m+1:]):
            if n.find('Start') != -1 or n == '\n':
                break
            e = n.strip().replace('\n','')
            e = e.split(':')

            if e[0] in string.uppercase:
                curvar = e[0]
                Data[curvar] = []
                continue
            else:
                e = e[1].strip()
                e = re.split(' +',e)
                e = [float(x) for x in e]
                if curvar == 'X' and e == [0,0,0,0,0]:
                    break
                else:
                    Data[curvar]=Data[curvar]+e

        if  GetMapping(Data) in [1,2]:

            TableEvent = BehMapping(Data)
            TableEventName = TableEvent[:,0]
            TableEventCode = TableEvent[:,1]

            Data['X'] = np.array(Data['X'])
            Data['X'] = Data['X'][np.nonzero(Data['X'])]
            ECodes = np.floor(Data['X']/1000000)
            Data['X'] = Data['X']-ECodes*1000000
            X = np.unique(ECodes)
            X.sort()
            Data['X'] = Data['X']-Data['X'][0]

            #create a np.object to be saved as a cell array in matlab
            Data['EventTS']   = []
            Data['EventCode'] = []
            Data['EventName'] = []

            for l in X:
                Data['EventTS'].append(Data['X'][ECodes == l])
                Data['EventCode'].append(int(l))
                Data['EventName'].append(TableEventName[l == TableEventCode][0])

            Data['EventTS']   = np.array(Data['EventTS']  , dtype = np.object, ndmin=1)
            Data['EventCode'] = np.array(Data['EventCode'], dtype = np.object, ndmin=1)
            Data['EventName'] = np.array(Data['EventName'], dtype = np.object, ndmin=1)
            Data.pop('X')

            if calcParams: Data = GetBhvParams(Data)

        savefile = os.path.join(savedir, filename.replace('*','_'+Data['MSN']))

        if os.path.isfile(savefile):
            print '\t...deleting'
            os.remove(savefile)
        else: print

        spio.savemat(savefile, {'Data':Data}, format = '5',
                     appendmat=False, oned_as='row')

########################################################################################################################

def GetMapping(Data):

    MSN = str(Data['MSN'])
    if [x for x in MSN1+MSN2 if MSN.find(x)!=-1]:
        if [x for x in MSN1 if MSN.find(x)!=-1] and not [x for x in MSN2 if MSN.find(x)!=-1]:
            Mapping=1
        elif [x for x in MSN2 if MSN.find(x)!=-1]:
            Mapping=2
        return Mapping
    else:
        return -1

########################################################################################################################

def BehMapping(Data=None):

    TableEvent1 = np.array([['SessionStart', 17],
            ['IRLickOn', 21],
            ['IRLickOff', 22],
            ['IRCTROn', 23],
            ['IRCTROff', 24],
            ['IRRTOn', 25],
            ['IRRTOff', 26],
            ['IRLTOn', 27],
            ['IRLTOff', 28],
            ['OUT1HL', 31],
            ['OUT2SL', 32],
            ['OUT3CTRRED', 33],
            ['OUT4CTRWHTOn', 34],
            ['OUT5RTRED', 35],
            ['OUT6RTGRN', 36],
            ['OUT7RTYLW', 37],
            ['OUT8EMPTY', 38],
            ['OUT9LTRED', 39],
            ['OUT10LTGRN', 40],
            ['OUT11LTYLW', 41],
            ['OUT12NOISE', 42],
            ['OUT13SOL1', 43],
            ['OUT14SOL2', 44],
            ['OUT15SOL3', 45],
            ['OUT16SOL4',46],
            ['SOUND1', 47],
            ['SOUND2', 48],
            ['SOUND3', 49],
            ['OUT4CTRWHTOff',50]], dtype = np.object)

    TableEvent2=np.array([['SessionStart',17],
            ['RightLickOn',21],
            ['RightlickOff',22],
            ['RightPokeOn',23],
            ['RightPokeOff',24],
            ['CentPokeOn',25],
            ['CentPokeOff',26],
            ['LeftLickOn',27],
            ['LeftLickOff',28],
            ['LeftPokeOn',29],
            ['LeftPokeOff',30],
            ['HouseLightOn',31],
            ['HouseLightOff',32],
            ['RightSipperLightOn',33],
            ['RightSipperLightOff',34],
            ['RedFrontLightOn',35],
            ['RedFrontLightOff',36],
            ['WhiteFrontLightOn',37],
            ['WhiteFrontLightOff',38],
            ['NPRedLightOn',39],
            ['NPRedLightOff',40],
            ['NPGreenLightOn',41],
            ['NPGreenLightOff',42],
            ['NPYellowLightOn',43],
            ['NPYellowLightOff',44],
            ['LeftSipperLightOn',45],
            ['LeftSipperLightOff',46],
            ['NoiseOn',47],
            ['NoiseOff',48],
            ['Solnd1',49],
            ['Solnd2',50],
            ['Solnd3',51],
            ['Solnd4',52],
            ['Sound1',53],
            ['Sound2',54],
            ['Catch',55]],dtype=np.object)

    if Data:
        if GetMapping(Data)==1:
            return TableEvent1
        elif GetMapping(Data)==2:
            return TableEvent2

    else:
        return TableEvent1, TableEvent2

########################################################################################################################

def GetBhvParams(Data):
    '''Get behavioral parameters\n
    Input: behavioral Data dictionary or file\n
    Output: behavioral data dict with all behavioral parameters calculated'''

    try:
        if os.path.isfile(Data):
            Data=LoadBehFile(Data)
    except TypeError:
        if type(Data)!=dict:
            raise NameError('Data is neither file nor dictionary !')

    #===========================================================================================================#
        
    if GetMapping(Data)==1:
        
        EvtName1=['SOUND1','SOUND2','IRRTOn','IRRTOff','IRCTROn','OUT4CTRWHTOn','OUT13SOL1','IRLickOn']
        EvtName2=['Tone1', 'Tone2', 'NpIn',  'NpOut',  'RpIn',   'CLight',      'Solnd1',   'Lick']
        Vars={}
        for k,l in zip(EvtName1,EvtName2):
            indx = np.flatnonzero(Data['EventName']==k)
            if indx:
                Vars[l] = Data['EventTS'][indx][0]
        
        Stims=['Tone1','Tone2','CLight']
        
        for j,k in enumerate(Stims):
            if k in Vars.keys() and 'Lick' in Vars.keys():

                CurStim='Stim'+str(j)
                Data[CurStim]={}
                Stim = Vars[k]
                # First find the hits and the misses
                # Then find the closest event to each one of the Hits
                HitsParams = GetHits(Stim, Vars['Lick'])
                RTT = HitsParams['ThirdLickHitTS']-HitsParams['StimHitsTS']
                RT4, _ = DistXY(HitsParams['StimHitsTS'], Vars['Lick'])
                RT4 = RTT - RT4

                if 'NpIn' in Vars.keys():
                    RT0, RT0Indx = DistYX( HitsParams['StimHitsTS'], Vars['NpIn'] )
                    Data[CurStim]['RT0'] = RT0

                if 'NpOut' in Vars.keys():
                    RT1, RT1Indx = DistXY( HitsParams['StimHitsTS'], Vars['NpOut'] )
                    Data[CurStim]['RT1'] = RT1

                if  ('RpIn' in Vars.keys()) and ('NpOut' in Vars.keys()):
                    RT2, RT2Indx = DistXY(Vars['NpOut'][RT1Indx], Vars['RpIn'])
                    RT3, RT3Indx = DistXY(Vars['RpIn'][RT2Indx], Vars['Lick'])
                    Data[CurStim]['RT2'] = RT2
                    Data[CurStim]['RT3'] = RT3
                    Data[CurStim]['RT4'] = RT4

                Data[CurStim]['Descr']    = str(k)
                Data[CurStim]['HitsTS']   = HitsParams['StimHitsTS']
                Data[CurStim]['HitsIndx'] = HitsParams['StimHitsIndx']
                Data[CurStim]['MissTS']   = HitsParams['StimMissTS']
                Data[CurStim]['MissIndx'] = HitsParams['StimMissIndx']
                Data[CurStim]['RTT']      = RTT
                Data[CurStim]['StimTS']   = Stim
                Data[CurStim]['RT4']      = HitsParams['ThirdLickHitTS'] - HitsParams['FirstLickHitTS']

                #Copy all the time stamps of the events to make easier calculations for rasters
                for x in ['NpIn', 'NpOut', 'RpIn', 'Solnd1', 'Lick']:
                    if Vars.has_key(x):
                        Data[CurStim][x] = Vars[x]
                if Vars.has_key('Solnd1'):
                    Data[CurStim]['Solnd'] = Vars['Solnd1']

    #===========================================================================================================#
                    
    elif GetMapping(Data)==2:
        
        #Obtain Timestamps of events
        TableEvent = BehMapping(Data)
        EvtCodes   = np.array([33,45,37,35,53,54,55,26,25,23,29,21,27,49,50,51,52])
        indx       = [np.flatnonzero(TableEvent[:,1]==k)[0] for k in EvtCodes if any(k==TableEvent[:,1])]
        EvtNames   = ['RSipLight','LSipLight','WhtFLight','RedFLight','Tone1','Tone2','Catch',
                      'NpOut','NpIn','RRpIn','LRpIn','RLick','LLick','Solnd1','Solnd2','Solnd3','Solnd4']

        Vars={}
        for j,k in zip(TableEvent[indx,0],EvtNames):
            if np.any(Data['EventName']==j):
                Vars[k] = Data['EventTS'][Data['EventName']==j][0]
                Vars[k] = np.array(Vars[k], ndmin=1)

        #Check the association between stimuli and reward
        if Data.has_key('S'):           #If data has the 'S' field it builds an array of Stim - Resp association
            Stims = np.array([['Tone1',    53],
                              ['Tone2',    54],
                              ['Catch',    55],
                              ['RSipLight',33],
                              ['LSipLight',45],
                              ['RedFLight',35],
                              ['WhtFLight',37]], dtype = np.object)

            StimResp = np.array(np.zeros((0,4)), dtype = np.object)

            for k in range(11,20):
                if np.flatnonzero(Data['S'][k]==Stims[:,1]).size>0:
                    if Data['S'][k-10]==1:
                        CLick='RLick'
                        ILick='LLick'
                    elif Data['S'][k-10]==2:
                        CLick='LLick'
                        ILick='RLick'
                    elif Data['S'][k-10]==3:
                        CLick='Miss'
                        ILick='Miss'

                    s=np.flatnonzero(Data['S'][k]==Stims[:,1])[0]
                    StimResp=np.append(StimResp,[[Stims[s,0],Stims[s,1],CLick,ILick]],axis=0)
        
        else:                           #If not, makes one by default
            StimResp = np.array([['Tone1',53,'RLick','LLick'],
                                 ['Tone2',54,'LLick','RLick'],
                                 ['Catch',55,'Miss' ,'Miss' ]],dtype=np.object)
        
        for j, k in enumerate(StimResp):    #main loop that extracts and calculates all the parameters for a given stimuli
            CurStim = 'Stim' + str(j)
            Data[CurStim] = {}

            if Vars.has_key(k[0]) and Vars.has_key(k[2]): #check whether the StimTS and the lickTS are present
                Stim  = Vars[k[0]]
                CLick = Vars[k[2]]

                # First find the hits and the misses
                HitsParams = GetHits(Stim, CLick)
                
                #Routine to calculate the incorrect responses parameters
                #Code to find the true Hits
                #strg='6K[R,L]_9K[R,L]|6K[R,L]_12K[R,L]|UNBIAS|2T|LeftSip|MAP2'
                
                if HitsParams['StimHitsTS'].size>0 and Vars.has_key(k[3]):
                    
                    ILick      = Vars[k[3]]
                    IncParams  = GetHits(Stim, ILick)
                    # Look for the Stim indices that are shared by Hits and Errors
                    indx = []       # HitsIndx that are equal to IncIndx 
                    for x in IncParams['StimHitsIndx']:
                        if np.flatnonzero( HitsParams['StimHitsIndx'] == x ) :
                            indx.append(x)
                    indx = np.array(indx)

                    # Get the total reaction time fot hits
                    RTT,  _ = DistXY( HitsParams['StimHitsTS'], HitsParams['ThirdLickHitTS'] )

                    # Get total reaction time for incorrects
                    if IncParams['StimHitsTS'].size > 0:
                        RTTe, _ = DistXY( IncParams ['StimHitsTS'], IncParams ['ThirdLickHitTS'] )

                        indx2 = []
                        if np.any(indx):
                            for y,x in enumerate(indx):
                                if not RTTe[IncParams['StimHitsIndx']==x] < RTT[HitsParams['StimHitsIndx']==x]:
                                    indx2.append(y)
                            indx = np.delete(indx, indx2)

                        if np.any(indx):
                            indx     = [np.flatnonzero(HitsParams['StimHitsIndx']==x)[0] for x in indx]
                            HitsParams['StimHitsIndx'] = np.delete( HitsParams['StimHitsIndx'], indx)
                            HitsParams['StimHitsTS'] = Stim[HitsParams['StimHitsIndx']]

                            HitsTmp = GetHits(HitsParams['StimHitsTS'], CLick)
                            HitsParams['FirstLickHitIndx'] = HitsTmp['FirstLickHitIndx']
                            HitsParams['ThirdLickHitIndx'] = HitsTmp['ThirdLickHitIndx']
                            HitsParams['FirstLickHitTS']   = HitsTmp['FirstLickHitTS']
                            HitsParams['ThirdLickHitTS']   = HitsTmp['ThirdLickHitTS']

                        # Calculate the true Misses
                        tmp = np.concatenate((HitsParams['StimHitsIndx'], IncParams['StimHitsIndx']))
                        MissIndx = np.array([o for o in range(len(Stim)) if o not in tmp])

                        if MissIndx.size > 0:
                            HitsParams['StimMissTS']   = Stim[MissIndx]
                            HitsParams['StimMissIndx'] = MissIndx
                        else:
                            HitsParams['StimMissTS']   = np.array([])
                            HitsParams['StimMissIndx'] = np.array([])
                        
                        if np.any(HitsParams['StimMissIndx']):
                            HitsParams['StimMissTS'] = Stim[HitsParams['StimMissIndx']]
                            
                        # Now that we know the true Incorrects, calculate parameters
                        if 'NpIn'  in Vars.keys():
                            RT0e, RT1Indx = DistYX(IncParams['StimHitsTS'], Vars['NpIn'])
                            
                        if 'NpOut' in Vars.keys():
                            RT1e, RT1Indx = DistXY(IncParams['StimHitsTS'], Vars['NpOut'])
        
                        if   k[3]=='LLick': RPort = 'LRpIn'
                        elif k[3]=='RLick': RPort = 'RRpIn'

                        if  RPort in Vars.keys() and 'NpOut' in Vars.keys():
                            RT2e, RT2Indx = DistXY( Vars['NpOut'][RT1Indx], Vars[RPort] )
                            if any(RT2e):
                                RT3e, RT3Indx = DistXY( Vars[RPort][RT2Indx], ILick )
                            else:
                                RT3e = np.array([])

                        RT4e    = IncParams['ThirdLickHitTS'] - IncParams['FirstLickHitTS']
                        
                    else:
                        RTTe = np.array([])
                        RT0e = np.array([])
                        RT1e = np.array([])
                        RT2e = np.array([])
                        RT3e = np.array([])
                        RT4e = np.array([])

                    Data[CurStim]['ErrTS']   = IncParams['StimHitsTS']
                    Data[CurStim]['ErrIndx'] = IncParams['StimHitsIndx']
                    Data[CurStim]['RTTe']    = RTTe
                    Data[CurStim]['RT0e']    = RT0e
                    Data[CurStim]['RT1e']    = RT1e
                    Data[CurStim]['RT2e']    = RT2e
                    Data[CurStim]['RT3e']    = RT3e
                    Data[CurStim]['RT4e']    = RT4e

                # Recalculate the Hits Timestamps and the RTT for corrects with the updated HitsInx
                if HitsParams['StimHitsTS'].size>0:
                
                    RTT , _ = DistXY( HitsParams['StimHitsTS'], HitsParams['ThirdLickHitTS'] )
                    
                    if Vars.has_key('NpIn'):      #Calculate the foreperiod
                            RT0, RT0Indx = DistYX(HitsParams['StimHitsTS'], Vars['NpIn'])
                            
                    if Vars.has_key('NpOut'):      #Calculate RT1
                        RT1, RT1Indx = DistXY(HitsParams['StimHitsTS'], Vars['NpOut'])

                    if   k[2]=='LLick': RPort='LRpIn'
                    elif k[2]=='RLick': RPort='RRpIn'

                    if  RPort in Vars.keys() and 'NpOut' in Vars.keys():
                        RT2, RT2Indx = DistXY(Vars['NpOut'][RT1Indx], Vars[RPort])
                        if any(RT2):
                            RT3, RT3Indx = DistXY(Vars[RPort][RT2Indx], CLick)
                        else:
                            RT3 = np.array([])
                    else: RT2 = np.array([])
                else:
                    RTT = np.array([])
                    RT0 = np.array([])
                    RT1 = np.array([])
                    RT2 = np.array([])
                    RT3 = np.array([])
                    RT4 = np.array([])
                    
                Data[CurStim]['Descr']    = str(k[0])
                Data[CurStim]['HitsTS']   = HitsParams['StimHitsTS']
                Data[CurStim]['HitsIndx'] = HitsParams['StimHitsIndx']
                Data[CurStim]['MissTS']   = HitsParams['StimMissTS']
                Data[CurStim]['MissIndx'] = HitsParams['StimMissIndx']
                Data[CurStim]['StimTS']   = Vars[k[0]]

                if Vars.has_key('NpIn'):  Data[CurStim]['NpIn']  = Vars['NpIn']
                if Vars.has_key('NpOut'): Data[CurStim]['NpOut'] = Vars['NpOut']

                Data[CurStim]['RTT']      = RTT
                
                if 'RT0' in locals(): Data[CurStim]['RT0'] = RT0
                if 'RT1' in locals(): Data[CurStim]['RT1'] = RT1
                if 'RT2' in locals(): Data[CurStim]['RT2'] = RT2
                if 'RT3' in locals(): Data[CurStim]['RT3'] = RT3
                    
                Data[CurStim]['RT4']      = HitsParams['ThirdLickHitTS'] - HitsParams['FirstLickHitTS']
                
                if k[2] in Vars.keys(): Data[CurStim]['HLick'] = Vars[k[2]]
                if k[3] in Vars.keys(): Data[CurStim]['ELick'] = Vars[k[3]]

                if  k[2] == 'RLick':
                    Data[CurStim]['HRpIn'] = Vars['RRpIn']
                    if Vars.has_key('Solnd1'):   Data[CurStim]['Solnd'] = Vars['Solnd1']
                    if Vars.has_key('Solnd2'): Data[CurStim]['Solnd'] = Vars['Solnd2']
                        
                elif k[2]=='LLick':
                    Data[CurStim]['HRpIn'] = Vars['LRpIn']
                    if Vars.has_key('Solnd3'):   Data[CurStim]['Solnd'] = Vars['Solnd3']
                    if Vars.has_key('Solnd4'): Data[CurStim]['Solnd'] = Vars['Solnd4']
                
                if   k[3]=='LLick': Data[CurStim]['ERpIn'] = Vars['LRpIn']                    
                elif k[3]=='RLick': Data[CurStim]['ERpIn'] = Vars['RRpIn']


            #~Calculate the parameters for the catch trials
            elif k[0]=='Catch' and k[2]=='Miss' and Vars.has_key('Catch'):    
                Stim = Vars[k[0]]
                Lick = np.concatenate((Vars['LLick'],Vars['RLick']))
                Lick.sort()
                CatchParams = GetHits(Stim, Lick)
                
                if 'NpIn'  in Vars.keys(): RT0, RT1Indx = DistYX(CatchParams['StimHitsTS'], Vars['NpIn'])
                if 'NpOut' in Vars.keys(): RT1, RT1Indx = DistXY(CatchParams['StimHitsTS'], Vars['NpOut'])
                
                Vars['RpIn'] = np.concatenate( (Vars['RRpIn'], Vars['LRpIn']) )
                Vars['RpIn'].sort()

                if  'RpIn' in Vars.keys() and 'NpOut' in Vars.keys():
                    RT2, RT2Indx = DistXY( Vars['NpOut'][RT1Indx], Vars['RpIn'])
                    if any(RT2):
                        RT3, RT3Indx = DistXY(Vars['RpIn'][RT2Indx], Lick)
                    else:
                        RT3 = np.array([])
                        
                Data[CurStim]['Descr']    = k[0]
                Data[CurStim]['HitsTS']   = CatchParams['StimHitsTS']
                Data[CurStim]['HitsIndx'] = CatchParams['StimHitsIndx']
                Data[CurStim]['MissTS']   = CatchParams['StimMissTS']
                Data[CurStim]['MissIndx'] = CatchParams['StimMissIndx']
                Data[CurStim]['StimTS']   = Vars[k[0]]
                Data[CurStim]['RLick']    = Vars['RLick']
                Data[CurStim]['LLick']    = Vars['LLick']
                
                if Vars.has_key('RpIn'): Data[CurStim]['HRpIn'] = Vars['RpIn']
                
                Data[CurStim]['HLick'] = Lick
                
                if Vars.has_key('NpIn'):  Data[CurStim]['NpIn']  = Vars['NpIn']
                if Vars.has_key('NpOut'): Data[CurStim]['NpOut'] = Vars['NpOut']

                Data[CurStim]['RTT'] = RTT
                
                if 'RT0' in locals(): Data[CurStim]['RT0'] = RT0
                if 'RT1' in locals(): Data[CurStim]['RT1'] = RT1
                if 'RT2' in locals(): Data[CurStim]['RT2'] = RT2
                if 'RT3' in locals(): Data[CurStim]['RT3'] = RT3
                
                Data[CurStim]['RT4']      = []
                
                if CatchParams['StimHitsTS'].size>0 and CatchParams['ThirdLickHitTS'].size>0:
                    Data[CurStim]['RTT'], _ = DistXY( CatchParams['StimHitsTS'], CatchParams['ThirdLickHitTS'] )
                else:
                    Data[CurStim]['RTT'] = np.array([])
                
            if not Data[CurStim]:
                Data.pop(CurStim)

    return Data

########################################################################################################################

def GetHits(StimTS, LickTS, RespWin = 3.0, WetLick = 3):
    '''If len(StimTS)=m and len(LickTS)=n
    create a two matrices of m x n and (n x m)'''
    #pdb.set_trace()
    LickTS     = np.array(LickTS, ndmin = 2)
    StimTS     = np.array(StimTS, ndmin = 2)
    xx         = np.tile(StimTS, (LickTS.size,1) )
    yy         = np.tile(LickTS, (StimTS.size,1) ).transpose()
    Dif        = np.round(yy - xx, 3)
    Dif[Dif<0.00] = 1e6

    #Find the indices
    LickIndx   = Dif.argmin(0)[Dif.min(0)<1e6]

    Res = {}
    #Get the time of the third lick
    ValidStims    = np.flatnonzero( LickIndx+WetLick-1 <= LickTS.size-1 )
    ValidStims
    DryLicks      = Dif[LickIndx[ValidStims], ValidStims]
    WetLicksIndx  = LickIndx[ValidStims] + WetLick -1
    WetLicks      = Dif[WetLicksIndx, ValidStims]
    
    Res['StimHitsIndx']     = np.flatnonzero( WetLicks <= RespWin )
    Res['StimMissIndx']     = np.delete( range( len(StimTS) ), Res['StimHitsIndx'] )
    Res['StimHitsTS']       = StimTS[0, Res['StimHitsIndx'] ]
    Res['StimMissTS']       = StimTS[0, Res['StimMissIndx'] ]
    Res['FirstLickHitIndx'] = LickIndx[ValidStims][ Res['StimHitsIndx'] ]
    Res['ThirdLickHitIndx'] = WetLicksIndx[ Res['StimHitsIndx'] ]
    Res['FirstLickHitTS']   = LickTS[0, Res['FirstLickHitIndx'] ]
    Res['ThirdLickHitTS']   = LickTS[0, Res['ThirdLickHitIndx'] ]

    return Res

########################################################################################################################

def DistXY(x, y):
    ''' Obtain the minimum distances between the events in two vectors
    of timestamps of different lenght. Note: Make sure that y happens after x'''
    x   = np.array(x, ndmin=1)
    y   = np.array(y, ndmin=1)
    xx  = np.tile(x, (np.size(y), 1))
    yy  = np.tile(y, (np.size(x), 1)).transpose()
    Dif = np.round_(yy - xx, 3)
    Dif[Dif<0.00] = 1e6
    iDif = Dif.argmin(0)
    Dif  = Dif.min(0)
    e    = Dif != 1e6
    Dif  = Dif[e]
    iDif = iDif[e]

    return Dif, iDif

########################################################################################################################

def DistYX(x, y):
    ''' Obtain the distances between the events in two vectors
    of timestamps of different lenght. Note: Make sure that x happens after y'''
    x    = np.array(x, ndmin=1)
    y    = np.array(y, ndmin=1)
    xx   = np.tile(x,(np.size(y),1))
    yy   = np.tile(y,(np.size(x),1)).transpose()
    Dif  = np.round(xx - yy,3)
    Dif[Dif<0.00] = 1e6
    iDif = Dif.argmin(0)
    Dif  = Dif.min(0)
    e    = Dif != 1e6
    Dif  = Dif[e]
    iDif = iDif[e]

    return Dif, iDif

########################################################################################################################

def SparseDistance(x, y, direction = 'xy', maxTime = 1e6):
    '''Sparse calculation of minimum distance between two vectors
    of different length.
    Inputs: x,y : vectors of timestamps of different length
            direction: "xy" if y happens after x
                       "yx" x happens after y
            maxTime: maximum time lag between the two events
    Output: touple with the distances and the index of the longest vector
            that give those distances'''
    
    x   = np.array(x, ndmin = 2)
    y   = np.array(y, ndmin = 2)
    xx  = np.tile(x, (y.size, 1))
    yy  = np.tile(y, (x.size, 1)).transpose()

    if direction == 'xy':
        Dif = np.round(yy - xx, 3)
    elif direction == 'yx':
        Dif = np.round(xx - yy, 3)

    Dif[Dif < 0.00] = 1e6
    xIndxDif = np.unique(Dif.argmin(axis = 1))
    yIndxDif = np.unique(Dif.argmin(axis = 0))
    Dif      = Dif[yIndxDif, xIndxDif]
    indx     = np.flatnonzero(Dif < maxTime)
    xIndxDif = xIndxDif[indx]
    yIndxDif = yIndxDif[indx]
    Dif      = Dif[indx]

    return (Dif, xIndxDif, yIndxDif)

########################################################################################################################

def LoadBehFile(filename = None, InitialDir=''):

    if not filename:     
        if InitialDir and os.path.isdir(InitialDir):
            p = InitialDir
        else:
            p = ''
            
        filename = QtGui.QFileDialog.getOpenFileNameAndFilter(caption='Select a *.beh file',
                                                      filter='*.beh',
                                                      directory = p)
        filename = str(filename[0])
        if len(filename) == 0: return

    if not os.path.isfile(filename): return
    
    Data               = loadmat(filename)
    Data               = Data['Data']
    Data['Subject']    = str(Data['Subject'])
    Data['File']       = str(Data['File'])
    Data['MSN']        = str(Data['MSN'])
    Data['Box']        = int(Data['Box'])
    Data['Experiment'] = int(Data['Experiment'])

    Stims = FindStims(Data)
    if not Stims:
        Data = GetBhvParams(Data)
    else:
        for k in Stims:
            if Data[k].has_key('Descr'):
                Data[k]['Descr'] = str(Data[k]['Descr'])
    return Data

########################################################################################################################

def FindStims(Data):
    Stims=[k for k in Data.keys() if k.find('Stim')!=-1]
    Stims.sort()
    return Stims

########################################################################################################################

def GetRatNames(prefix = 'HMV', pth = ''):

    if not pth:
        pth = str(QtGui.QFileDialog.getExistingDirectory(caption = 'Beh Files Directory'))
        if not pth: return
    
    files = glob(os.path.join(pth, '*'))
    names = []
    for k in files:
        r=re.search('%s[0-9]{1,2}(?=_)' % prefix, os.path.split(k)[1])
        if r:
            names.append(r.group())

    ratnames = np.unique(names)
    ratnums = []
    for k in ratnames:
        ratnums.append(int(re.search('(?<=%s)[0-9]{1,2}' % prefix, k).group()))
    ratnums = np.array(ratnums)
    return list(ratnames[ratnums.argsort()])

########################################################################################################################

def SplitMPCFiles(filename, outpth):

    import string
    import shutil

    fid=open(filename,'rU')
    w=fid.readlines()
    fid.close()

    filehdr=w[0]
    lines='\n\n\n'
    st=[]

    for x,y in enumerate(w):
        if y.find('Start Date')!=-1:
            st.append(x)

    st.append(len(w))
    pth,outfile=os.path.split(filename)

    if os.path.isdir(os.path.join(outpth,outfile.replace(' ','_'))):
        shutil.rmtree(os.path.join(outpth,outfile.replace(' ','_')))

    os.mkdir(os.path.join(outpth,outfile.replace(' ','_')))

    for x,y in enumerate(st):
        if x+2>len(st):
            break
        fid=open(os.path.join(outpth,outfile.replace(' ','_'),outfile+string.lowercase[x]),'w')
        temp=w[y:st[x+1]]
        temp.insert(0,lines)
        temp.insert(0,filehdr)
        fid.writelines(temp)
        fid.close()

########################################################################################################################

def GetFilenames(RatName, RegExp='1T_REW[0-9].beh|NP0[0-9]A?.beh', BhvDir = ''):

    if not os.path.isdir(BhvDir):
        print 'That directory does not exist !'
        return

    f = glob(os.path.join(BhvDir,RatName+'_*.beh'))
    if not f: return
    f.sort()
    files=[]; MSN=[]

    for k in f:
        if re.search(RegExp, k):
            files.append(k)
            if re.search('(?<=_[0-9]{4}_)[HMV]?.*(?=\.beh)',k):
                MSN.append(re.search('(?<=_[0-9]{4}_)[HMV]?.*(?=\.beh)',k).group())
                if re.search('HMV_',MSN[-1]):
                    MSN[-1]=re.search('(?<=HMV_).*',MSN[-1]).group()

    return (files, MSN)

########################################################################################################################

def rmBehFiles(pth = '', pattern = '*.beh'):
    # Deletes only the files that match a certain pattern

    if not pth:
        pth = str(QtGui.QFileDialog.getExistingDirectory(caption = 'Dir to Delete'))
        if not pth: return

    pattern=os.path.join(pth, pattern)    
    files = glob(pattern)
    
    if files:
        for k in files:
            if os.path.isfile(k):
                os.remove(k)
    else:
        print 'nothing to delete !'

########################################################################################################################

def PrintDataSummary(Data):

    print '================================================================================'
    for k in ['Subject','StartDate']:
        print '%s:\t%s' % (k, str(Data[k]))
    for k in ['Box','MSN']:
        print '%s:\t\t%s' % (k, str(Data[k]))

    Stims=[k for k in Data.keys() if k.find('Stim')!=-1]
    if not Stims:
        return
    else:
        Stims.sort()

    keys=['HitsTS','MissTS','ErrTS']

    for k in Stims:
        print 'Stim: %s\t' % Data[k]['Descr'],
        for j in keys:
            if Data[k].has_key(j):
                print '%s:\t%s\t' % ( j[0:-2], str(Data[k][j].size) ),
        if Data[k].has_key('RTT') and Data[k]['RTT'].size>0:
            print 'RTT:\t%0.2f\t' % np.mean(Data[k]['RTT'])
        else:
            print
    print '================================================================================'

########################################################################################################################

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

########################################################################################################################
