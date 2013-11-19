#bhvfuncs.py

import scipy.io as spio
import numpy as np
import os, re, string
from glob import glob
from PyQt4 import QtGui
import ipdb

########################################################################################################################

def mpc2beh(filesdir = '', savedir = '', overwrite = False, calcParams = True):

    '''Routine to transform the MPC files into matlab files\n
    Input:
            *filesdir=can be a direcotyr or a filename
            *savedir= must be a directory
            *overwrite= whether to overwrite the existing file or not
    Output: data.'''

    if not filesdir:
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

    pd = QtGui.QProgressDialog('Running MPC to beh ...', 'Cancel', 0, len(files))
    pd.setWindowTitle('Converting Med PC files to to beh ...')
    pd.setGeometry(500, 500, 500, 100)
    pd.show()

    for k in files:

        YY  = k[3:5]
        MM  = k[6:8]
        DD  = k[9:11]
        HH  = k[12:14]
        MI  = k[15:17]
        RAT = k[k.find('Subject')+8:]

        filename = '%s_%s%s%s_%s%s*.beh' % (RAT,YY,MM,DD,HH,MI)

        # update the progress bar and also check if the conversion was canceled
        pd.setLabelText('Processing... ' + k)
        pd.setValue(pd.value()+1)
        QtGui.QApplication.processEvents()
        if pd.wasCanceled():
            return

        if glob(os.path.join(savedir,filename)) != [] \
           and overwrite == 0 \
           or not os.path.isfile(os.path.join(filesdir,k)):
            continue

        fid = open(os.path.join(filesdir, k), 'rU')
        w   = fid.readlines()
        fid.close()

        Data = {}

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
            
            elif e[0].find('\\') != -1:
                Data['Comments'] = e[0]
                continue
                
            else:
                e = e[1].strip()
                e = re.split(' +', e)
                e = [float(x) for x in e]
                if curvar == 'X' and e == [0,0,0,0,0]:
                    continue
                else:
                    Data[curvar] = Data[curvar] + e

        Data['X'] = np.array(Data['X'])
        if Data['X'].size < 10: continue
        Data['X'] = Data['X'][ np.nonzero(Data['X']) ]
        ECodes    = np.floor( Data['X']/1000000 )
        Data['X'] = Data['X'] - ECodes*1000000
        X         = np.unique(ECodes)
        X.sort()
        Data['X'] = Data['X'] - Data['X'][0]

        #create a np.object to be saved as a cell array in matlab
        Data['EventTS']   = []
        Data['EventCode'] = X
        Data['EventName'] = []

        TableEvent     = BehMapping(GetMapping(Data))
        TableEventName = TableEvent[:,0]
        TableEventCode = TableEvent[:,1]

        for l in X:
            indx = np.flatnonzero(ECodes == l)
            if np.any(indx):
                Data['EventTS'].append(Data['X'][ECodes == l])
                Data['EventName'].append(TableEventName[l == TableEventCode][0])

        Data['EventTS']   = np.array(Data['EventTS']  , dtype = np.object, ndmin=1)
        Data['EventCode'] = np.array(Data['EventCode'], dtype = np.object, ndmin=1)
        Data['EventName'] = np.array(Data['EventName'], dtype = np.object, ndmin=1)
        
        if calcParams:
            Data = GetBhvParams(Data)
            
        #eliminate the variables coming from MedPC
        Data.pop('X')

        savefile = os.path.join(savedir, filename.replace('*', '_' + Data['MSN']))

        if os.path.isfile(savefile):
            os.remove(savefile)

        spio.savemat(savefile, {'Data':Data}, format = '5',
                     appendmat=False, oned_as='row')

########################################################################################################################

def GetMapping(Data):

    if re.search('(CWHT|CENTWHT|[0-9]{1,2}K)[R,L]_[0-9]{1,2}K[R,L]|[R,L]Sip', string.upper(Data['MSN'])) or \
       re.search('UNBIAS|DRIFT|LEFT|DU2|LITI|REPCATCH|2T|NOISE|CLICK', string.upper(Data['MSN'])) or \
       np.any(Data['EventCode'] > 50):
        return 2
    else:
        return 1

########################################################################################################################

def BehMapping(Mapping):
    '''
    This function returns a table of names and codes, depending on the mapping
    Mapping is an integer that can be 1 or 2
    '''

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

    if Mapping == 1:
        return TableEvent1
    elif Mapping == 2:
        return TableEvent2

########################################################################################################################

def GetBhvParams(Data):
    '''Get behavioral parameters\n
    Input: behavioral Data dictionary or file\n
    Output: behavioral data dict with all behavioral parameters calculated'''

    try:
        if os.path.isfile(Data):
            Data = LoadBehFile(Data)
    except TypeError:
        if type(Data)!=dict:
            raise NameError('Data is neither file nor dictionary !')

    #===========================================================================================================#

    if GetMapping(Data) == 1:

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

                CurStim = 'Stim' + str(j)
                Data[CurStim] = {}
                Stim = Vars[k]
                # First find the hits and the misses
                # Then find the closest event to each one of the Hits
                HitsParams = GetHits(Stim, Vars['Lick'])
                RTT = HitsParams['ThirdLickHitTS'] - HitsParams['StimHitsTS']
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

    elif GetMapping(Data) == 2:

        # First get the table mapping
        TableEvent = BehMapping( GetMapping(Data) )
        
        # create a table to remap names and codes
        EvtNames   = np.array([['NPLed',     43],
                               ['RSipLight', 33],
                               ['LSipLight', 45],
                               ['WhtFLight', 37],
                               ['RedFLight', 35],
                               ['Noise',     47],
                               ['Tone1',     53],
                               ['Tone2',     54],
                               ['Catch',     55],
                               ['NpOut',     26],
                               ['NpIn',      25],
                               ['RRpIn',     23],
                               ['LRpIn',     29],
                               ['RLick',     21],
                               ['LLick',     27],
                               ['Solnd1',    49],
                               ['Solnd2',    50],
                               ['Solnd3',    51],
                               ['Solnd4',    52]], dtype = np.object)

        # get the variables present in the data structure into a dictionary
        # to make its manipulation easier
        Vars = {}
        for name, code in EvtNames:
            indx = np.flatnonzero(TableEvent[:,1] == code)
            if TableEvent[indx,0] in Data['EventName']:
                Vars[name] = Data['EventTS'][ Data['EventName'] == TableEvent[indx,0] ][0]

        #Check the association between stimuli and reward
        #If data has the 'S' field it builds an array of Stim - Resp association
        if Data.has_key('S'):
            
            # fill the S array and reshape it
            S = np.array(Data['S'])
            if S.size%10 == 0:
                S = S.reshape([S.size/10, 10])
            else:
                n = int(np.ceil(S.size/10.0)*10)
                S = np.concatenate([S, np.zeros([n-S.size])]).reshape([n/10,10])
                
            Data['S'] = S
            
            # iterate over the columns that contain stimuli -resp associations:
            # element 1 contains the association: 1-->right, 2-->left, 3-->miss (catch)
            # The idea is to create a table with the following columns:
            # stim-resp association; stim name, stim code, correct lick, incorrect lick
            StimResp = []
            for k in np.flatnonzero(S[0,:]):
                # stim -response association
                if S[:,k][0] == 1:
                    StimResp.append([1, EvtNames[EvtNames[:,1] == S[:,k][1], 0][0],
                                     S[:,k][1], 'RLick', 'LLick'])
                
                elif S[:,k][0] == 2:
                    StimResp.append([2, EvtNames[EvtNames[:,1] == S[:,k][1], 0][0],
                                     S[:,k][1], 'LLick', 'RLick'])
                
                elif S[:,k][0] == 3:
                    StimResp.append([3, EvtNames[EvtNames[:,1] == S[:,k][1], 0][0],
                                     S[:,k][1], 'Miss', 'Miss'])

        #If not, create one by default
        else:
            StimResp = [[1, 'Tone1', 53, 'RLick', 'LLick'],
                        [2, 'Tone2', 54, 'LLick', 'RLick'],
                        [3, 'Catch', 55, 'Miss' , 'Miss' ]]

        #main loop that extracts and calculates all the parameters for a given stimuli
        for j, k in enumerate(StimResp):

            # set current stimulus name
            CurStim = 'Stim' + str(j)

            # create an empty dictionary to hold the data
            Data[CurStim] = {}

            # check whether the StimTS and the lickTS are present and that this
            # is not a catch stimuli
            if k[1] in Vars and k[3] in Vars and k[0] != 3:
                
                # only add the stimuli that were presented after the first yellow LED in the nose poke
                if 'NPLed' in Vars:
                    StimTS = Vars[ k[1] ][ Vars[ k[1] ] > Vars['NPLed'][0] ]
                else:
                    StimTS = Vars[k[1]]
                
                ### ADD ONLY THE VALID STIMULI --> THOSE THAT HAVE AT LEAST
                ### 3 SECONDS FROM THE LAST TIMESTAMP
                if (Data['X'][-1] - StimTS[-1]) < 3.0:
                    StimTS = StimTS[0:-1]
                    
                hLick  = Vars[ k[3] ]

                # First find the hits
                HitsParams = GetHits(StimTS, hLick)

                # check whether are there any incorrects in the vars dictionary
                if k[4] in Vars:

                    # get the parameters for the incorrects
                    eLick  = Vars[ k[4] ]
                    IncParams  = GetHits(StimTS, eLick)

                    # for the case of no errors
                    if HitsParams['StimHitsTS'].size > 0 and IncParams['StimHitsIndx'].size == 0:
                        HitsTS   = HitsParams['StimHitsTS']
                        HitsIndx = HitsParams['StimHitsIndx']
                        MissTS   = HitsParams['StimMissTS']
                        MissIndx = HitsParams['StimMissIndx']
                        ErrTS    = np.array([])
                        ErrIndx  = np.array([])

                    # for the case of Hits and Errors
                    elif HitsParams['StimHitsTS'].size > 0 and IncParams['StimHitsIndx'].size > 0:

                        # Look for the Stim indices that are shared by Hits and Errors
                        indx = np.intersect1d(HitsParams['StimHitsIndx'],
                                              IncParams ['StimHitsIndx'])

                        # When there is intersection between hits and errors
                        # this can happen if the animal licks in one side and then the other
                        # obtain the true hits and the true errors
                        if indx.size > 0:

                            # Get the total reaction time for hits and errors
                            hRTT = HitsParams['ThirdLickHitTS'] - HitsParams['StimHitsTS']
                            eRTT = IncParams ['ThirdLickHitTS'] - IncParams ['StimHitsTS']

                            # check the reaction time for each intersection case and
                            # get the indices of the minimum
                            hIndx = np.searchsorted(HitsParams['StimHitsIndx'], indx)
                            eIndx = np.searchsorted(IncParams ['StimHitsIndx'], indx)
                            minIndx = np.flatnonzero(np.argmin([ hRTT[hIndx], eRTT[eIndx] ], axis = 0))

                            # eliminate those indices that are shared
                            HitsParams['StimHitsIndx'] = np.delete( HitsParams['StimHitsIndx'], indx[minIndx])
                            IncParams['StimHitsIndx']  = np.delete( IncParams['StimHitsIndx'],  indx[minIndx] == False)

                            # Get the stimulus timestamps again
                            HitsParams['StimHitsTS'] = StimTS[HitsParams['StimHitsIndx']]
                            IncParams['StimHitsTS']  = StimTS[IncParams['StimHitsIndx']]

                            # With the true hit indices recalculate all the hits and errors
                            HitsParams = GetHits(HitsParams['StimHitsTS'], hLick)
                            IncParams  = GetHits(IncParams ['StimHitsTS'], eLick)

                        # now get the misses
                        MissIndx = np.arange(StimTS.size)
                        tmp      = np.concatenate([HitsParams['StimHitsIndx'],
                                                   IncParams['StimHitsIndx']])
                        MissIndx = np.delete(MissIndx, tmp)
                        MissTS   = StimTS[MissIndx]

                        # ... and the rest of the parameters
                        HitsTS   = HitsParams['StimHitsTS']
                        HitsIndx = HitsParams['StimHitsIndx']
                        ErrTS    = IncParams['StimHitsTS']
                        ErrIndx  = IncParams['StimHitsIndx']

                    # for the case of no hits and errors > 0
                    elif HitsParams['StimHitsTS'].size == 0 and IncParams ['StimHitsIndx'].size > 0:

                        # get the incorrects
                        ErrTS   = IncParams['StimHitsTS']
                        ErrIndx = IncParams['StimHitsIndx']

                        # Calculate the true Misses
                        MissIndx = np.delete(range(len(StimTS)), IncParams['StimHitsIndx'])
                        MissTS = StimTS[MissIndx]

                        # the hits are simply empty arrays
                        HitsTS   = np.array([])
                        HitsIndx = np.array([])

                    # for the case of no Hits and no errors
                    elif HitsParams['StimHitsTS'].size == 0 and IncParams ['StimHitsIndx'].size == 0:
                        HitsTS   = np.array([])
                        HitsIndx = np.array([])
                        ErrTS    = np.array([])
                        ErrIndx  = np.array([])
                        MissTS   = StimTS
                        MissIndx = np.arange(StimTS.size)

                # in case there are no incorrects in the Vars dictionary
                else:
                    HitsTS   = HitsParams['StimHitsTS']
                    HitsIndx = HitsParams['StimHitsIndx']
                    MissTS   = HitsParams['StimMissTS']
                    MissIndx = HitsParams['StimMissIndx']
                    ErrTS    = np.array([])
                    ErrIndx  = np.array([])

                # fill the data structure with the information we have calculated
                Data[CurStim]['Descr']    = str(k[1])
                Data[CurStim]['HitsTS']   = HitsTS
                Data[CurStim]['HitsIndx'] = HitsIndx
                Data[CurStim]['ErrTS']    = ErrTS
                Data[CurStim]['ErrIndx']  = ErrIndx
                Data[CurStim]['MissTS']   = MissTS
                Data[CurStim]['MissIndx'] = MissIndx
                Data[CurStim]['StimTS']   = StimTS
                
                # add nose poke information
                if Vars.has_key('NpIn'):  Data[CurStim]['NpIn']  = Vars['NpIn']
                if Vars.has_key('NpOut'): Data[CurStim]['NpOut'] = Vars['NpOut']
                if   k[3]=='LLick': RpIn = Vars['LRpIn']
                elif k[3]=='RLick': RpIn = Vars['RRpIn']
                Data[CurStim]['RpIn'] = RpIn

                # add the lick information
                if k[3] in Vars.keys():
                    Lick = Vars[k[3]]
                    Data[CurStim]['Lick'] = Vars[k[3]]
                    
                # add the appropriate solenoid timestamps
                if k[0] == 1:
                    if 'Solnd1' in Vars and Vars['Solnd1'].size > 0:
                        Data[CurStim]['Solnd'] = Vars['Solnd1']
                    elif 'Solnd2' in Vars and Vars['Solnd2'].size > 0:
                        Data[CurStim]['Solnd'] = Vars['Solnd2']
                elif k[0] == 2:
                    if 'Solnd3' in Vars and Vars['Solnd3'].size > 0:
                        Data[CurStim]['Solnd'] = Vars['Solnd3']
                    elif 'Solnd4' in Vars and Vars['Solnd4'].size > 0:
                        Data[CurStim]['Solnd'] = Vars['Solnd4']

                # Calculate Reaction Times if are there any hits
                if HitsTS.size > 0:

                    # get the total reaction time and add it to the data structure
                    RTT = HitsParams['ThirdLickHitTS'] - HitsParams['StimHitsTS']
                    Data[CurStim]['RTT'] = RTT

                    # check the presence of the nosepoke variables as well as
                    # that the training protocol is a NosePoke task
                    if Vars.has_key('NpIn') and Vars.has_key('NpOut') and \
                       re.search('NP(?=0[0-9]A?)[0]', Data['MSN']):

                        # calculate the foreperiod
                        RT0, _, _ = SparseDistance(HitsTS, Vars['NpIn'], direction = 'yx')

                        #Calculate RT1
                        RT1, RT1x, RT1y = SparseDistance(HitsTS, Vars['NpOut'], direction = 'xy')
                        #pdb.set_trace()
                        #Calculate RT2 (from NP exit to Resp Port In)
                        RT2, RT2x, RT2y = SparseDistance(Vars['NpOut'][RT1y], RpIn, direction = 'xy')

                        #Calculate RT3 (from Resp Port In to First Lick)
                        RT3, RT3x, RT3y = SparseDistance(RpIn[RT2y], Lick, direction = 'xy')

                        # RT4 is the time from the first lick to the third lick
                        RT4 = HitsParams['ThirdLickHitTS'] - HitsParams['FirstLickHitTS']

                    else:
                        RT0 = np.array([])
                        RT1 = np.array([])
                        RT2 = np.array([])
                        RT3 = np.array([])
                        RT4 = np.array([])

                else:
                        RT0 = np.array([])
                        RT1 = np.array([])
                        RT2 = np.array([])
                        RT3 = np.array([])
                        RT4 = np.array([])
                        RTT = np.array([])

                Data[CurStim]['RT0'] = RT0
                Data[CurStim]['RT1'] = RT1
                Data[CurStim]['RT2'] = RT2
                Data[CurStim]['RT3'] = RT3
                Data[CurStim]['RT4'] = RT4
                Data[CurStim]['RTT'] = RTT

            # Calculate the parameters for the catch trials
            elif k[0] == 3 and k[1] in Vars:

                # only add the stimuli that were presented after the first
                # yellow LED in the nose poke
                if 'NPLed' in Vars:
                    StimTS = Vars[ k[1] ][ Vars[ k[1] ] > Vars['NPLed'][0] ]
                else:
                    StimTS = Vars[k[1]]
                                    
                # eliminate those that are not valid
                if (Data['X'][-1] - StimTS[-1]) < 3.0:
                    StimTS = StimTS[0:-1]

                if Vars.has_key('LLick') and Vars.has_key('RLick'):
                    Lick = np.concatenate((Vars['LLick'], Vars['RLick']))
                elif 'LLick' not in Vars.keys():
                    Lick = Vars['RLick']
                elif 'RLick' not in Vars.keys():
                    Lick = Vars['LLick']

                Lick.sort()

                # get the hits for the catch trials
                CatchParams = GetHits(StimTS, Lick)

                Data[CurStim]['Descr']    = 'Catch'
                Data[CurStim]['HitsTS']   = CatchParams['StimHitsTS']
                Data[CurStim]['HitsIndx'] = CatchParams['StimHitsIndx']
                Data[CurStim]['MissTS']   = CatchParams['StimMissTS']
                Data[CurStim]['MissIndx'] = CatchParams['StimMissIndx']
                Data[CurStim]['StimTS']   = StimTS
                Data[CurStim]['Lick']     = Lick

                if Vars.has_key('LRpIn') and Vars.has_key('RRpIn'):
                    RpIn = np.concatenate([Vars['LRpIn'], Vars['RRpIn']])
                elif Vars.has_key('LRpIn') and not Vars.has_key('RRpIn'):
                    RpIn = Vars['RRpIn']
                elif Vars.has_key('RRpIn') and not Vars.has_key('LRpIn'):
                    RpIn = Vars['LRpIn']

                RpIn.sort()

                Data[CurStim]['RpIn'] = RpIn

                if Vars.has_key('NpIn'):  Data[CurStim]['NpIn']  = Vars['NpIn']
                if Vars.has_key('NpOut'): Data[CurStim]['NpOut'] = Vars['NpOut']

                #pdb.set_trace()
                # now calculate the catch trial reaction times
                if CatchParams['StimHitsIndx'].size > 0:

                    # first get the total reaction time
                    RTT = CatchParams['ThirdLickHitTS'] - CatchParams['StimHitsTS']
                    Data[CurStim]['RTT'] = RTT

                    # get the foreperiods
                    if 'NpIn'  in Vars.keys():
                        RT0, RT0x, RT0y = SparseDistance(CatchParams['StimHitsTS'], Vars['NpIn'], direction  = 'yx')

                    # get RT1 (from stim to nose poke exit)
                    if 'NpOut' in Vars.keys():
                        RT1, RT1x, RT1y = SparseDistance(CatchParams['StimHitsTS'], Vars['NpOut'], direction = 'xy')

                        # then RT2 (from nose poke exito to response port in)
                        RT2, RT2x, RT2y = SparseDistance( Vars['NpOut'][RT1y], RpIn, direction = 'xy')

                        # then RT3 (from resp port in to first lick)
                        RT3, RT3x, RT3y = SparseDistance(RpIn[RT2y], Lick, direction = 'xy')

                         # RT4 is the time from the first lick to the third lick
                        RT4 = CatchParams['ThirdLickHitTS'] - CatchParams['FirstLickHitTS']

                else:
                    RT0 = np.array([])
                    RT1 = np.array([])
                    RT2 = np.array([])
                    RT3 = np.array([])
                    RT4 = np.array([])
                    RTT = np.array([])

                Data[CurStim]['RT0'] = RT0
                Data[CurStim]['RT1'] = RT1
                Data[CurStim]['RT2'] = RT2
                Data[CurStim]['RT3'] = RT3
                Data[CurStim]['RT4'] = RT4
                Data[CurStim]['RTT'] = RTT


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
    #DryLicks      = Dif[LickIndx[ValidStims], ValidStims]
    WetLicksIndx  = LickIndx[ValidStims] + WetLick -1
    WetLicks      = Dif[WetLicksIndx, ValidStims]

    Res['StimHitsIndx']     = np.flatnonzero( WetLicks <= RespWin )
    Res['StimMissIndx']     = np.delete( np.arange(StimTS.size), Res['StimHitsIndx'] )
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

    Inputs:
        x,y:
            vectors of timestamps of different length
        direction:
            "xy" if y happens after x "yx" x happens after y
        maxTime:
            maximum time lag between the two events

    Output:
        Dif:
            distances between the vectors.
        xIndxDif:
            indices of the first vector
        yIndxDif:
            indices of the second vector that give those differeces'''

    x   = np.array(x, ndmin = 2)
    y   = np.array(y, ndmin = 2)

    if x.size ==0 or y.size ==0:
        return np.array([]), np.array([]), np.array([])

    xx  = np.tile(x, (y.size, 1))
    yy  = np.tile(y, (x.size, 1)).transpose()

    if direction == 'xy':
        Dif = np.round(yy - xx, 3)
    elif direction == 'yx':
        Dif = np.round(xx - yy, 3)

    Dif[Dif < 0.00] = maxTime

    if x.size > y.size:
        xIndx = Dif.argmin(1)
        yIndx = Dif.argmin(0)[xIndx]
    else:
        yIndx = Dif.argmin(0)
        xIndx = Dif.argmin(1)[yIndx]

    Dif   = Dif[yIndx, xIndx]
    indx  = np.flatnonzero(Dif < maxTime)
    xIndx = xIndx[indx]
    yIndx = yIndx[indx]
    Dif   = Dif[indx]

    return (Dif, xIndx, yIndx)

########################################################################################################################

def SparseDistance2(x, y):

    xy     = np.concatenate([x,y],1)
    mainIndx   = np.argsort(xy)
    xySort = np.sort(xy)
    if x.size < y.size:
        Indx  = np.flatnonzero(mainIndx < x.size)
    else:
        Indx = np.flatnonzero(mainIndx > x.size)

    if (Indx[-1]+1) < xySort.size:
        dist = xySort[Indx+1]-xySort[Indx]
    else:
        pass
        # find the biggest Indx followed by

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

    else:
        filename = os.path.join(InitialDir, filename)

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
            for n in Data[k].keys():
                if Data[k].has_key('Descr'):
                    Data[k]['Descr'] = str(Data[k]['Descr'])
                Data[k][n] = np.array(Data[k][n])

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

def GetFilenames(RatName, RegExp = '1T_REW[0-9].beh|NP0[0-9]A?.beh', BhvDir = ''):

    if not os.path.isdir(BhvDir):
        print 'That directory does not exist !'
        return

    filesList = glob(os.path.join(BhvDir,RatName+'_*.beh'))
    if not filesList: return
    filesList.sort()
    files=[]; MSN=[]

    # iterate over the list of files
    for f in filesList:
        
        # if the regular expression is found 
        if re.search(RegExp, f, re.IGNORECASE):
            files.append(f)

            # Try to find the MSN in the file name
            match = re.search('(?<=_[0-9]{4}_)[HMV]?.*(?=\.beh)', f)
            if match:
                MSN.append(match.group())
                # try to eliminate the 'HMV_' at the beggining of the MSN name
                if re.search('HMV_', MSN[-1]):
                    MSN[-1] = re.search('(?<=HMV_).*', MSN[-1]).group()

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
                if type(Data[k][j]) in [np.ndarray, list]:
                    print '%s:\t%s\t' % ( j[0:-2], str(Data[k][j].size) ),
                else:
                    print '%s:\t%s\t' % ( j[0:-2], '' ),
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
