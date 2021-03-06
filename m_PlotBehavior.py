#PlotBehavior.py
# Module containing al the routines to plot behavioral events

########################################################################################################################
## Module Initialization Routine

import os, sys, re
from datetime import datetime as dtm
from m_NeuralPlot import rasterPlot

from scipy.stats import probplot, norm, gaussian_kde
from scipy.optimize import curve_fit

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

from PyQt4 import QtGui, QtCore

import guidata
import guidata.dataset.datatypes as dt
import guidata.dataset.dataitems as di
app = guidata.qapplication()

import m_bhvfuncs as bhv

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavToolbar
from matplotlib.figure import Figure

import ipdb

# MPL Widget Class to embed in Qt
class MplWidget(FigCanvas):
    def __init__(self, parent=None):
        self.fig = Figure()
        self.fig.set_facecolor('w')
        
        FigCanvas.__init__(self, self.fig)
        if parent: self.setParent(parent)

        FigCanvas.setSizePolicy(self,
                                QtGui.QSizePolicy.Expanding,
                                QtGui.QSizePolicy.Expanding)
        FigCanvas.updateGeometry(self)
        
########################################################################################################################

def PlotLATER(Data=[], ax=[], OutlierCalc = 'pctile', PlotMiss=False,
              ShowFit=True, ShowTitle=False, ShowEq=False):

    if not Data: return

    # create a function to pass as parameter for the line fit
    def LineFuncFit(x, m, n): return m*x+n

    # x values vector for the fit
    xVal = np.linspace(-50,10)

    rc('font', size=10)
    rc('font', family='monospace')
    rc('font', serif='Bitstream Vera Sans Mono')

    # check if ax is an axes instance; if not, create one 
    if ax == [] or not isinstance(ax, Axes):
        plt.figure(facecolor='w')
        ax=plt.subplot(111)

    # color list for different stimuli
    col=['b','g','r','c','m']

    # get the simuli in the data structure
    Stims = bhv.FindStims(Data)

    # iterate over stimuli
    for j,k in enumerate(Stims):

        # get the RTs
        RT1,_ = bhv.DistXY(Data[k]['StimTS'], Data[k]['NpOut'])

        # replace zeros with 0.001
        RT1[RT1==0] = 1e-3

        # do a probability plot
        a ,_ = probplot(-1/RT1, dist='norm')
        p,   = ax.plot(a[1], a[0], '+', ms=7, mew=2, color=col[j],
                       alpha=0.5, label=str(Data[k]['Descr'])+' Hits')

        # calculate a fit taking the central part of the distribution
        if ShowFit and len(RT1) > 10:
            
            # choose the method to calculate outliers
            if OutlierCalc == 'pctile':
                l, u   = np.percentile(a[1], (15,95))
            elif OutlierCalc == 'mc':
                m,l,u = Outliers4Skew(a[1])

            # get indices of values and plot a line fit    
            indx   = np.flatnonzero((a[1] >l) & (a[1]<u ))
            popt, _= curve_fit(LineFuncFit, a[1][indx], a[0][indx])
            ax.plot(xVal, LineFuncFit(xVal, popt[0], popt[1]),
                     lw=2, color=col[j], alpha=0.5, label='_nolegend_')

        if ShowEq:
            bp = dict(boxstyle="round", alpha=0.7, fc='w', ec=[.5,.5,.5])
            ax.text(-9,0,r'$y_{Hits}=%0.2fx + %0.2f$' % (m, n), fontsize = 20, bbox = bp)

    # set limits
    ax.set_xlim(-8,0)
    ax.set_ylim(-3.5,3.5)
    
    # get the percentiles of a normal distribution for the y axes ticks
    NDist  = norm(0, 1)
    Pctils = np.array([1,5,10,25,50,75,90,95,99])/100.0
    ytcks   = [NDist.ppf(k) for k in Pctils]
    ax.set_yticks(ytcks)
    ax.set_yticklabels([str(k) for k in Pctils])

    # set x axes ticks and tick labels
    xtcks  = np.array([0.1, 0.2, 0.3,0.5, 1, 5])
    ax.set_xticks(-1/xtcks)
    ax.set_xticklabels([str(k) for k in xtcks])

    # add legend and grid and axes labels
    ax.legend(loc='upper left', fancybox=True, prop={'size':10})
    ax.grid(alpha=0.6)
    ax.set_ylabel('Probability')
    ax.set_xlabel('RT1 (msec)')

    # if desired show title
    if ShowTitle:
        ax.set_title(Data['Subject']+'  '+'%02d/%02d/%02d' % tuple(Data['StartDate'])+'\n'+Data['MSN'],
                     fontsize=10)
        ax.figure.tight_layout()

########################################################################################################################

def PlotDensityFuncs(pth = ''):
    '''Examine the CDFs and/or PDFs of each RT'''

    if not pth:
        pth = str(QtGui.QFileDialog.getExistingDirectory(caption = 'Beh Files Directory'))
        if not pth: return
        
    ratNames = bhv.GetRatNames(prefix = 'HMV', pth = pth)
    if not ratNames: return
    
    class SetParams(dt.DataSet):
        RatName = di.ChoiceItem('Choose a Rat',tuple(ratNames))
        Task    = di.StringItem('Regular Expr','NP0[0-9]A?.beh')
        PDF     = di.BoolItem('Plot PDF ?',default=True).set_pos(col=0)
        CDF     = di.BoolItem('Plot CDF ?',default=True).set_pos(col=1)
        Norm    = di.BoolItem('Normalize ?',default=False).set_pos(col=0)

    Params=SetParams()
    if Params.edit()==1:
        
        RatName = ratNames[Params.RatName]
        Task    = Params.Task
        files,_ = bhv.GetFilenames(RatName, RegExp=Params.Task, BhvDir = pth)
        if not 'files' in locals(): return
        
        Stim0=[]; Stim1=[]
        for k in files:
            Data = bhv.LoadBehFile(k)
            if not [x for x in Data.keys() if x.find('Stim')!=-1]:
                Data = bhv.GetBhvParams(Data)
            if Data.has_key('Stim0'):
                Stim0.append(Data['Stim0'])
            if Data.has_key('Stim1') and\
             Data['Stim1'].has_key('RTT') and\
              Data['Stim1'].has_key('RT1') and\
               Data['Stim1'].has_key('RT2') and\
                Data['Stim1'].has_key('RT3'):
                Stim1.append(Data['Stim1'])
        
        rc('font',size=9, family='monospace', serif='Bitstream Vera Sans Mono')

        r=np.linspace(0,1,len(files))
        b=np.linspace(1,0,len(files))

        keys = ['RTT','RT1','RT2','RT3']
        nS   = len(files)

        if Params.CDF:
            plt.figure(dpi=100,facecolor='w')
            
            for m,n in enumerate(keys):
                m = m+1
                plt.subplot(2, 2, m)
                
                for j,k in enumerate(Stim0):
                    if k[n].size>1:
                        kernel = gaussian_kde(k[n])
                        x = np.linspace(k[n].min(), k[n].min(), 100)
                        y = np.cumsum(kernel(x))
                        y = y/y[-1]
                        #if Params.Norm==1: y=y/sum(y)
                        plt.plot(x,y,color=[r[j],0,b[j]])
                        
                plt.title(n+' n='+str(nS)+' sessions')
                plt.xlim(-0.2,3)

                if m in [1,3]:
                    if Params.Norm==1: plt.ylabel('Probabilty')
                    else: plt.ylabel('Count')
                if m in [3,4]:
                    plt.xlabel('Time (sec)')

                plt.grid()

            plt.title(RatName + ' ' + Params.Task)
            plt.tight_layout()

        if Params.PDF:
            plt.figure(RatName, dpi=100, facecolor='w')
            for m,n in enumerate(keys):
                m=m+1
                plt.subplot(2,2,m)
                for j,k in enumerate(Stim0):
                    if k[n].size>1:
                        kernel = gaussian_kde(k[n])
                        x = np.linspace(k[n].min(), k[n].min(), 100)
                        y = kernel(x)
                        
                        #if Params.Norm==1: y=y/sum(y)
                        plt.plot(x,y,color=[r[j],0,b[j]])
                plt.title(n+' n='+str(nS)+' sessions')
                plt.xlim(-0.2,3)
                if m in [1,3]:
                    if Params.Norm==1: plt.ylabel('Probabilty')
                    else: plt.ylabel('Count')
                if m in [3,4]:
                    plt.xlabel('Time (sec)')
                plt.grid()
            plt.title(RatName + ' ' + Params.Task)
            plt.tight_layout()

########################################################################################################################

def RxTimeEvolution(pth = ''):
    """Load the required files for a given rat and task and extract the stimuli"""
    
    import scipy.stats as st
    
    ratNames = bhv.GetRatNames(prefix = 'HMV', pth = pth)
    class SetParams(dt.DataSet):
        RatName   = di.ChoiceItem('Choose a Rat', tuple(ratNames))
        Task      = di.StringItem('Regular Expr','(NP0[0-9]A?.beh)')
        SaveFig   = di.BoolItem(text='', label='SaveFig ?')

    Params=SetParams()
    
    if Params.edit()==1:
        RatName  = ratNames[Params.RatName]
        Task     = Params.Task
        files, _ = bhv.GetFilenames(RatName, RegExp=Task, BhvDir = pth)
        FigHandles=[]
        Dates=[]; Stims={}
        for k in range(10):
            Stims['Stim'+str(k)]=[]
        
        for k in files:
            Data=bhv.LoadBehFile(k)
            Dates.append(Data['StartDate'])
            if not [x for x in Data.keys() if x.find('Stim')!=-1]:
                Data=bhv.GetBhvParams(Data)
            s=[k for k in Data.keys() if k.find('Stim')!=-1]
            s.sort()
            
            for x in s:
                Stims[x].append(Data[x])

        for k in Stims.keys():
            if not any(Stims[k]):
                Stims.pop(k)

        for x in Stims:
            Stim=Stims[x]
            ##  Calculate the average RxTimes across sessions and plot them in a bar graph

            mRTT=[]; mRT1=[]; mRT2=[]; mRT3=[];
            eRTT=np.zeros([len(Stim),2]); eRT1=np.zeros([len(Stim),2]);
            eRT2=np.zeros([len(Stim),2]); eRT3=np.zeros([len(Stim),2])

            for j,k in enumerate(Stim):

                if k.has_key('RTT'):
                    k['RTT']=k['RTT'].flatten()
                    mRTT=np.append(mRTT, st.nanmean(k['RTT']))
                    if k['RTT'].size>1:
                        eRTT[j,:]=bootci(k['RTT'],100)
                    else:
                        eRTT[j,:]=[0,0]
                else:
                    mRTT = np.append(mRTT,0)

                if k.has_key('RT1'):
                    k['RT1']=k['RT1'].flatten()
                    mRT1 = np.append(mRT1, st.nanmean(k['RT1']))
                    if k['RT1'].size>1:
                        eRT1[j,:]=bootci(k['RT1'],100)
                    else:
                        eRT1[j,:]=[0,0]
                else:
                    mRT1 = np.append(mRT1,0)

                if k.has_key('RT2'):
                    k['RT2']=k['RT2'].flatten()
                    mRT2 = np.append(mRT2, st.nanmean(k['RT2']))
                    if k['RT2'].size>1:
                        eRT2[j,:]=bootci(k['RT2'],100)
                    else:
                        eRT2[j,:]=[0,0]
                else:
                    mRT2 = np.append(mRT2,0)

                if k.has_key('RT3'):
                    k['RT3']=k['RT3'].flatten()
                    mRT3 = np.append(mRT3, st.nanmean(k['RT3']))
                    if k['RT3'].size>1:
                        eRT3[j,:]=bootci(k['RT3'],100)
                    else:
                        eRT3[j,:]=[0,0]
                else:
                    mRT3 = np.append(mRT3,0)

            mRTT=np.nan_to_num(mRTT)
            mRT1=np.nan_to_num(mRT1)
            mRT2=np.nan_to_num(mRT2)
            mRT3=np.nan_to_num(mRT3)
            eRT1=[mRT1-eRT1[:,0],eRT1[:,1]-mRT1]
            eRT2=[mRT2-eRT2[:,0],eRT2[:,1]-mRT2]
            eRT3=[mRT3-eRT3[:,0],eRT3[:,1]-mRT3]

            # Set some figure properties
            rc('font',size=9)
            rc('font',family='monospace')
            rc('font',serif='Bitstream Vera Sans Mono')

            f=plt.figure(dpi=100,facecolor='w')
            FigHandles.append(f)
            plt.subplot2grid((1,5),(0,0),colspan=4)

            # Stacked barplots of the data
            b1=plt.bar(range(1,len(Stim)+1),mRT1,bottom=0,
                  yerr=eRT1,color=[.3,.5,.8],align='center',edgecolor='',
                  error_kw={'elinewidth':2.5})

            b2=plt.bar(range(1,len(Stim)+1),mRT2,bottom=mRT1,
                  yerr=eRT2,
                  color=[.3,.8,.5],align='center',edgecolor='',
                  error_kw={'elinewidth':2.5})

            x=mRT1+mRT2
            b3=plt.bar(range(1,len(Stim)+1),mRT3,bottom=x,
                  yerr=eRT3,color=[.8,.5,.3],align='center',edgecolor='',
                  error_kw={'elinewidth':2.5})

            plt.legend([b1,b2,b3],['RT1','RT2','RT3'],mode='none',ncol=1,fancybox=True,loc=0)
            plt.title(RatName+' '+Task)
            plt.ylim(0,1.1*mRTT.max())
            plt.xticks(range(1,len(Stim)+1),range(1,len(Stim)+1))
            plt.xlim(0,len(Stim)+1)
            plt.xlabel('Session Number')
            plt.ylabel('Time (sec)')
            plt.grid()

            # To plot the average Reaction times

            plt.subplot2grid((1,5),(0,4),colspan=4)

            mmRT1=st.nanmean(mRT1);
            eeRT1=bootci(mRT1,100)
            eeRT1=np.array([[mmRT1-eeRT1[0]],[eeRT1[1]-mmRT1]])
            b1=plt.bar(1, mmRT1, bottom=0, color=[.3,.5,.8], edgecolor='',
                  align='center', yerr=eeRT1, error_kw={'elinewidth':2.5})

            mmRT2=st.nanmean(mRT2);
            eeRT2=bootci(mRT2,100)
            eeRT2=np.array([[mmRT2-eeRT2[0],eeRT2[1]-mmRT2]])
            b2=plt.bar(1, mmRT2, bottom=mmRT1, color=[.3,.8,.5], edgecolor='',
                  align='center', yerr=eeRT2, error_kw={'elinewidth':2.5})

            b=mmRT1+mmRT2
            mmRT3=st.nanmean(mRT3);
            eeRT3=bootci(mRT3,100)
            eeRT3=np.array([[mmRT3-eeRT3[0],eeRT3[1]-mmRT3]])
            b3=plt.bar(1, mmRT3, bottom=b, color=[.8,.5,.3], edgecolor='',
                  align='center', yerr=eeRT3, error_kw={'elinewidth':2.5})

            plt.title('Average')
            plt.legend([b1,b2,b3],['RT1','RT2','RT3'],fancybox=True)
            plt.xlim(0,2)
            plt.ylim(0,1.1*max(mRTT))
            plt.xticks([1],['Avg'])
            plt.grid()
            plt.tight_layout()

        if Params.SaveFig:
            SaveFigure(FigHandles,FigName=RatName+'_'+Params.Task+'_RxTimesEvolution')
        

########################################################################################################################

def Raster2(Stim, Event, TWin=[1,4], axes=[], yAxVar = 'Trial', color='k', mew=1, alpha=1):

    if not axes:
        axes=plt.subplot(111)
    elif not isinstance(axes, Axes):
        axes=plt.subplot(111)
    
    Event2 = []
    y = []
    if yAxVar == 'Trial':
        for j,k in enumerate(Stim):
            tmp = Event[( Event > k - TWin[0] ) & ( Event < k + TWin[1] )] -k
            Event2.extend(tmp)
            y.extend(j*np.ones_like(tmp))
    elif yAxVar == 'Time':
        for k in Stim:
            tmp = Event[( Event > k - TWin[0] ) & ( Event < k + TWin[1] )] -k
            Event2.extend(tmp)
            y.extend(k*np.ones_like(tmp))
    
    axes.plot(Event2, y,'|', color=color, mew=mew, alpha=alpha, rasterized = True)
    axes.set_xlim(-TWin[0], TWin[1])
    
########################################################################################################################

def Raster3(Stim, EventTS, TWin=[1.0,4.0], ax=[], yAxVar = 'Trial', color='k', lw=1, alpha=1):
    '''Makes use of the new eventplot function in matplotlib 1.3'''

    if not ax:
        ax=plt.subplot(111)
    elif not isinstance(ax, Axes):
        ax=plt.subplot(111)
    #pdb.set_trace()
    Event2 = []
    y = []
    if yAxVar == 'Trial':
        for j,k in enumerate(Stim):
            tmp = EventTS[( EventTS > k - TWin[0] ) & ( EventTS < k + TWin[1] )] -k
            if tmp.size > 0:
                Event2.append(tmp)
            else:
                Event2.append(np.array(100, ndmin = 1))
                #Event2.append(None)
            y.append(j)
    elif yAxVar == 'Time':
        for k in Stim:
            tmp = EventTS[( EventTS > k - TWin[0] ) & ( EventTS < k + TWin[1] )] -k
            Event2.append(tmp)
            y.append(k)
    
    rasterPlot(Event2, ax=ax, color = color, lw = lw, alpha = alpha)
    #ax.eventplot(Event2, lineoffsets = y, colors=[color], alpha = alpha, lw = lw)
    ax.set_xlim(-TWin[0], TWin[1])

########################################################################################################################

class Settings(dt.DataSet):
    '''Select the paramters to save the currently active figure'''
    BhvDir  = di.DirectoryItem(label='Beh Dir')
    ImgDir  = di.DirectoryItem(label='Img Dir')
    Format  = di.ChoiceItem(label='Format', choices=[('.jpg','.jpg'),('.png','.png'),('.svg','.svg'),('.pdf','.pdf')])
    dpi     = di.IntItem(label='dpi', default=300, min=50, max=600, nonzero=True, slider=True)
    
settings = Settings()

class BhvRasters(QtGui.QMainWindow):
      
    def __init__(self, workingDir = None):
        if not workingDir:
            self.BhvDir = str(QtGui.QFileDialog.getExistingDirectory())
            if not self.BhvDir: return
        else:
            self.BhvDir = workingDir
            
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle('Raw Behavioral Data Browser')
        self.main_widget = QtGui.QWidget(self)
        self.main_layout = QtGui.QHBoxLayout(self.main_widget)
        self.main_layout.setMargin(0)
        self.main_layout.setSpacing(0)
        self.oldMSN = ''
        self.StimColors = ['b','g','r','c','m','y']
        
        ########### STATUS BAR ##############################################
        
##        self.StBar=self.statusBar()
##        self.StBar.showMessage('')
        
        ###################################################################################################

        self.left_panel = QtGui.QWidget(self)
        self.left_panel.setMaximumWidth(400)
        self.left_panel.setMinimumWidth(400)
        left_panel_lay = QtGui.QVBoxLayout(self.left_panel)
        left_panel_lay.setMargin(2)
        
        grp1 = QtGui.QGroupBox('Get Files', self.left_panel)
        vlay = QtGui.QVBoxLayout()
        vlay.setSpacing(2)
        vlay.setMargin(2)
        hlay = QtGui.QHBoxLayout()
        self.RegExpEdit = QtGui.QLineEdit(self.main_widget)
        self.RegExpEdit.setText('NP04A_(6K|CentWht)R_12KL')
        hlay.addWidget(QtGui.QLabel('RegExp',self))
        hlay.addWidget(self.RegExpEdit)
        vlay.addLayout(hlay)
        
        hlay = QtGui.QHBoxLayout()
        self.selRatName= QtGui.QComboBox(self.main_widget)
        self.selRatName.addItems(bhv.GetRatNames(pth = self.BhvDir))
        self.selRatName.currentIndexChanged.connect(self.GetFiles)
        hlay.addWidget(self.selRatName)
        
        self.GetFilesBtn = QtGui.QPushButton('GetFiles',self.main_widget)
        self.GetFilesBtn.clicked.connect(self.GetFiles)
        hlay.addWidget(self.GetFilesBtn)
        vlay.addLayout(hlay)

        self.SettingsBtn = QtGui.QPushButton('Settings',self.main_widget)
        self.SettingsBtn.clicked.connect(self.editSettings)
        vlay.addWidget(self.SettingsBtn)
        
        grp1.setLayout(vlay)
        left_panel_lay.addWidget(grp1)
##        left_panel_lay.addStretch(1)

        ###################################################################################################
        ## Parameters Group BOx
        
        grp2 = QtGui.QGroupBox('Parameters', self.left_panel)
        
        vlay = QtGui.QVBoxLayout()
        vlay.setSpacing(2)
        vlay.setMargin(2)
        
        hlay = QtGui.QHBoxLayout()
        self.AxNRows = QtGui.QSpinBox()
        self.AxNRows.setMinimum(1)
        self.AxNRows.setMaximum(3)
        self.AxNRows.setValue(2)
        hlay.addWidget(QtGui.QLabel('AxNRows'))
        hlay.addWidget(self.AxNRows)
        hlay.addStretch(1)
        self.AxNCols = QtGui.QSpinBox()
        self.AxNCols.setMinimum(1)
        self.AxNCols.setMaximum(5)
        self.AxNCols.setValue(3)
        hlay.addWidget(QtGui.QLabel('AxNCols'))
        hlay.addWidget(self.AxNCols)
        self.SetNAxesBtn = QtGui.QPushButton('Set Axes')
        self.SetNAxesBtn.clicked.connect(self.SetNAxesProc)
        hlay.addStretch(1)
        hlay.addWidget(self.SetNAxesBtn)
        vlay.addLayout(hlay)

        self.naxes = self.AxNRows.value()*self.AxNCols.value()
        
        hlay = QtGui.QHBoxLayout()
        self.YAxesVarCombo = QtGui.QComboBox()
        self.YAxesVarCombo.addItems(['Trial','Time'])
        hlay.addWidget(QtGui.QLabel('Y Axes Var'))
        hlay.addWidget(self.YAxesVarCombo)
        hlay.addStretch(1)
        
        self.TWin1 = QtGui.QDoubleSpinBox(self.main_widget)
        self.TWin1.setRange(0.0, 4.0)
        self.TWin1.setSingleStep(0.1)
        self.TWin1.setValue(1.0)
        self.TWin1.setFixedWidth(60)
        hlay.addWidget(QtGui.QLabel('TWIN1'))
        hlay.addWidget(self.TWin1)
        hlay.addStretch(1)
        
        self.TWin2 = QtGui.QDoubleSpinBox(self.main_widget)
        self.TWin2.setRange(0.1, 10.0)
        self.TWin2.setSingleStep(0.1)
        self.TWin2.setValue(4.0)
        self.TWin2.setFixedWidth(60)
        hlay.addWidget(QtGui.QLabel('TWIN2'))
        hlay.addWidget(self.TWin2)
        vlay.addLayout(hlay)

        hlay = QtGui.QHBoxLayout()
        self.KernelType = QtGui.QComboBox()
        self.KernelType.addItems(['Hamming', 'Square'])
        hlay.addWidget(QtGui.QLabel('Kernel Type'))
        hlay.addWidget(self.KernelType)
        self.KernelSize = QtGui.QSpinBox()
        self.KernelSize.setRange(60, 1000)
        self.KernelSize.setValue(240)
        hlay.addStretch(1)
        hlay.addWidget(QtGui.QLabel('Kern Size'))
        hlay.addWidget(self.KernelSize)
        vlay.addLayout(hlay)

        hlay = QtGui.QHBoxLayout()
        self.RTPlotTypeCombo = QtGui.QComboBox()
        self.RTPlotTypeCombo.addItems(['Hist StepF','Hist Step','Hist Cum','PDFs','CDFs'])
        self.RTPlotTypeCombo.currentIndexChanged.connect(self.HitsTypeChanged_Proc)
        hlay.addWidget(QtGui.QLabel('RTs Plot Type'))
        hlay.addWidget(self.RTPlotTypeCombo)
        hlay.addStretch(1)
        self.HistTWin = QtGui.QDoubleSpinBox()
        self.HistTWin.setRange(0.1, 10.0)
        self.HistTWin.setValue(2.0)
        self.HistTWin.setSingleStep(0.1)
        hlay.addWidget(QtGui.QLabel('HistTWin'))
        hlay.addWidget(self.HistTWin)
        hlay.addStretch(1)
        self.HistNBins = QtGui.QSpinBox()
        self.HistNBins.setRange(0, 100)
        self.HistNBins.setValue(20)
        hlay.addWidget(QtGui.QLabel('NBins'))
        hlay.addWidget(self.HistNBins)
        vlay.addLayout(hlay)
        
        self.PlotTypeTable = QtGui.QTableWidget(self.naxes,5,self)
##        self.PlotTypeTable.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
##        self.PlotTypeTable.updateGeometry()
##        self.PlotTypeTable.setMinimumHeight(350)
        self.PlotTypeTable.setVerticalHeaderLabels(['Ax%d' % k for k in range(1,self.naxes+1)])
        self.PlotTypeTable.setColumnWidth(0,65)
        self.PlotTypeTable.setColumnWidth(1,60)
        self.PlotTypeTable.setColumnWidth(2,60)
        self.PlotTypeTable.setColumnWidth(3,50)
        self.PlotTypeTable.setColumnWidth(4,90)
        self.PlotTypeTable.setHorizontalHeaderLabels(['PlotType','Stim','Ref','OutCm','Sort By'])
        
        self.PlotTypeCombo  = []
        self.Stim2PlotCombo = []
        self.ReferenceCombo = []
        self.SortByCombo    = []
        self.OutcomeCombo   = []
        
        
        for k in range(self.naxes):
            self.PlotTypeTable.setRowHeight(k,20)
            self.PlotTypeCombo.append(QtGui.QComboBox(self))
            self.PlotTypeCombo[-1].addItems(['Later', 'RT Dist', 'Perf Evol', 'Behavior'])
            self.PlotTypeCombo[-1].currentIndexChanged.connect(self.PlotTypeStatus)
            self.PlotTypeTable.setCellWidget(k,0,self.PlotTypeCombo[-1])

            self.Stim2PlotCombo.append(QtGui.QComboBox(self))
            for n in ['Stim0', 'Stim1', 'Stim2']:
                self.Stim2PlotCombo[-1].addItem(n, QtCore.QVariant(n))
            self.Stim2PlotCombo[-1].currentIndexChanged.connect(self.Stim2PlotChanged_Proc)
            self.PlotTypeTable.setCellWidget(k,1,self.Stim2PlotCombo[-1])

            self.ReferenceCombo.append(QtGui.QComboBox(self))
            self.ReferenceCombo[-1].addItems(['CentNP','Stim','NpExit','RpIn','1stLick'])
            self.ReferenceCombo[-1].setCurrentIndex(1)
            self.PlotTypeTable.setCellWidget(k,2,self.ReferenceCombo[-1])
            
            self.OutcomeCombo.append(QtGui.QComboBox(self))
            self.OutcomeCombo[-1].addItems(['All','Hits','Error','Miss'])
            self.PlotTypeTable.setCellWidget(k,3,self.OutcomeCombo[-1])
            
            self.SortByCombo.append(QtGui.QComboBox(self))
            self.SortByCombo[-1].addItems(['No sort','CentNP','Stim','NpExit','RpIn','1stLick','3rdLick',
                                           'RT0/CentNP','RT0/Stim','RT0/NpExit','RT0/RpIn','RT0/1stLick','RT0/3rdLick'])
            self.PlotTypeTable.setCellWidget(k,4,self.SortByCombo[-1])

        self.PlotTypeCombo[0].setCurrentIndex(1)
        self.PlotTypeCombo[0].setCurrentIndex(0)
        self.PlotTypeCombo[3].setCurrentIndex(1)
        for k in [1,2,4,5]:
            self.PlotTypeCombo[k].setCurrentIndex(3)
        
        self.PlotTypeTable.setFont(QtGui.QFont(self.PlotTypeTable.font().family(),8))
        vlay.addWidget(self.PlotTypeTable)

        grp2.setLayout(vlay)
        left_panel_lay.addWidget(grp2)
        grp2.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        grp2.updateGeometry()

##        left_panel_lay.addStretch(1)
        ###################################################################################################
        ## Rat Meta info Panel

        grp4 = QtGui.QGroupBox('MetaInfo', self.left_panel)
        vlay = QtGui.QVBoxLayout()
        vlay.setSpacing(2)
        vlay.setMargin(2)
        
        self.MetaInfoTable = QtGui.QTableWidget(6,1,self)
        self.MetaInfoTable.setVerticalHeaderLabels(['RatName', 'File', 'Date', 'MSN', 'Box', 'Group'])
        self.MetaInfoTable.setHorizontalHeaderLabels([''])
        self.MetaInfoTable.horizontalHeader().setVisible(False)
        self.MetaInfoTable.setFont(QtGui.QFont(self.MetaInfoTable.font().family(),8))
        self.MetaInfoTable.setColumnWidth(0,350)
        for k in range(6): self.MetaInfoTable.setRowHeight(k,20)
        vlay.addWidget(self.MetaInfoTable)
        
        grp4.setLayout(vlay)
        left_panel_lay.addWidget(grp4)
        
        ###################################################################################################

        grp3 = QtGui.QGroupBox('Plot Control', self.left_panel)
        vlay = QtGui.QVBoxLayout()
        vlay.setSpacing(2)
        vlay.setMargin(2)
        
        self.PlotBtn = QtGui.QPushButton('Plot Figure',self.main_widget)
        self.PlotBtn.clicked.connect(self.Plot)
        vlay.addWidget(self.PlotBtn)
        
        hlay = QtGui.QHBoxLayout()
        self.PrevBtn = QtGui.QPushButton('<< Prev',self.main_widget)
        self.PrevBtn.clicked.connect(self.PrevPlot)
        hlay.addWidget(self.PrevBtn)
        self.NextBtn = QtGui.QPushButton('Next >>',self.main_widget)
        self.NextBtn.clicked.connect(self.NextPlot)
        hlay.addWidget(self.NextBtn)
        vlay.addLayout(hlay)

        hlay = QtGui.QHBoxLayout()
        self.SaveFigCheck = QtGui.QCheckBox('Save each fig ?', self.main_widget)
        hlay.addWidget(self.SaveFigCheck)
        self.SaveFigBtn = QtGui.QPushButton('Save Fig',self.main_widget)
        self.SaveFigBtn.clicked.connect(self.SaveFig)
        hlay.addWidget(self.SaveFigBtn)
        vlay.addLayout(hlay)

        self.PanelBtn = QtGui.QPushButton('Panel')
        self.PanelBtn.clicked.connect(self.showPanel_Proc)
        vlay.addWidget(self.PanelBtn)
        grp3.setLayout(vlay)
        left_panel_lay.addWidget(grp3)
        left_panel_lay.addStretch(1)
        self.main_layout.addWidget(self.left_panel)

        
        
        ###################################################################################################

        self.right_panel = QtGui.QWidget(self)
        
        self.main_fig = MplWidget(self.right_panel)
        self.fig = self.main_fig.figure
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.ntb = NavToolbar(self.main_fig, self.main_widget)
        self.ntb.setIconSize(QtCore.QSize(12,12))
        self.ntb.setFloatable(True)
        self.ntb.setMovable(True)

        vlay = QtGui.QVBoxLayout(self.right_panel)
        vlay.setMargin(0)
        vlay.setSpacing(0)
        self.selFiles = QtGui.QComboBox(self.main_widget)
##        QObj.connect(self.selFiles, QtCore.SIGNAL('currentIndexChanged(int)'), self.Plot)
        vlay.addWidget(self.selFiles)
        vlay.addWidget(self.main_fig)
        vlay.addWidget(self.ntb)

        self.right_panel.setLayout(vlay)
        self.main_layout.addWidget(self.right_panel)
##        self.main_fig.figure.canvas.mpl_connect('draw_event', self.draw_callback)

        self.ax = []
        for k in range(1,self.naxes+1):
            self.ax.append(self.main_fig.figure.add_subplot(self.AxNRows.value(),
                                                            self.AxNCols.value(), k))

        if sys.platform == 'linux2':
            QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Plastique'))
            
        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)
        self.show()

        self.PanelWidget = QtGui.QTableWidget(10,2)
        for k in range(self.PanelWidget.rowCount()):
            self.PanelWidget.setRowHeight(k,20)
            self.PanelWidget.setCellWidget(k,0,QtGui.QLabel('Hola'))

    def showPanel_Proc(self):
        self.PanelWidget.show()

    def editSettings(self):
        settings.edit()
        
    def HitsTypeChanged_Proc(self):
        if self.RTPlotTypeCombo.currentText() not in ['Hist StepF','Hist Step', 'Hist Cum']:
            self.HistNBins.setDisabled(True)
        else:
            self.HistNBins.setEnabled(True)
        
    def PlotTypeStatus(self, indx1, indx2=None):

        if indx2 is not None:
            indx = indx2
        else:
            sender = self.sender()
            indx = self.PlotTypeCombo.index(sender)

        if self.PlotTypeCombo[indx].currentText() == 'Behavior':
            if hasattr(self,'Data'):
                if re.search('RT', str(self.ReferenceCombo[indx].currentText()))\
                   or self.oldMSN != self.Data['MSN'] or str(self.Stim2PlotCombo[indx].itemText(0)) == 'All':
                    self.ReferenceCombo[indx].clear()
                    self.ReferenceCombo[indx].addItems(['CentNP','Stim','NpExit','RpIn','1stLick'])
                    self.ReferenceCombo[indx].setCurrentIndex(1)
                    self.Stim2PlotCombo[indx].clear()
                    for l in bhv.FindStims(self.Data):
                        self.Stim2PlotCombo[indx].addItem(self.Data[l]['Descr'], QtCore.QVariant(l))
            self.Stim2PlotCombo[indx].setEnabled(True)
            self.SortByCombo[indx].setEnabled(True)
            self.OutcomeCombo[indx].setEnabled(True)
            self.ReferenceCombo[indx].setEnabled(True)
            
        elif self.PlotTypeCombo[indx].currentText() == 'RT Dist':
            if hasattr(self,'Data'):
                if not re.search('RT', str(self.ReferenceCombo[indx].currentText())) or self.oldMSN != self.Data['MSN']:
                    self.ReferenceCombo[indx].clear()
                    self.ReferenceCombo[indx].addItems(['RT1','RT2','RT3','RT4','RTT'])
                    self.Stim2PlotCombo[indx].clear()
                    self.Stim2PlotCombo[indx].addItem('All', QtCore.QVariant('All'))
                    for k in bhv.FindStims(self.Data):
                        self.Stim2PlotCombo[indx].addItem(self.Data[k]['Descr'], QtCore.QVariant(k))
            self.Stim2PlotCombo[indx].setEnabled(True)
            self.SortByCombo[indx].setEnabled(False)
            self.OutcomeCombo[indx].setEnabled(False)
            self.ReferenceCombo[indx].setEnabled(True)
            
        elif self.PlotTypeCombo[indx].currentText() == 'Perf Evol':
            if str(self.Stim2PlotCombo[indx].itemText(0)) != 'All':
                self.Stim2PlotCombo[indx].clear()
                self.Stim2PlotCombo[indx].addItem('All', QtCore.QVariant('All'))
                for l in bhv.FindStims(self.Data):
                    self.Stim2PlotCombo[indx].addItem(self.Data[l]['Descr'], QtCore.QVariant(l))
                    
            self.Stim2PlotCombo[indx].setEnabled(True)
            self.SortByCombo[indx].setEnabled(False)
            self.OutcomeCombo[indx].setEnabled(False)
            self.ReferenceCombo[indx].setEnabled(False)

        elif self.PlotTypeCombo[indx].currentText() == 'Later':
            self.Stim2PlotCombo[indx].setEnabled(False)
            self.SortByCombo[indx].setEnabled(False)
            self.ReferenceCombo[indx].setEnabled(False)
            self.OutcomeCombo[indx].setEnabled(False)

    def Stim2PlotChanged_Proc(self):
        sender = self.sender()
        indx = self.Stim2PlotCombo.index(sender)
        if re.search('RT1',sender.currentText()):
            self.OutcomeCombo[indx].setCurrentIndex(1)
        else:
            self.OutcomeCombo[indx].setCurrentIndex(0)
            
    def SetNAxesProc(self):
        self.main_fig.figure.clear()

        self.naxes = self.AxNRows.value()*self.AxNCols.value()
        
        if self.naxes > len(self.ax):
            for k in range(len(self.ax), self.naxes):
                self.PlotTypeTable.insertRow(k)
                self.PlotTypeTable.setRowHeight(k,20)
                self.PlotTypeCombo.append(QtGui.QComboBox(self))
                self.PlotTypeCombo[-1].addItems(['Later', 'RT Dist', 'Perf Evol', 'Behavior'])
                self.PlotTypeCombo[-1].currentIndexChanged.connect(self.PlotTypeStatus)
                self.PlotTypeTable.setCellWidget(k,0,self.PlotTypeCombo[-1])

                self.Stim2PlotCombo.append(QtGui.QComboBox(self))
                if hasattr(self, 'Data'):
                    for l in bhv.FindStims(self.Data):
                        self.Stim2PlotCombo[-1].addItem(self.Data[l]['Descr'], QtCore.QVariant(l))
                self.PlotTypeTable.setCellWidget(k, 1, self.Stim2PlotCombo[-1])

                self.ReferenceCombo.append(QtGui.QComboBox(self))
                self.ReferenceCombo[-1].addItems(['CentNP','Stim','NpExit','RpIn','1stLick'])
                self.ReferenceCombo[-1].setCurrentIndex(1)
                self.PlotTypeTable.setCellWidget(k,2,self.ReferenceCombo[-1])
                
                self.OutcomeCombo.append(QtGui.QComboBox(self))
                self.OutcomeCombo[-1].addItems(['All','Hits','Error','Miss'])
                self.PlotTypeTable.setCellWidget(k,3,self.OutcomeCombo[-1])

                self.SortByCombo.append(QtGui.QComboBox(self))
                self.SortByCombo[-1].addItems(['No sort','CentNP','Stim','NpExit','RpIn','1stLick','3rdLick',
                                               'RT0/CentNP','RT0/Stim','RT0/NpExit','RT0/RpIn','RT0/1stLick','RT0/3rdLick'])
                self.PlotTypeTable.setCellWidget(k,4,self.SortByCombo[-1])
                
        elif self.naxes < len(self.ax):
            for k in range(len(self.ax)-1, self.naxes-1, -1):
                self.PlotTypeTable.removeRow(k)
                self.PlotTypeCombo.pop(k)
                self.Stim2PlotCombo.pop(k)
                self.SortByCombo.pop(k)
                self.OutcomeCombo.pop(k)
                self.ReferenceCombo.pop(k)
  
        self.PlotTypeTable.setVerticalHeaderLabels(['Ax%d' % k for k in range(1,self.naxes+1)])
        self.ax = []
        for k in range(1,self.naxes+1):
            self.ax.append(self.main_fig.figure.add_subplot(self.AxNRows.value(), self.AxNCols.value(),k))
            
        self.Plot()

    def GetFiles(self):
        self.selFiles.clear()        
        self.regExp  = str(self.RegExpEdit.text())
        self.ratname = str(self.selRatName.currentText())
        files,_= bhv.GetFilenames(self.ratname, self.regExp, self.BhvDir)
        if files: self.selFiles.addItems(files)

    def LoadData(self):
        self.Data = bhv.LoadBehFile(str(self.selFiles.currentText()))
        self.Map = bhv.GetMapping(self.Data)

    def SaveFig(self):
        FigName = os.path.split(str(self.selFiles.currentText()))[1]
        FigName = FigName.replace('.beh',settings.Format)
        self.fig.savefig(os.path.join(settings.ImgDir,FigName),
                         dpi = settings.dpi)
    
    def on_key(self, event):        
        if event.key=='alt': self.NextPlot()
        elif event.key=='control': self.PrevPlot()
        
    def PrevPlot(self):
        if self.selFiles.count() > 0 and self.selFiles.currentIndex() > 0:
            self.selFiles.setCurrentIndex(self.selFiles.currentIndex()-1)
##            self.StBarOldMsg = self.StBar.currentMessage()
##            self.StBar.clearMessage()
            self.Plot()
        else:
            pass
##            self.StBar.showMessage('You have reached the first file !')

    def NextPlot(self):
        c = self.selFiles.count()
        if c > 0 and self.selFiles.currentIndex() < c-1:
            self.selFiles.setCurrentIndex(self.selFiles.currentIndex()+1)
##            self.StBarOldMsg = self.StBar.currentMessage()
##            self.StBar.clearMessage()
            self.Plot()
        else:
            pass
##            self.StBar.showMessage('You have reached the last file !')
        
    def Plot(self):
        
        if self.selFiles.count()>0:
            Data  = bhv.LoadBehFile(str(self.selFiles.currentText()))
            Map   = bhv.GetMapping(Data)
            stims = bhv.FindStims(Data)
            self.Data = Data
            for j,k in enumerate(self.PlotTypeCombo):
                self.PlotTypeStatus(None, j)
            self.oldMSN = Data['MSN']
        else:
            return
        
        if Map==1:
            Resp=['Lick','RpIn']
        elif Map==2:
            Resp=['Lick','RpIn']

        self.MetaInfoTable.setItem(0, 0, QtGui.QTableWidgetItem(Data['Subject']))
        self.MetaInfoTable.setItem(1, 0, QtGui.QTableWidgetItem(Data['File'].split('\\')[-1]))
        self.MetaInfoTable.setItem(2, 0, QtGui.QTableWidgetItem('%02d/%02d/%02d' % tuple(Data['StartDate'])))
        self.MetaInfoTable.setItem(3, 0, QtGui.QTableWidgetItem(Data['MSN']))
        self.MetaInfoTable.setItem(4, 0, QtGui.QTableWidgetItem(str(Data['Box'])))
        self.MetaInfoTable.setItem(5, 0, QtGui.QTableWidgetItem(str(Data['Group'])))

        self.TWin = [self.TWin1.value(), self.TWin2.value()]

        for k in range(self.naxes):
            plottype = str(self.PlotTypeCombo[k].currentText())
            stim     = str(self.Stim2PlotCombo[k].itemData(self.Stim2PlotCombo[k].currentIndex()).toString())
            curax    = self.main_fig.figure.axes[k]
            sortby   = str(self.SortByCombo[k].currentText())
            ref      = str(self.ReferenceCombo[k].currentText())
            curax.cla()
            
            if   self.OutcomeCombo[k].currentText() == 'All':
                TS = 'StimTS'
            elif self.OutcomeCombo[k].currentText() == 'Hits':
                TS = 'HitsTS'
            elif self.OutcomeCombo[k].currentText() == 'Error':
                TS = 'ErrTS'
            elif self.OutcomeCombo[k].currentText() == 'Miss':
                TS = 'MissTS'
            
            if plottype == 'Later':
                PlotLATER(Data, ShowFit=False, OutlierCalc = 'mc', ax = curax, ShowEq=0)
                curax.set_title('LATER')
                if curax.yaxis_inverted():
                    curax.invert_yaxis()
                
            elif plottype == 'RT Dist':
                if curax.yaxis_inverted():
                    curax.invert_yaxis()
                curax.set_autoscale_on(True)
                
                if self.Stim2PlotCombo[k].currentText() == 'All':
                    stims  = bhv.FindStims(Data)
                else:
                    stims = [str(self.Stim2PlotCombo[k].currentText())]
                whatRT = self.ReferenceCombo[k].currentText()

                for l in stims:
                    
                    if whatRT == 'RT1':
                        RT, _ = bhv.DistXY(Data[l]['StimTS'], Data[l]['NpOut'])
                    elif whatRT == 'RT2':
                        RT = Data[l]['RT2']
                    elif whatRT == 'RT3':
                        RT = Data[l]['RT3']
                    elif whatRT == 'RT4':
                        RT = Data[l]['RT4']
                    elif whatRT == 'RTT':
                        RT = Data[l]['RTT']
                    else:
                        continue

                    color = self.StimColors[int(re.search('[0-9]',l).group())]
                    if self.RTPlotTypeCombo.currentText()=='Hist StepF':
                        if RT.size > 3:
                            curax.hist(RT, self.HistNBins.value(), range=(0, self.HistTWin.value()), lw = 0,
                                       histtype='stepfilled', alpha = 0.8, label = Data[l]['Descr'],
                                       color=color)
                        else:
                            curax.plot([], label='_nolegend_')

                    elif self.RTPlotTypeCombo.currentText()=='Hist Step':
                        if RT.size>3:
                            curax.hist(RT, self.HistNBins.value(), range=(0, self.HistTWin.value()), lw = 3,
                                       histtype='step', label = Data[l]['Descr'],
                                       color=color)
                        else:
                            curax.plot([], label='_nolegend_')

                    elif self.RTPlotTypeCombo.currentText()=='Hist Cum':
                        if RT.size>3:
                            curax.hist(RT, self.HistNBins.value(), cumulative=True, range=(0, self.HistTWin.value()),
                                       lw = 3, histtype='step', label = Data[l]['Descr'], color=color)
                        else:
                            curax.plot([], label='_nolegend_')

                            
                    elif self.RTPlotTypeCombo.currentText()=='PDFs':
                        if RT.size>3:
                            pdf = gaussian_kde(RT)
                            x = np.linspace(0, self.HistTWin.value(), 100)
                            y = pdf(x)
                            y = y/np.sum(y)
                            curax.plot(x, y, lw=5, label = Data[l]['Descr'], color=color)
                        else:
                            curax.plot([], label='_nolegend_')
                    elif self.RTPlotTypeCombo.currentText()=='CDFs':
                        if RT.size>3:
                            pdf = gaussian_kde(RT)
                            x = np.linspace(0, self.HistTWin.value(), 100)
                            y = np.cumsum(pdf(x))
                            y = y/y[-1]
                            curax.plot(x, y, lw=5, label = Data[l]['Descr'], color=color)
                        else:
                            curax.plot([], label='_nolegend_')
                        
                curax.set_title(whatRT)
                curax.set_xlim(0, self.HistTWin.value())
                curax.set_ylabel('Count')
                curax.legend(fancybox=True, prop={'size':10})
                curax.grid()
                            
            elif plottype == 'Behavior':
                
                if ref == 'CentNP':
                    _, indx = bhv.DistYX(Data[stim][TS], Data[stim]['NpIn'])
                    Ref = Data[stim]['NpIn'][indx]
                elif ref == 'Stim':
                    Ref = Data[stim][TS]
                elif ref == 'NpExit':
                    _, indx = bhv.DistXY(Data[stim][TS], Data[stim]['NpOut'])
                    Ref = Data[stim]['NpOut'][indx]
                elif ref == 'RpIn':
                    _, indx = bhv.DistXY(Data[stim][TS], Data[stim][Resp[1]])
                    Ref =  Data[stim][Resp[1]][indx]
                elif ref == '1stLick':
                    _, indx = bhv.DistXY(Data[stim][TS], Data[stim][Resp[0]])
                    Ref =  Data[stim][Resp[0]][indx]

                l = ['No sort','CentNP','Stim','NpExit','RpIn','1stLick','3rdLick']  

                # in case of no sort
                if sortby == 'No sort' or ref==sortby:
                    RT = range(len(Ref))
                                   
                elif re.search('CentNP', sortby):
                    if l.index(ref) < 1:
                        RT, indx = bhv.DistXY(Ref, Data[stim]['NpIn'])
                    elif l.index(ref) > 1:
                        RT, indx = bhv.DistYX(Ref, Data[stim]['NpIn'])
                    
                elif re.search('Stim', sortby):
                    if l.index(ref) < 2:
                        RT , indx = bhv.DistXY(Ref, Data[stim]['StimTS'])
                    elif l.index(ref)>2:
                        RT , indx = bhv.DistYX(Ref, Data[stim]['StimTS'])
                    
                elif re.search('NpExit', sortby):
                    if l.index(ref) < 3:
                        RT, indx = bhv.DistXY(Ref, Data[stim]['NpOut'])
                    elif l.index(ref) > 3:
                        RT, indx = bhv.DistYX(Ref, Data[stim]['NpOut'])

                elif re.search('RpIn', sortby):
                    if l.index(ref) < 4:
                        RT, indx = bhv.DistXY(Ref, Data[stim][Resp[1]])
                    elif l.index(ref) > 4:
                        RT, indx = bhv.DistYX(Ref, Data[stim][Resp[1]])
                        
                elif re.search('1stLick', sortby):
                    if l.index(ref) < 5:
                        RT, indx = bhv.DistXY(Ref, Data[stim][Resp[0]])
                    elif l.index(ref) > 5:
                        RT, indx = bhv.DistYX(Ref, Data[stim][Resp[0]])
                        
                elif re.search('3rdLick', sortby):
                    res = bhv.GetHits(Data[stim][TS], Data[stim][Resp[0]])
                    if l.index(ref) < 6:         
                        RT, indx = bhv.DistXY(Ref, res['ThirdLickHitTS'])
                    elif l.index(ref) > 6:
                        RT, indx = bhv.DistYX(Ref, res['ThirdLickHitTS'])
                        
                s = np.argsort(RT)
		
                '''if sortby == 'RT0/NpExit':
                    ipdb.set_trace()'''
                    
                if re.search('RT0', sortby):
                    RT0, _, _ = bhv.SparseDistance(Data[stim]['StimTS'], Data[stim]['NpIn'],
                                                   direction = 'yx')
                    RT0  = np.round_(RT0, 3)
                    uRT0 = np.sort(np.unique(RT0))
                    s=[]

                    for k in uRT0:
                        indx = np.flatnonzero(RT0 == k)
                        s.extend(indx[np.argsort(RT[indx])])
                        
                yvar = self.YAxesVarCombo.currentText()
                
                Raster3(Ref[s], Data[stim][Resp[0]], yAxVar=yvar,
                        TWin=self.TWin, ax=curax, color=[.5,.5,.5], alpha=0.5, lw = 1) 

                Raster3(Ref[s], Data[stim]['NpIn'],  yAxVar=yvar,
                        TWin=self.TWin, ax=curax, color='b', lw = 2)

                Raster3(Ref[s], Data[stim]['NpOut'], yAxVar=yvar,
                        TWin=self.TWin, ax=curax, color='c', lw = 2)

                Raster3(Ref[s], Data[stim][Resp[1]], yAxVar=yvar,
                        TWin=self.TWin, ax=curax, color='g', lw = 2)
                        
                if Data[stim].has_key('Solnd'):
                    Raster3(Ref[s], Data[stim]['Solnd'], yAxVar=yvar,
                            TWin=self.TWin, ax=curax, color='r', lw=2)

                Raster3(Ref[s], Data[stim]['StimTS'], yAxVar=yvar,
                        TWin=self.TWin, ax=curax, color='k', lw=2)

                curax.set_title(Data[stim]['Descr'])

                if yvar == 'Trial':
                    curax.set_ylim(len(Ref),0)
                    curax.set_ylabel('Trial No')
                elif yvar == 'Time':
                    curax.set_ylim(Ref[-1],0)
                    curax.set_ylabel('Session Time (min)')
                    curax.set_yticks(range(0,3700,600))
                    curax.set_yticklabels(range(0,70,10))
                
            elif plottype=='Perf Evol':
                '''Performs a circular convolution to obtain the hit rate'''
                if curax.yaxis_inverted(): curax.invert_yaxis()
                                
                ktype = str(self.KernelType.currentText())
                ksize = self.KernelSize.value()
                
                if ktype == 'Hamming':
                    kernel  = np.hamming(ksize)
                elif ktype == 'Square':
                    kernel  = np.ones(ksize)

                if self.Stim2PlotCombo[k].currentText() == 'All':
                    stims  = bhv.FindStims(Data)
                else:
                    stims = [self.Stim2PlotCombo[k].itemData(self.Stim2PlotCombo[k].currentIndex()).toString()]
                
                for l in stims:
                    color = self.StimColors[int(re.search('[0-9]',l).group())]
                    
                    if self.Data[l].has_key('HitsTS'):
                        if self.Data[l]['HitsTS'].any() and self.Data[l]['HitsTS'].shape:
                            TS = np.int32(np.round(self.Data[l]['HitsTS']))
                            a  = np.zeros(TS[-1]+1)
                            a[TS]=1
                            b  = np.concatenate((kernel, np.zeros(len(a)-len(kernel)))) # pad with zeros
                            resp = np.fft.ifft(np.fft.fft(a) * np.fft.fft(b)).real # perform the multiplication of the fourier transforms
                            time = np.arange(len(a))/60.0
                            resp = 60*resp/sum(kernel)
                            if time[-1] < 60:
                                time = np.append(time, time[-1]+0.1)
                                time = np.append(time, 60)
                                resp = np.append(resp, 0)
                                resp = np.append(resp, 0)
                            curax.plot(time, resp, color=color, lw=3, label=self.Data[l]['Descr']) # plot
                        else:
                            curax.plot(range(60), np.zeros(60), color=color, lw=3, label=self.Data[l]['Descr'])
                            
                curax.grid(axis='y')
                curax.legend(loc=0, fancybox=True, prop={'size':10})
                curax.set_title('Performance')
                curax.set_xlim(0, 61)
                curax.set_xlabel('Time (min)')
                curax.set_ylabel('Resp/min')

        self.main_fig.figure.tight_layout()
        self.main_fig.figure.canvas.draw()

        if self.SaveFigCheck.checkState()>0:
            self.SaveFig()
            
########################################################################################################################

class LearningCurveGUI(QtGui.QWidget):
    
    def __init__(self, BhvDir = ''):
        self.BhvDir = BhvDir
        QtGui.QWidget.__init__(self)
        self.setWindowTitle("Learning Curve Explorer")
        
        mainLay = QtGui.QHBoxLayout(self)
        
        splitter1 = QtGui.QSplitter(QtCore.Qt.Horizontal)
        frame1 = QtGui.QFrame()
        frame1.setFrameStyle(QtGui.QFrame.StyledPanel)
        
        vLayFrame1 = QtGui.QVBoxLayout(frame1)
        
        # regular expression combo box
        grp = QtGui.QGroupBox('Search Files')
        vLay = QtGui.QVBoxLayout(grp)
        # rat name combobox        
        self.RatNameCombo = QtGui.QComboBox(self)
        self.RatNameCombo.addItems(bhv.GetRatNames(pth = BhvDir))
        self.RatNameCombo.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        hLay = QtGui.QHBoxLayout()
        hLay.setMargin(0)
        hLay.addWidget(QtGui.QLabel('Select Rat Name'))
        hLay.addWidget(self.RatNameCombo)
        vLay.addLayout(hLay)
        
        self.RegExpText = QtGui.QTextEdit()
        self.defaultRegExp = 'NP0[0-4](A)?_(RSip|LSip|CentWht|CWht|[0-9]{1,2}K|NOISE|CLICK)R_(RSip|LSip|CentWht|CWht|[0-9]{1,2}K|NOISE|CLICK)L'
        self.RegExpText.setText(self.defaultRegExp)
        self.RegExpText.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        self.RegExpText.setMaximumHeight(60)
        vLay.addWidget(self.RegExpText)
         # add a 'set default reg exp' button
        self.setDefaultRegExpBtn = QtGui.QPushButton('Reset Reg Exp')
        self.setDefaultRegExpBtn.clicked.connect(self.setDefaultRegExpProc)
        vLay.addWidget(self.setDefaultRegExpBtn)
        
        self.dateCheck = QtGui.QCheckBox('Include Sessions After:')
        self.dateEdit = QtGui.QDateEdit()
        self.dateEdit.setDate(QtCore.QDate(2011,1,1))
        hLay = QtGui.QHBoxLayout()
        hLay.setMargin(0)
        hLay.addWidget(self.dateCheck)
        hLay.addWidget(self.dateEdit)
        vLay.addLayout(hLay)
        
        self.searchFilesBtn = QtGui.QPushButton('Search Files')
        self.searchFilesBtn.clicked.connect(self.searchFilesProc)
        self.searchFilesBtn.setStyleSheet('QPushButton{background-color: rgba(243,134,48)}')
        vLay.addWidget(self.searchFilesBtn)
                
        
        self.filesTable = QtGui.QTableWidget(0, 1)
        self.filesTable.verticalHeader().setVisible(False)
        self.filesTable.setHorizontalHeaderLabels(['Include Files'])
        self.filesTable.setColumnWidth(0, 300)
        #self.filesTable.horizontalHeader().setVisible(False)
        self.filesTable.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        vLay.addWidget(self.filesTable)
        
        '''self.selectBtn = QtGui.QPushButton('Include All')
        self.selectBtn.setCheckable(True)
        self.selectBtn.clicked.connect(self.selectProc)

        self.selectInvertBtn = QtGui.QPushButton('Toggle Sel')
        self.selectInvertBtn.clicked.connect(self.selectInvertProc)
        hLay = QtGui.QHBoxLayout()
        hLay.setMargin(0)
        hLay.addWidget(self.selectBtn)
        hLay.addWidget(self.selectInvertBtn)
        vLay.addLayout(hLay)'''
        
        vLayFrame1.addWidget(grp)
        
        # include RTs checkbox
        grp = QtGui.QGroupBox('Plot options')
        
        grpLay = QtGui.QGridLayout(grp)
        
        # Raw vs proportion plot
        self.PlotOptionCheck = QtGui.QCheckBox('Raw / Proportion')
        self.PlotOptionCheck.setChecked(True)
        grpLay.addWidget(self.PlotOptionCheck, 0, 0, 1, 4)
        
        self.IncRTCheck = []
        
        self.IncRTCheck.append(QtGui.QCheckBox('RT1'))
        self.IncRTCheck[-1].setChecked(True)
        grpLay.addWidget(self.IncRTCheck[-1], 1, 0)
        
        self.IncRTCheck.append(QtGui.QCheckBox('RT2'))
        self.IncRTCheck[-1].setChecked(True)
        grpLay.addWidget(self.IncRTCheck[-1], 1, 1)
        
        self.IncRTCheck.append(QtGui.QCheckBox('RT3'))
        self.IncRTCheck[-1].setChecked(True)
        grpLay.addWidget(self.IncRTCheck[-1], 1, 2)
        
        self.IncRTCheck.append(QtGui.QCheckBox('RT4'))
        self.IncRTCheck[-1].setChecked(True)
        grpLay.addWidget(self.IncRTCheck[-1], 1, 3)
        
        # add a plot button
        self.PlotBtn = QtGui.QPushButton('Plot')
        self.PlotBtn.clicked.connect(self.Plot)
        self.PlotBtn.setStyleSheet('QPushButton{background-color: rgba(0,190,0)}')
        grpLay.addWidget(self.PlotBtn, 2, 0, 1, 2)
        
        self.saveFigBtn = QtGui.QPushButton('Save Figure')
        self.saveFigBtn.clicked.connect(self.saveFigProc)
        self.saveFigBtn.setStyleSheet('QPushButton{background-color: rgba(8,107,154)}')
        grpLay.addWidget(self.saveFigBtn, 2, 2, 1, 2)
        
        vLayFrame1.addWidget(grp)
        #vLayFrame1.addStretch(1)
        
        # add a settings button
        self.settingsBtn = QtGui.QPushButton('Settings')
        self.settingsBtn.clicked.connect(self.settingsProc)
        vLayFrame1.addWidget(self.settingsBtn)
        
        splitter1.addWidget(frame1)        
        
        frame2 = QtGui.QFrame()
        frame2.setFrameStyle(QtGui.QFrame.StyledPanel)
        vLayFrame2 = QtGui.QVBoxLayout(frame2)
        vLayFrame2.setMargin(0)
        vLayFrame2.setSpacing(0)
        
        # add a figure and a toolbar
        self.MainFig = MplWidget()
        self.Ntb = NavToolbar(self.MainFig, self)
        self.Ntb.setIconSize(QtCore.QSize(15,15))
        #vLayFrame2.addLayout(gLay)
        vLayFrame2.addWidget(self.MainFig)
        vLayFrame2.addWidget(self.Ntb)
        
        splitter1.addWidget(frame2)
        
        mainLay.addWidget(splitter1)
        # if running in linux set a certain style for the buttons and widgets
        if sys.platform == 'linux2':
            QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Plastique'))

        # finally show the entire gui            
        self.show()
        
    def setDefaultRegExpProc(self):
        self.RegExpText.setText(self.defaultRegExp)
        
    def selectProc(self):
        if self.filesTable.rowCount() == 0:
            return
            
        state = self.selectBtn.isChecked()
        
        if state:
            self.selectBtn.setText('All Selected')
        else:
            self.selectBtn.setText('None Selected')
        
        for n in self.includeFileChecks:
            n.setChecked(state)
            
    def selectInvertProc(self):
        if self.filesTable.rowCount() == 0:
            return
        
        for n in self.includeFileChecks:
            n.setChecked(not n.checkState())
        
    def searchFilesProc(self):
        
        from datetime import date
        
        RatName = str(self.RatNameCombo.currentText())
        Task    = str(self.RegExpText.toPlainText())
        self.files, MSN = bhv.GetFilenames(RatName, Task, self.BhvDir)
        
        for n in reversed(range(self.filesTable.rowCount())):
            self.filesTable.removeRow(n)
            
        # delete checkboxes
        if hasattr(self, 'includeFileChecks') and len(self.includeFileChecks) > 0:
            for n in self.includeFileChecks:
                n.deleteLater()
                
        self.includeFileChecks = []
        
        if self.dateCheck.checkState():
            dateRef = self.dateEdit.date().toPyDate()
        else:
            dateRef = date(2000,1,1)
            
        for n, f in enumerate(self.files):
            self.filesTable.insertRow(n)            
            f = os.path.split(f)[1]
            self.includeFileChecks.append(QtGui.QCheckBox(f))
            self.includeFileChecks[-1].setFont(QtGui.QFont('', pointSize=8))
            self.filesTable.setCellWidget(n,0, self.includeFileChecks[-1])
            self.filesTable.setRowHeight(n, 20)
            
            dateStr = re.search('[0-9]{6}', f).group()
            dateFile = date(int(dateStr[0:2]) + 2000, int(dateStr[2:4]), int(dateStr[4:]))
                
            if dateFile >= dateRef:
                self.includeFileChecks[-1].setChecked(True)            
            
        self.filesTable.setAlternatingRowColors(True)
    
    def saveFigProc(self):
        fname = str(QtGui.QFileDialog.getSaveFileName(caption = 'Save a copy of the figure ...'))
        
        if fname:
            self.MainFig.figure.savefig(fname, dpi = 300, format = 'svg')
        
    def settingsProc(self):
        if settings.edit():
            pass
        
    def Plot(self):
        
        if self.filesTable.rowCount() == 0:
            self.searchFilesProc()
            
        RatName = str(self.RatNameCombo.currentText())
        Task    = str(self.RegExpText.toPlainText())
        files, MSN = bhv.GetFilenames(RatName, Task, self.BhvDir)
        CalcType   = self.PlotOptionCheck.checkState()
        IncRTs     = self.IncRTCheck[0].checkState()
        
        files = []
        for k in self.includeFileChecks:
            if k.checkState():
                files.append(os.path.join(self.BhvDir, str(k.text())))
        
        # create dictionaries to hold hits, miss and error, and colors
        keys = ['Tone1', 'Tone2', 'RSipLight', 'LSipLight', 'WhtFLight', 'Catch']
        col  = ['b', 'g', 'b', 'g', 'r', [.5, .5, .5]]
        tmp = {}; colors = {}        
        for c, k in zip(col, keys):
            tmp[k]    = []
            colors[k] = c

        # create a dictionary for every type of rection time
        RTs = ['RT1','RT2','RT3']
        RTDict = {}
        for k in RTs: RTDict[k] = {}
            
        # create lists to hold data
        MSN = []; dates = []; R = []; L = []
        
        if RatName == '':
            initialDate = dtm(year = 2012, month = 11, day = 29)
            exclude = [dtm(2012, 12, 13), dtm(2013, 4, 1)]
        else:
            initialDate = dtm(1,1,1)
            exclude = []
        
        # iterate over files to extract the information
        for j,k in enumerate(files):
            
            Data  = bhv.LoadBehFile(k)  # load the data
            Stims = bhv.FindStims(Data) # find stimuli
            date  = Data['StartDate']   # get session date and transform it to a date format
            date  = dtm(year = date[2]+2000, month = date[0], day = date[1])
        
            # this is if one wants to exlude some sessions, or include sessions
            # after a certain date
            if date < initialDate or date in exclude: continue
        
            # store date and trining protocol
            dates.append(Data['StartDate'])
            MSN.append(Data['MSN'])
        
            # iterate over stimuli to extract RTs and performance info
            for m, n in enumerate(Stims):
                
                # get the number of errors
                if Data[n].has_key('ErrTS'):
                    eSize = np.array(Data[n]['ErrTS']).size
                else:
                    eSize = 0
        
                if Data[n]['Descr'] not in tmp:
                    tmp[Data[n]['Descr']] = []
                # append the number of hits, miss and errors
                tmp[Data[n]['Descr']].append([j,
                                              np.array(Data[n]['HitsTS']).size,
                                              eSize,
                                              np.array(Data[n]['MissTS']).size])
        
                # Extract RTs
                # create a key inside each dictionary if it doesn't exists
                key = Data[n]['Descr']
                for r in RTs:
                    if not RTDict[r].has_key(key):
                        RTDict[r][key] = []
        
                    # create a key to hold the session number for a particular RT
                    if not RTDict[r].has_key('x_'+key):
                        RTDict[r]['x_'+key] = []
        
                    if key != 'Catch':
                        # get the reaction times
                        # this is to eliminate weird RTs coming from a bug in my Med-PC code
                        Data[n][r] = np.array(Data[n][r])
                        if Data[n][r].size > 1 and np.any( Data[n][r] > 0.3 ):
                            RT = Data[n][r][ Data[n][r] < 3.0 ]
                        else:
                            RT = Data[n][r]
                    else:
                        RT, x,y = bhv.SparseDistance(Data[n]['StimTS'], Data[n]['NpOut'])
                        
                    RTDict[r][key].append(RT)
        
                    # get the session number
                    RTDict[r]['x_'+key].append(j)

        
            # get the catchs for the left and the right
            lCatch = [k for k in bhv.FindStims(Data) if Data[k]['Descr'] == 'Catch']
            if len(lCatch) > 0:
                n = lCatch[0]
                RLick = Data['EventTS'][np.flatnonzero(Data['EventName'] == 'RightLickOn')][0]
                LLick = Data['EventTS'][np.flatnonzero(Data['EventName'] == 'LeftLickOn')][0]
                
                R.append([j,
                          Data[n]['StimTS'].size,
                          bhv.GetHits(Data[n]['StimTS'], RLick)['StimHitsIndx'].size])
                L.append([j,
                          Data[k]['StimTS'].size,
                          bhv.GetHits(Data[n]['StimTS'], LLick)['StimHitsIndx'].size])
            else:
                R.append([j, np.nan, np.nan])
                L.append([j, np.nan, np.nan])
        
        
        #v = 1000
        for k in tmp.keys():
            tmp[k] = np.array(tmp[k])
        '''    if tmp[k].size > 0 and v > tmp[k][:,0].min():
                v = tmp[k][:,0].min()'''
                
        for k in tmp.keys():
            if tmp[k].size > 0:
                tmp[k][:,0] = tmp[k][:,0]# - v
        
        R = np.array(R); R[:,0] = R[:,0]#-v
        L = np.array(L); L[:,0] = L[:,0]#-v
        
        # set some parameters
        s = range(len(dates)); lw = 6; ms = 10; a = 0.7
        
        # clear the figure and create new axes
        fig = self.MainFig.figure
        fig.clf()
        ax  = fig.add_subplot(111)
        
        
        # change the labels to reflect the stimulus - outcome association                
        for k in keys[0:-1]:
            t = tmp[k]
            if not np.all(np.isnan(t)):
                if k == 'Tone1':
                    label = k+'--> Right'
                elif k == 'Tone2':
                    label = k+'--> Left'
                elif k == 'WhtFLight':
                    label = k+'--> Right'
                else:
                    label = k
        
                if CalcType:
                    ax.plot(t[:,0], t[:,1]/np.float32((t[:,1]+t[:,2])+t[:,3]),'-o',
                            ms=ms, mew = 0, lw=lw, alpha=a, label=label)
                else:
                    ax.plot(t[:,0], t[:,1],'-o',
                            ms=ms, mew = 0, lw=lw, alpha=a, label=label)

        
        # plot the catch trials as a proportion
        if CalcType:
            
            ax.plot(R[:,0], R[:,2]/np.float32(R[:,1]), '--o', color = 'b',
                    ms=ms, mew = 0, lw=lw, alpha=a, label='catch R')
            ax.plot(L[:,0], L[:,2]/np.float32(L[:,1]), '--o', color = 'g',
                    ms=ms, mew = 0, lw=lw, alpha=a, label='catch L')
        
        # ... or as raw data
        else:
            ax.plot(R[:,0], R[:,2], '--o', color = 'b',
                    ms=ms, mew = 0, lw=lw, alpha=a, label='catch R')
            ax.plot(L[:,0], L[:,2], '--o', color = 'g',
                    ms=ms, mew = 0, lw=lw, alpha=a, label='catch L')
        
        # put a dotted line before the first session of the white light
        if not np.all(np.isnan(tmp['WhtFLight'])):
            ax.axvline(tmp['WhtFLight'][0][0]-0.5, color='k', alpha=0.5, linestyle='--', lw=3)
            
        
        # title and labels
        ax.set_title('%s %s' % (RatName, Task))
        ax.set_xlabel('Session Date')
        
        # y label for proportion of trials
        if CalcType:
            ax.set_ylabel('Proportion of Trials')
            ax.set_ylim(-0.05, 1.05)
            ax.axhline(0.5, 0, 10, lw = 2, linestyle = '--', color='k', alpha = 0.5)
            ax.axhline(1.0, 0, 10, lw = 2, linestyle = '--', color='k', alpha = 0.5)
            ax.axhline(0.0, 0, 10, lw = 2, linestyle = '--', color='k', alpha = 0.5)
        else:
            ax.set_ylabel('Number of Trials')
        
        # set the x axis tick marks to be the training dates
        ax.set_xticks(s)
        ax.set_xticklabels(['%02d-%02d-%d' % (k[0],k[1],k[2]) for k in dates],
                           rotation = 30, horizontalalignment = 'right', fontsize = 9)
        
        # iterate over tick labels, if monday change to bold red
        for j,k in zip(ax.get_xmajorticklabels(), dates):
            if dtm(k[2],k[0],k[1]).isoweekday() == 1: # check whether it is a monday
                j.set_color('r')
                j.set_weight('bold')

        
        # if selected, plot RT1 for each stimuli
        # ... I should add an option to select different RTs ...
        if IncRTs:
            ax2 = ax.twinx() # create a y axis
            r = 'RT1'
            for k in RTDict[r]:
                if k.find('x_') != 0:
                    x = RTDict[r]['x_'+k]# - v
                    ax2.plot(x ,[np.array(n).mean() for n in RTDict[r][k]],
                                  color = colors[k], lw = 3, label = k)
                    ci = np.array([bootci(n, 100) for n in RTDict[r][k]])
                    ax2.fill_between(x, ci[:,0], ci[:,1], alpha = 0.5, color = colors[k])
            ax2.set_ylabel('Time (sec)')
            ax2.set_ylim(0, 1.0)
        
        # grid and limits
        ax.grid(axis = 'both')
        ax.set_xlim(-1, s[-1]+1)
        
        # make the first axes' background transparent ad raise it
        #ax.set_axis_bgcolor('none')
        #ax.set_zorder(10)
        
        # handler for the legend. The legend corresponds to the first axes
        #leg = ax.legend(loc='lower left', shadow=True, prop={'size':12}, fancybox = True)
        #leg.set_zorder(100)
        #leg.draggable()
        
        fig.tight_layout()
        fig.canvas.draw()
        
########################################################################################################################        


def BilateralTask(pth):

    import datetime as dtime
    
    ratNames = bhv.GetRatNames(prefix = 'HMV', pth = pth)

    class SetParams(dt.DataSet):
        RatName  = di.ChoiceItem('Choose a Rat', ratNames )
        Task     = di.StringItem('Regular Expr', 'UNBIAS|6KR_9KL|2T')
        PlotMSN  = di.BoolItem('PlotMSN?').set_pos(col=0)
        PlotEq   = di.BoolItem('Plot Equations?').set_pos(col=1)
        SaveFig  = di.BoolItem('Save Figure', default=False)
        BiasIndx = di.ChoiceItem('Discr Indx', [(1,'Simple'),(2,'Complex'),])
        
    Params=SetParams()
    if Params.edit()==1:
        RatName    = ratNames[Params.RatName]
        Task       = Params.Task
        files, MSN = bhv.GetFilenames(RatName, RegExp=Task, BhvDir = pth)

        Dates=[]
        # build a dictionary to hold all the data
        Stims = {}; Stims['Stim0'] = {}; Stims['Stim1'] = {}
        params = ['mRTT', 'mRT1', 'mRT2', 'mRT3', 'eRTT', 'eRT1', 'eRT2', 'eRT3', 'n']
        for p in params:
            Stims['Stim0'][p] = []
            Stims['Stim1'][p] = []

        # iterate over files
        for k in files:
            Data = bhv.LoadBehFile(k)
            Dates.append(np.int32(Data['StartDate']))
            if not bhv.FindStims(Data):
                Data = bhv.GetBhvParams(Data)

            for s in ['Stim0', 'Stim1']:
                if Data.has_key(s):
                    for p in params:
                        if p == 'n':
                            if Data[s].has_key('RTTe'):
                                rE = Data[s]['RTTe'].size
                            else:
                                rE = 0
                            Stims[s][p].append([Data[s]['RTT'].size, rE])
                        if Data[s]['RTT'].size > 1:
                            if p[0] == 'm':
                                Stims[s][p].append(Data[s][p[1:]].mean())
                            elif p[0] =='e':
                                Stims[s][p].append(bootci(Data[s][p[1:]], 100))
                        else:
                            if p[0] == 'm':
                                Stims[s][p].append(0)
                            elif p[0] =='e':
                                Stims[s][p].append([0,0])

        # transform everything into an array
        for p in params:
            Stims['Stim0'][p] = np.array(Stims['Stim0'][p])
            Stims['Stim1'][p] = np.array(Stims['Stim1'][p])

        mn = ['mRT1','mRT2','mRT3']
        er = ['eRT1','eRT2','eRT3']
        for m, e in zip(mn,er):
            Stims['Stim0'][e] = [Stims['Stim0'][m]-Stims['Stim0'][e][:,0]\
                                ,Stims['Stim0'][e][:,1]-Stims['Stim0'][m]]
            Stims['Stim1'][e] = [Stims['Stim1'][m]-Stims['Stim1'][e][:,0]\
                                ,Stims['Stim1'][e][:,1]-Stims['Stim1'][m]]

        # Calculate task dates and MSNs
        for j, k in enumerate(MSN):
            if re.search('(?<=NP04A_).*',k):
                MSN[j]=re.search('(?<=NP04A_).*', k).group()

        Dates   = np.array(Dates)
        Mondays = np.array([], dtype = np.int32)

        # determine if that date is a monday
        for j,k in enumerate(Dates):
            if dtime.date(k[2],k[0],k[1]).isoweekday()==1:
                Mondays = np.append(Mondays, j)
                
        Dates2 = ['%02d/%02d' % (k[0],k[1]) for k in Dates]
 
        gAlpha = 0.6; orig = 0; a = 0.6; h = 0.6; ErrLine = 3
        col = [ [0,.2,.7] , [0,.4,0] , [.8,.1,.1] ]
        yl  = [-1, max([Stims['Stim0']['mRTT'].size, Stims['Stim0']['mRTT'].size])]

        plt.figure()
        x  = np.arange(1, Stims['Stim0']['mRT1'].size+1)
        for m, e, c in zip(mn, er, col):
            plt.barh(x, Stims['Stim0'][m], left=orig, height=h, xerr=Stims['Stim0'][e],
                     color=c, alpha=a, align='center', edgecolor='', label = m,
                     error_kw={'elinewidth':ErrLine,'ecolor':c})
            orig = orig+Stims['Stim0'][m]

        x  = np.arange(1, Stims['Stim1']['mRT1'].size+1)
        orig = 0
        for m, e, c in zip(mn, er, col):
            plt.barh(x, -Stims['Stim1'][m], left=orig, height=h, xerr=Stims['Stim1'][e],
                     color=c, alpha=a, align='center', edgecolor='', label = '_nolegend_',
                     error_kw={'elinewidth':ErrLine,'ecolor':c})
            orig = orig-Stims['Stim1'][m]

        plt.plot(np.zeros(len(x[Mondays])),x[Mondays],'o',color=[1,1,0], label='Monday')
        plt.legend(fancybox = True)#,  mode='expand', ncol=4, prop = {'size':9},loc=8)
        
        plt.title(RatName+' '+Task)
        xl = 1.25*(max(max(Stims['Stim0']['mRTT']), max(Stims['Stim1']['mRTT'])))
        plt.xlim(-xl,xl)    

        # Draw a textbox with the name of the task
        if Params.PlotMSN == 1:
            bp = dict(boxstyle="round", alpha=0.7, fc='w', ec=[.5,.5,.5])
            for j,k in enumerate(MSN):
                plt.text(-xl+0.1, x[j], k, va='center', fontsize=7, bbox=bp)
                
        # Add 'L' and 'R' arrows
        bp = dict(boxstyle="LArrow", alpha=0.7, fc='w', ec=[.5,.5,.5])
        plt.text(-0.85*xl, 0.5+yl[1]/2.0, 'L', ha='center', va='center', fontsize=15, bbox=bp)
        bp['boxstyle']="RArrow"
        plt.text( 0.85*xl, 0.5+yl[1]/2.0, 'R', ha='center', va='center', fontsize=15, bbox=bp)
        plt.grid(alpha = gAlpha)
        plt.ylim(0, len(x)+1)
        plt.ylabel('Session Date')
        plt.yticks(x, Dates2, fontsize=8)
        plt.xlabel('Time (sec)')

        #######################################################################################################
        ##  Plots the number of trials for corrects, incorrects for each session

        col=[[0,.2,.7],[0,.4,0],[.8,.1,.1]]
        a=0.8

        plt.subplot(2,2,2)

        x  = np.arange(1, len(Stim0)+1)
        b1 = plt.barh(x, nR[:,0], height=h, left=0,       align='center', color=col[0], alpha=a)
        b2 = plt.barh(x, nR[:,1], height=h, left=nR[:,0], align='center', color=col[2], alpha=a)

        x  = np.arange(1, len(Stim1)+1)
        b1 = plt.barh(x, -nL[:,0], height=h, left=0,        align='center', color=col[0], alpha=a)
        b2 = plt.barh(x, -nL[:,1], height=h, left=-nL[:,0], align='center', color=col[2], alpha=a)

        plt.plot([0,0],[-1,yl[1]+2],'k',linewidth=2)
        p1, = plt.plot(np.zeros(len(x[Mondays])), x[Mondays], 'o', color=[1,1,0])

        plt.legend([b1,b2,p1],['Corrects','Incorrects','Monday'],mode='expand',
                   loc=8, ncol=3, fancybox=True, prop={'size':9})

        plt.title('Number of trials for '+RatName+' '+Task, fontdict={'fontsize':10})
        plt.xlabel('Number of trials')
        plt.yticks(x, '')
        plt.grid(alpha = gAlpha)
        plt.tight_layout()
        xl = max(np.append(np.sum(nR,axis=1), np.sum(nL,axis=1)))

        plt.ylim(yl)
        plt.xlim(-1.1*xl,1.1*xl)

        # Add 'L' and 'R' arrows
        bp = dict(boxstyle="LArrow", alpha=0.7, fc='w', ec=[.5,.5,.5])
        plt.text(-0.85*1.1*xl, 0.5+yl[1]/2., 'L', ha='center', va='center', fontsize=15, bbox=bp)
        bp = dict(boxstyle="RArrow", alpha=0.7, fc='w', ec=[.5,.5,.5])
        plt.text( 0.85*1.1*xl, 0.5+yl[1]/2., 'R', ha='center', va='center', fontsize=15, bbox=bp)

        #######################################################################################################
        ##  Plots the discriminability index

        plt.subplot(2,2,3)

        dR=nR[:,0]/(nR[:,0]+nR[:,1])
        dL=nL[:,0]/(nL[:,0]+nL[:,1])

        x   = np.arange(1, max([len(Stim0), len(Stim1)])+1)
        b1, = plt.plot(dR, range(1,dR.size+1), linewidth=2, color=col[0], marker='o', mec=col[0], ms=6)
        b2, = plt.plot(dL, range(1,dL.size+1), linewidth=2, color=col[1], marker='o', mec=col[1], ms=6)
        p1, = plt.plot(np.zeros(len(x[Mondays])), x[Mondays], 'o', color=[1,1,0])

        plt.legend([b2,b1,p1], ['Left','Right','Monday'], mode='expand',
                   ncol=3, loc=8, fancybox=True, prop={'size':9})

        if Params.PlotEq:
            bp = dict(boxstyle="round", alpha=0.7, fc='w', ec=[.5,.5,.5])
            plt.text(0.5, yl[1]*0.9, r'$\frac{Corrrects}{Corrects+Incorrects}$',
                     ha='center', va='center', fontsize=14, bbox=bp)

        if PlotMSN==1:
            bp = dict(boxstyle="round", alpha=0.7, fc='w', ec=[.5,.5,.5])
            for j,k in enumerate(MSN):
                plt.text(0.025,x[j],k,va='center',fontsize=7, bbox=bp)

        plt.xlabel('Discriminability Index')
        plt.xlim(-0.1,1.1)
        plt.ylim(yl)
        plt.ylabel('Session Date')
        plt.yticks(x,Dates2,fontsize=8)
        plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0] , [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        plt.grid(alpha = gAlpha)

        #######################################################################################################
        ##  PLots the bias index

        plt.subplot(2,2,4)

        if Params.PlotEq:
            if Params.BiasIndx==2:
                eq=r'$\frac{Cor_{Right}-Cor_{Left}}{Cor_{Right}+Inc_{Right}+Cor_{Left}+Inc_{Left}}$'
            elif Params.BiasIndx==1:
                eq=r'$\frac{Cor_{Right}}{Cor_{Right}+Inc_{Right}}-\frac{Cor_{Left}}{Cor_{Left}+Inc_{Left}}$'
                
            bp = dict(boxstyle="round", alpha=0.7, fc='w', ec=[.5,.5,.5])
            plt.text(-0.5, yl[1]*0.9, eq, ha='center', va='center', fontsize=14, bbox=bp)

        if Params.BiasIndx==2:
            dI=((nR[:,0]-nR[:,1])-(nL[:,0]-nL[:,1]))/(nR[:,0]+nR[:,1]+nL[:,0]+nL[:,1]+1e-3)
        elif Params.BiasIndx==1:
            dI=dR-dL
            
        b1  = plt.barh(x,  dI, height=h, align='center',color=col[0],alpha=a)
        plt.plot([0,0],[-1,len(x)+2],'k',linewidth=2)
        p1, = plt.plot(np.zeros(len(x[Mondays])),x[Mondays],'o',color=[1,1,0])
        plt.xlim(-1,1)
        plt.ylim(yl)
        plt.yticks(x,'')
        plt.grid(alpha = gAlpha)

        # Add 'L' and 'R' arrows
        bp = dict(boxstyle="LArrow", alpha=0.7, fc='w', ec=[.5,.5,.5])
        plt.text(-0.8, 0.5+yl[1]/2., 'L', ha='center', va='center', fontsize=15, bbox=bp)
        bp = dict(boxstyle="RArrow", alpha=0.7, fc='w', ec=[.5,.5,.5])
        plt.text( 0.8, 0.5+yl[1]/2., 'R', ha='center', va='center', fontsize=15, bbox=bp)

        plt.tight_layout()

        #######################################################################################################
        ##  Save the figure

        if Params.SaveFig:
            SaveFigure(FigHandles)

########################################################################################################################

def Sessions_NHits_NMiss_NError(Stims=[], horizontal=True, labels=[], axes=[],
                                XYTicks=[True,True], XYLabels=[True,True], colors=[ [0,.2,.7], [0,.4,0], [.8,.1,.1] ]):

    '''Helper function to plot the Sessions numbers of corrects, incorrects and misses as an stacked barplot
       INPUT:
           -Sesssions:  [numpy array] containing one column for each number of hits, miss and errors,
                        and rows for each session
           -horizontal: whether to plot vertical or horizontal bars.
           -labels:     labels of each column
           -axes:       provide an axes handle if you want to embed it in a figure
           -XYTicks:    whether to draw the x and y tickmarks. Default: True
           -XYLabels:   whether to draw the x and y labels. Default: True'''

    BarAlpha  = 0.8
    BarHeight = 0.8
    GridAlpha = 0.6

    
    if not axes or not isinstance(axes, matplotlib.axes.Axes):
        ax=plt.subplot(111)
        rc('font', size=9, family='monospace', serif='Bitstream Vera Sans Mono')

    if horizontal:
        if Stim0:
            x  = np.arange(1, len(Stim0)+1)
            b1 = plt.barh(x, nR[:,0], height=BarHeight, left=0,       align='center',color=col[0], alpha=BarAlpha)
            b2 = plt.barh(x, nR[:,1], height=BarHeight, left=nR[:,0], align='center',color=col[2], alpha=BarAlpha)

        if Stim1:
            x  = np.arange(1, len(Stim1)+1)
            b1 = plt.barh(x, -nL[:,0], height=BarHeight, left=0,        align='center', color=col[0], alpha=BarAlpha)
            b2 = plt.barh(x, -nL[:,1], height=BarHeight, left=-nL[:,0], align='center', color=col[2], alpha=BarAlpha)

        plt.plot([0,0],[-1,yl[1]+2],'k',linewidth=2)
        p1, = plt.plot(np.zeros(len(x[Mondays])), x[Mondays], 'o', color=[1,1,0])

        # Add 'L' and 'R' arrows
        bp = dict(boxstyle="LArrow", alpha=0.7, fc='w', ec=[.5,.5,.5])
        plt.text(-0.85*1.1*xl, 0.5+yl[1]/2., 'L', ha='center', va='center', fontsize=15, bbox=bp)

        bp = dict(boxstyle="RArrow", alpha=0.7, fc='w', ec=[.5,.5,.5])
        plt.text( 0.85*1.1*xl, 0.5+yl[1]/2., 'R', ha='center', va='center', fontsize=15, bbox=bp)
        
    else:
        pass
    
    plt.legend([b1,b2,p1],['Corrects','Incorrects','Monday'],mode='expand',
               loc=8, ncol=3, fancybox=True, prop={'size':9})

    plt.title('Number of trials for '+RatName+' '+Task, fontdict={'fontsize':10})
    plt.xlabel('Number of trials')
    plt.yticks(x, '')
    plt.grid(alpha=GridAlpha)
    plt.tight_layout()
    xl = max(np.append(np.sum(nR,axis=1), np.sum(nL,axis=1)))

    plt.ylim(yl)
    plt.xlim(-1.1*xl,1.1*xl)

########################################################################################################################
        
def PlotBeh(Data):
    import matplotlib.pyplot as p

    Stims = [k for k in Data.keys() if k.find('Stim')!=-1]

    rc('font',size=8)
    rc('font',family='serif')

    if bhv.GetMapping(Data)==1:

        p.figure(facecolor='w',dpi=100,figsize=(8,3.5))
        RTT = Stim['RTT']
        sRT = np.argsort(RTT)
        RT0 = Stim['RT0'][sRT]
        RT1 = Stim['RT1'][sRT]
        RT2 = Stim['RT2'][sRT]
        RT3 = Stim['RT3'][sRT]

        nTrials = len(Stim['RTT'])
        y = range(nTrials)

        p.plot(-RT0,y,'c.',markeredgecolor='c')
        p.plot(np.zeros(nTrials),y,'k.',markeredgecolor='k')
        p.plot(RT1,y,'b.',markeredgecolor='b')
        p.plot(RT1+RT2,y,'r.',markeredgecolor='r')
        p.plot(RT1+RT2+RT3,y,'g.',markeredgecolor='g')
        p.xlim(-1,np.max(RTT)+0.5)
        p.title('Reaction times')
        p.show()

    elif GetMapping(Data)==2:
        p.figure(facecolor='w',dpi=100,figsize=(8,3.5))
        j=1
        for k in Stims:

            p.subplot(2,2,j)

            if Data[k].has_key('RTTc') and Data[k].has_key('RT0c') and Data[k].has_key('RT1c') and Data[k].has_key('RT2c') and Data[k].has_key('RT3c'):
                RTT=Data[k]['RTTc']
                sRT=np.argsort(RTT)
                nTrials=len(Data[k]['RTTc'])
                y=range(nTrials)
                p.plot(np.zeros(nTrials),y,'k.',markeredgecolor='k')

                RT0 = Data[k]['RT0c'][sRT]
                RT1 = Data[k]['RT1c'][sRT]
                RT2 = Data[k]['RT2c'][sRT]
                RT3 = Data[k]['RT3c'][sRT]

                p.ylim=(-1,RTT.size*1.1)
                p.plot(-RT0,y,'c.',markeredgecolor='c')
                p.plot(RT1,y,'b.',markeredgecolor='b')
                p.plot(RT1+RT2,y,'r.',markeredgecolor='r')
                p.plot(RT1+RT2+RT3,y,'g.',markeredgecolor='g')
                p.xlim(-1,np.max(RTT)+0.5)
                p.grid(True)
                p.title('Reaction times: Corrects '+str(Data[k]['Descr']))
                p.xlabel('Time (sec)',fontsize=10)
            j=j+1

            p.subplot(2,2,j)

            if Data[k].has_key('RTTi') and Data[k].has_key('RT0i') and Data[k].has_key('RT1i') and Data[k].has_key('RT2i') and Data[k].has_key('RT3i'):
                RTT=Data[k]['RTTi']
                sRT=np.argsort(RTT)
                nTrials=len(Data[k]['RTTi'])
                y=range(nTrials)
                p.plot(np.zeros(nTrials),y,'k.',markeredgecolor='k')

                RT0=Data[k]['RT0i'][sRT]
                RT1=Data[k]['RT1i'][sRT]
                RT2=Data[k]['RT2i'][sRT]
                RT3=Data[k]['RT3i'][sRT]

                p.ylim=(-1,RTT.size*1.1)
                p.plot(-RT0,y,'c.',markeredgecolor='c')
                p.plot(RT1,y,'b.',markeredgecolor='b')
                p.plot(RT1+RT2,y,'r.',markeredgecolor='r')
                p.plot(RT1+RT2+RT3,y,'g.',markeredgecolor='g')
                p.xlim(-1,np.max(RTT)+0.5)
                p.grid(True)
                p.title('Reaction times: Incorrects '+str(Data[k]['Descr']))
                p.xlabel('Time (sec)',fontsize=10)
            j=j+1

########################################################################################################################

def PlotRTStats(Stim):
    import matplotlib.pyplot as p

    rc('font',size=8)
    rc('font',family='serif')

    for x in Stim.keys():
        if re.search('RT[0-9][c,i]', x):
            m=2
            break
        else:
            m=1

    if m==1:
        keys = ['RT1','RT2','RT3']
    elif m == 2:
        keys = ['RT1c','RT2c','RT3c']
        
    fig = p.figure(facecolor='w', dpi=100, figsize=(8,3.5))

    for j, k in enumerate(keys):
        ax = fig.add_subplot(3, 2, 2*(j+1)-1)
        ax.hist(Stim[k], 50, histtype='stepfilled')
        ax.set_xlim(-.5,1.5)
        ax.grid(True)
        
        ax = fig.add_subplot(3, 2, 2*(j+1))
        ax.boxplot(Stim[k],notch=1, sym='')
        ax.grid(True)
                
    fig.show()

########################################################################################################################

def SaveFigure(FigHandles=[],FigName='Fig_',dpi=300, SaveDir='', Format='.jpg'):

    '''Helper function to save figures'''
    if not isinstance(FigHandles,list): FigHandles=[FigHandles]
    if not FigHandles or not [isinstance(k,Figure) for k in FigHandles]:
        import matplotlib._pylab_helpers
        FigHandles=[manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        if FigHandles:
            class SetParams(dt.DataSet):
                '''Select the paramters to save the currently active figure'''
                Figures = di.MultipleChoiceItem(label='Select a Figure',
                                                choices=[str(k.number) for k in FigHandles])
                
                FigName = di.StringItem(label='FigName',default='Fig_').set_pos(col=0)
                
                Format  = di.ChoiceItem(label='Format',
                                     choices=['.jpg','.png','.svg','.pdf']).set_pos(col=1)
                
                ImgDir  = di.DirectoryItem(label='Directory 2 save').set_pos(col=0)
                
                dpi     = di.IntItem(label='dpi',default=300,min=50,
                                     max=600,nonzero=True,slider=True).set_pos(col=0)
                
            Params = SetParams()
            
            if Params.edit() == 1 and len(Params.Figures)>0:
                FigHandles = [FigHandles[k] for k in Params.Figures]
                Format  = ['.jpg','.png','.svg','.pdf'][Params.Format]
                SaveDir = Params.ImgDir
                dpi     = Params.dpi
                
                if not Params.FigName:
                    FigName='Fig_'
                else:
                    FigName=Params.FigName
                    
            else:
                return
        else:
            return

    for k in FigHandles:
        plt.figure(k.number)
        plt.savefig(os.path.join(SaveDir,FigName+'_'+str(k.number)+Format),
                    dpi = dpi)

########################################################################################################################

def Outliers4Skew(x):
    ''' Taken from:
    G. Brys; M. Hubert; A. Struyf (2004). A Robust Measure of Skewness.
    % Journal of Computational and Graphical Statistics 13(4), 996-1017.'''

    x_med = np.median(x)
    xi = x[x<=x_med]
    xj = x[x>=x_med]

    h=[]

    for i in xi:
        for j in xj:
            if (j-i)==0:
                h.append(0)
                continue
            else:
                h.append(((j-x_med)-(x_med-i))/(j-i))

    MedCouple = np.median(h)
    p         = np.percentile(x,[25,75])
    IQR       = p[1]-p[0]
    Lower     = p[0]-1.5*np.exp(-3.5*MedCouple)*IQR
    Upper     = p[1]+1.5*np.exp( 4.0*MedCouple)*IQR

    return MedCouple, Lower, Upper

########################################################################################################################

def bootci(Data, nSamples = 1000, Stat=['mean','median'][0], alpha = 5):
    '''Calculates the confidence interval by generating
    nSamples with replacement from a population'''

    # imports
    from scipy.stats import scoreatpercentile, nanmean

    # convert Data into an array
    Data  = np.array(Data)

    # get its length
    lSamp = Data.size
    
    # return if smaller than 5 samples
    if lSamp <= 5:
        return np.array([Data.mean(), Data.mean()])

    # generate n random indexes
    index = np.array([np.random.randint(0,lSamp,lSamp) for y in range(nSamples)])

    # calculate the statistic
    if Stat=='mean':
        mSamp = nanmean(Data[index], axis=1)
    elif Stat=='median':
        mSamp = np.median(Data[index], axis=1)

    # calculate the confidence interval
    CI    = np.array([scoreatpercentile(mSamp, alpha),
                      scoreatpercentile(mSamp, 100-alpha)])
    return CI
