#!/usr/bin/env python3
# -*- coding : utf-8 -*- 

import numpy as np
import matplotlib.pyplot as plt 
from qualityPlot.qualityplot import *
from subprocess import call
from os.path import exists
import argparse 
import sys
#from TurbulentPolymerSimulation 
            



class PlotReLambda:

    def __init__(self,args,update=False):
        '''

        '''
        if(len(args)):
            self.edge = 24              # EDGE OF COMPUTATIONAL DOMAIN [cm]
            self.getArgs() 
            self.createLocalFolders()
            self.createVectors()
            self.LET = self.computeLargeEddyTime()
            self.scaleTime()
            self.plotVector()      

        else:
            pass

#------------------------------------------------------------------------------------
    def getArgs(self):
        ''' 
            get file , parameter , options from command line 
            
        '''
        
        parser = argparse.ArgumentParser(description='update local files from simulations directory \
                                                      then scale the simulation\'s time and plot the \
                                                      scaled time VS polymer scaled length')
        parser.add_argument('-d' , '--directory' ,  type=argparse.FileType(), required=True,
                              help='file contains list of local folders storing Simulations data')
        parser.add_argument('-s' , '--sources' ,  type=argparse.FileType(), required=True,
                              help='file contains list of simulation directory including path and server sources of simulations')
        parser.add_argument('-u' , '--update' , type=bool , default=False , required=False ,
                              help='update local files with the simulations output file')
        parser.add_argument('-p', '--plot_show' , type=bool , default=False , required=False,
                              help='show or not the plot produced by this script , if \'False\' save the image without open it')

        self.args = parser.parse_args()
        
        ##--------- assignment 
        
        self.local,self.remote,self.name  = [],[],[]
        
        with open(self.args.directory.name) as f:
            for line in f:
                self.local.append(line.split()[0])
                self.name.append(line.split()[1])
        
        with open(self.args.sources.name) as f:
            for line in f:
                self.remote.append(line.split()[0])
        
        if len(self.local) != len(self.remote):
            print('local folder and remote folder must be to the same size')
        
        if self.args.plot_show:
            self.show = True
        else:
            self.show = False
            
        ##-----

        if self.args.update:
            self.updateFile()

#------------------------------------------------------------------------------------
    def createLocalFolders(self): 
        ''' 
            verify if the directory already exist:
            if not create it , else do nothing ..
        '''
        
        for dirs in self.local:
            if not exists(dirs):
                call('mkdir ' + dirs, shell=True)

       
#------------------------------------------------------------------------------------------------------------------------------

    def updateFile(self):
        '''
            update the locals file with the last update from server's simulations directory
            in this case for the file contains time e polymer scaled length

        '''
        for local,remote in zip(self.local,self.remote):
            call('scp ' + remote + 'snf_taylorscal.dat ' + local, shell=True)

#------------------------------------------------------------------------------------------------------------------------------

    def createVectors(self):
        '''

        '''
        self.time, self.reynolds = [],[]    

        for path in self.local:
            self.time.append(np.genfromtxt( path + 'snf_taylorscal.dat', usecols=(0,) ))
            self.reynolds.append(np.genfromtxt( path + 'snf_taylorscal.dat', usecols=(2,) ))


#------------------------------------------------------------------------------------------------------------------------------

    def computeLargeEddyTime(self):
        ''' 
                catch the scaling factor used by the code and compute u' 
        '''
        unitt, unitl, unitv  = [],[],[]
        for path,i in zip(self.local,range(len(self.local))):
          # try:  
            with open( path + 'memo.dat') as f:
              for line in f:
                if line.strip().startswith('unitl'):
                    unitl.append(float(line.split('=')[1].split()[0]))
                    unitt.append(float(line.split('=')[2].split()[0]))
            unitv.append(unitl[i]/unitt[i])
        
        K , uPrimeMean,let = [],[],[] 
        for path,i in zip(self.local,range(len(self.local))):
            K.append(np.genfromtxt( path + 'snf_quads.dat', usecols=(1,))) 
            uPrimeMean.append(np.sqrt(2/3 * np.mean(K[i])) * unitv[i])
            let.append(self.edge/2/uPrimeMean[i])

        return let
#------------------------------------------------------------------------------------------------------------------------------
    
    def scaleTime(self):
        '''
            scaled the simulations real time with the Large Eddy time (large eddy turn over time)

        '''
        for i in range(len(self.local)):
            self.time[i] /= self.LET[i]
            #print(self.LET[i])

#------------------------------------------------------------------------------------------------------------------------------

    def plotVector(self):
        '''

        '''
        
        fig,axs = SansItalic_II(**{'lines.linewidth': 3, 'fontsize':24} )(4,1,(26.0,36.0))
        
        for i in range(0,4):
            axs[i].plot(self.time[i],self.reynolds[i]         ,label=self.name[i],color='C3')   # ,linestyle=(0, (1, 1)))
            axs[i].plot(self.time[i+4],self.reynolds[i+4]     ,label=self.name[i+4],color='C2')
            axs[i].plot(self.time[i+2*4],self.reynolds[i+2*4] ,label=self.name[i+2*4],color='C1')
            axs[i].plot(self.time[i+3*4],self.reynolds[i+3*4] ,label=self.name[i+3*4],color='C0') #,linestyle=(0, (1, 1)))
            axs[i].set(xlabel = r'$t/\tau_{LET} $')
            axs[i].set(ylabel = r'$Re_{/\Lambda} $')
            axs[i].set_title('m=3e7 - Polymer length', y=1.05)
            #axs[i].set_xlim([0,400])
            axs[i].legend()
        plt.savefig('timeVSreynolds.pdf')
        if self.show:
           plt.show()
            
        








#------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    plot1 = PlotPolLength(sys.argv)
