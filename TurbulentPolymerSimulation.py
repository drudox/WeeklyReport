#!/usr/bin/env python3
# -*- coding : utf-8 -*-

'''
    Class for store and elaborate 
    Turbulent - Polymer - Simulation Data

    @author Marco Ghiani 9/2019

'''

#---------------------------------------------------------------

import sys
import numpy as np
import matplotlib.pyplot as plt 
import argparse 
from subprocess import call
from os.path import exists

try:
    from qualityPlot.qualityplot import *
except FileNotFoundError as exc:
    
    print("Error: \n qualityPlot.qualityplot FILE NON FOUND")
    print("Program Terminate with EXIT 0")
    sys.exit(0)

#-----------------------------------------------------------------


DEBUG = True




class TurbulentPolyData :
    '''
        class defined for compute the quantities that define a 
        Polymer chain in a Homogeneous Isotropic incompressible Turbulent Flow 

    '''

    def __init__(self,args):
        
        if not len(args):      # check if list of args is empty
            self.usage()            # if empty print help .. 
            sys.exit(0)
        else:
            self.getValue()         # get value from command line 
            self.chceckLocalFolders()
####   -----------------------------------   Domain quantities   ---------------------------------       

        self.DomainEdge = 24    
        self.node  = 256
        self.cell_size = self.DomainEdge/self.node
        
####   ------------------------------------    scaling factor    ---------------------------------
        self.unitl , self.unitt, self.unitv = self.scalingFactor() 
 
####   -------------------- Compute and catch >>>>> FLUID <<<<< variable -------------------------
        
        self.scaledViscosity = self.getViscosity()

        self.viscosity       = self.unscaleViscosity()

        self.time, self.K, self.eps, self.uPrime, \
        self.uPrimeMean, self.epsilon, self.TurbKineticEnergy = self.flowQuantities()
        
        self.timeLambda, self.TaylorLambda, self.ReynoldsLambda, self.TaylorNumMean , \
        self.ReyLambdaMean = self.flowField()
       
#---------------------------------------  Turbulence Parameters ----------------------------------------

        self.eta, self.tau, self.LET = self.turbulenceParameters()
#------------------------------------------------------------------------------------------------------- 

#--------------------- Compute and catch  >>>> POLYMER <<<< quantities  --------------------------------
        
        self.molMass = self.getPolMolMass()
        self.kuhnMolMass = self.getKuhnMolMass()
        self.Zimm = self.getZimm()
        self.timePolymer, self.polymerLength ,self.DeborahNum, self.WeissenbergNum = self.polymerQuantities()
        self.lchain_mx , self.undim_lchain_mx = self.getLenChainmx()
        self.chainSize = self.getChainSize()
        '''
            for i in self.chainSize:
            print( i/self.lchain_mx[0])
        '''
#dilute limit dimensional chain size
#------------------------------------------------------------------------------------------------------

        '''' Scaling times respect Large Eddy Turn over Time''' 
        for i in range(len(self.local_paths)):
            
            self.time[i]            /= self.LET[i]
            self.timePolymer[i]     /= self.LET[i]
            self.timeLambda[i]      /= self.LET[i]

        
#-----------------------------------------------------------------------------------------------------------------
######-----------------------------------------       METHODS        ---------------------------------------######
#-----------------------------------------------------------------------------------------------------------------
    
    def getValue(self):
        ''' 
            get value from command line and initialize 
            variable !
        '''
        global DEBUG        
        

        parser = argparse.ArgumentParser( description="transfer and/or process Simulation's \
                                                       files stored in local folder, if \'-u\' (--update) \
                                                       option is present and set to \'True\' the local folders \
                                                       will be update with file presents in Simulation's folder" )
        parser.add_argument( '-d', '--paths', type=argparse.FileType(), required=True,\
                             help='file contains list of local folders storing Simulations data' )
        parser.add_argument( '-s', '--sources', type=argparse.FileType(), required=True,\
                             help='file contains list of (remote) folders storing Simulations result' )
        parser.add_argument( '-u', '--update' , default=False ,type=bool ,\
                             help= 'update files in local folder with files contained in folders of Simulations') 


        self.args = parser.parse_args()
        


## ----------------- SETTING VARIABLE PASSED BY COMMAND LINE ------------------ ##      
        
        self.local_paths = []
        self.simulation  = []
        with open(self.args.paths.name) as f:
            for line in f:
                self.local_paths.append(line.split()[0])


        with open(self.args.paths.name) as f:
            for line in f:
                self.simulation.append(line.split()[1])
        
        if DEBUG:
            for i in range(len(self.local_paths)):
                 print(self.local_paths[i] , '    '  ,self.simulation[i])
   
        #-----------------------------------------------------------      
        #with open(self.args.paths.name) as f:
        #    for line in f:
        #        self.local_paths.append(line.split()[0])
        
        #if DEBUG:
        #    for i in range(len(self.local_paths)):
        #         print(self.local_paths[i])     



        #-----------------------------------------------------------
        self.sources_paths = []
        with open(self.args.sources.name) as f:
            for line in f:
                self.sources_paths.append(line.split()[0])
        
        if DEBUG:
            for i in range(len(self.sources_paths)):
                 print(self.sources_paths[i])
   
        #-----------------------------------------------------------      
        if(self.args.update):
             self.update()
        #-----------------------------------------------------------
        if( len(self.local_paths) != len(self.sources_paths)):
            print( 'Error !\nlocal_paths MUST BE TO SAME DIMENSION of sources')
            sys.exit()
        else:
            print('local_paths size equal to sources_path size! \n ... proceed')
            
        

#--------------------------------------------------------------------------------------------------------
    def usage(self):
        print('usage!')

#---------------------------------------------------------------------------------------------------------
    
    def KolmogorovScale(self, viscosity , epsilon):
        eta = ((viscosity**3)/epsilon)**0.25
        return eta

    def KolmogorovTime(self,viscosity,epsilon):
        tau = (viscosity/epsilon)**0.5
        return tau

    def LargeEddy(self,edge,uPmean):
        let = edge/2/uPmean
        return let

#--------------------------------------------------------------------------------

    def scalingFactor(self):
        '''
            catch the scaling factor (from memo.dat) for length and time
        '''
        unitt, unitl, unitv  = [],[],[]
        for path,i in zip(self.local_paths,range(len(self.local_paths))):
            with open(path + 'memo.dat') as f:
                 for line in f:
                     if line.strip().startswith('unitl'):
                        unitl.append(float(line.split('=')[1].split()[0]))
                        unitt.append(float(line.split('=')[2].split()[0]))
                        #print('unitt %10.6E\nunitl %10.6E\n' %(unitt , unitl))
        
            unitv.append(unitl[i]/unitt[i])

        return unitl,unitt,unitv

#---------------------------------------------------------------------------------------------------------
    def flowQuantities(self):
        ''' 
              catch the flow quantities (output from code) 
              compute the velocity fluctuation u' 
              and finally perform the mean value of the quantities
        '''
        time, K, eps, uPrime = [],[],[],[]
        uPrimeMean, epsilon, TurbKineticEnergy = [],[],[]


        for path,i in zip(self.local_paths,range(len(self.local_paths))):
            time.append(np.genfromtxt( path + 'snf_quads.dat', usecols=(0,)))
            K.append(np.genfromtxt( path + 'snf_quads.dat', usecols=(1,))) 
            eps.append(np.genfromtxt( path + 'snf_quads.dat', usecols=(2,)))
            uPrime.append(np.sqrt(2/3 * K[i]) * self.unitv[i])
        
        for i in range(len(self.local_paths)):
            uPrimeMean.append(np.sqrt(2/3 * np.mean(K[i])) * self.unitv[i])
            epsilon.append(np.mean(eps[i]))
            TurbKineticEnergy.append(np.mean(K[i])) 

        return time, K, eps, uPrime ,uPrimeMean, epsilon, TurbKineticEnergy



#---------------------------------------------------------------------------------------------------------
    def getViscosity(self):
        visc = []
        for path in self.local_paths:
            with open(path + 'memo.dat') as f:
                check = False
                for line in f:
                    if line.strip().startswith('SCALING UNITS'):
                        check = True
                        pass
                    elif check & line.strip().startswith('nfvisc'):
                        visc.append(float(line.split('=')[1]))
                        #print('Unscaled viscosity : %13.6f' %nfvisc)
        return visc


#---------------------------------------------------------------------------------------------------------

    def unscaleViscosity(self):  
        '''
            knowing that viscosity = L^2/time compute the 
            physic viscosity unscaled it using the scaling factors 
        '''
        v = []
        for i in range(len(self.scaledViscosity)):
            v.append(self.scaledViscosity[i] * self.unitl[i]*self.unitl[i] / self.unitt[i])
        return v

#---------------------------------------------------------------------------------------------------------


    def flowField(self):

        timeLambda, TaylorLambda, ReynoldsLambda = [],[],[] 
        TaylorNumMean, ReyLambdaMean = [],[]

        for path in self.local_paths:
            timeLambda.append(np.genfromtxt(path +'snf_taylorscal.dat', usecols=(0,)))       
            TaylorLambda.append(np.genfromtxt(path +'snf_taylorscal.dat', usecols=(1,)))
            ReynoldsLambda.append(np.genfromtxt(path +'snf_taylorscal.dat', usecols=(2,)))
 #  
        
        for i in range(len(self.local_paths)):
            TaylorNumMean.append(np.mean(TaylorLambda[i])*self.unitl[i])
            ReyLambdaMean.append(np.mean(ReynoldsLambda[i]))

       
        return timeLambda, TaylorLambda, ReynoldsLambda, TaylorNumMean, ReyLambdaMean

 #------------------------------------------------------------------------------------------------------ 

    def turbulenceParameters(self):        

        eta, tau ,LET = [] , [] , []

        for i in range(len(self.local_paths)):
            #### KOLMOGOROV SCALE
            eta.append( self.KolmogorovScale(self.scaledViscosity[i],self.epsilon[i]) * self.unitl[i])
            #print('kolmogorov scale: ',self.eta[i])
            
            #### KOLMOGOROV TIME 
            tau.append(self.KolmogorovTime(self.scaledViscosity[i],self.epsilon[i]) * self.unitt[i])
            #print('kolmogorov time: ',self.tau[i])
            
            #### LARGE EDDY TURN OVER 
            LET.append(self.LargeEddy(self.DomainEdge,self.uPrimeMean[i])) 
            #print('Large Eddy Turn Over: ',self.LET[i])
            #print('-'*58)
        
        return eta, tau ,LET


#------------------------------------------------------------------------------------------------------ 

    def polymerQuantities(self):
        '''

        '''
        timePolymer, polymerLength, DeborahNum, WeissenbergNum = [], [], [], [] 

        for path in self.local_paths:
            timePolymer.append(np.genfromtxt( path +'polength.dat', usecols=(0,)))
            polymerLength.append(np.genfromtxt( path +'polength.dat', usecols=(1,)))
        
        #self.molMass     = self.getPolMolMass()
        #for i in range(len(self.local_paths)):
        #    print(self.molMass[i])
        # longest relaxation time (Zimm time)
        #Zimm        = self.getZimm()
        # Weissenberg Number
        for i in range(len(self.local_paths)):   
            WeissenbergNum.append(self.Zimm[i]/self.tau[i])
        # Deborah Number
        for i in range(len(self.local_paths)):
            DeborahNum.append(self.uPrime[i] * self.Zimm[i] / self.eta[i])     # up*Zimm/kolmogorov_length 


        return timePolymer, polymerLength, DeborahNum, WeissenbergNum 






#---------------------------------------------------------------------------------------------------------
    def getLenChainmx(self):
        '''
            catch the maximum chain length from memo.dat
        '''
        lchain_mx , undim_lchain_mx = [] , []

        for path in self.local_paths:
            with open(path + 'memo.dat') as f:
                 for line in f:
                     if line.strip().startswith('dimensional maximum single chain length, lchaimx'):
                          lchain_mx.append(float(line.split('=')[1]))
                          #print(lchain_mx)
                          break
                     elif line.strip().startswith('dimensionless maximum single chain length, lchaimx'):
                          undim_lchain_mx.append(float(line.split('=')[1]))
                          #print(undim_lchain_mx) 
        return lchain_mx , undim_lchain_mx

#---------------------------------------------------------------------------------------------------------
    
    def getChainSize(self):
        '''
        
        '''
        chainSize = []
        for path in self.local_paths:
                with open(path + 'memo.dat') as f:
                    for line in f:
                        if line.strip().startswith('dilute limit dimensional chain size'):
                            chainSize.append(float(line.split('=')[1]))
        return chainSize


#---------------------------------------------------------------------------------------------------------
    def getZimm(self):
        '''
            catch the Zimm longest relaxation time 
        '''
        zimm = []
        for path in self.local_paths:
            with open(path + 'memo.dat') as f:
                 for line in f:
                     if line.strip().startswith('dilute dimensional Zimm longest'):
                         zimm.append(float(next(f).split()[0]))
                         #zimm = float(zimm[0])
        return zimm

#------------------------------------------------------------------------------------------------
    def getPolMolMass(self):
        '''
        
        '''
        molMass = []
        for path in self.local_paths:
                with open(path + 'memo.dat') as f:
                    for line in f:
                        if line.strip().startswith('polymer molar mass'):
                            molMass.append(float(line.split('=')[1]))
        return molMass 
#------------------------------------------------------------------------------------------------
    def getKuhnMolMass(self):
        '''
        
        '''
        kuhn = []
        for path in self.local_paths:
                with open(path + 'memo.dat') as f:
                    for line in f:
                        if line.strip().startswith('Kuhn monomer molar mass'):
                            kuhn.append(float(line.split('=')[1]))
        return kuhn 
    
#------------------------------------------------------------------------------------------------


    def DeborahStatistics(self):
        '''
            perform several averages of De Number and of Polimer scaled length 
            then compute a mean of this average in order to validate steady state of simulation
            and plot the correct mean value

        '''
        polLength40, polLength80, polLength120,polLength160 = [],[],[],[]
        timePol40, timePol80, timePol120, timePol160   = [],[],[],[]
        Deborah40, Deborah80, Deborah120, Deborah160   = [],[],[],[]
        time40, time80, time120, time160               = [],[],[],[]
        
        
        listDeborah = [0,1,3,4,6,7,8,9,10,11,12,13,14,15]
        for i in listDeborah:   


            ####-----------------
            #Deborah0         = self.DeborahNum[(self.time > 0)  &  (self.time < 40)]
            #time0            =       self.time[(self.time > 0)  &  (self.time < 40)]
        
            ######----------------
        
            polLength40.append(self.polymerLength[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 80)]) 
            timePol40.append(self.timePolymer[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 80)]) 

            Deborah40.append(self.DeborahNum[i][(self.time[i] > 40)  &  (self.time[i] < 80)])
            time40.append(self.time[i][(self.time[i] > 40)  &  (self.time[i] < 80)])
        
            ######---------------- 
        
            polLength80.append(self.polymerLength[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 120)]) 
            timePol80.append(self.timePolymer[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 120)]) 
   
            Deborah80.append(self.DeborahNum[i][(self.time[i] > 40)  &  (self.time[i] < 120)])
            time80.append(self.time[i][(self.time[i] > 40)  &  (self.time[i] < 120)])

            ######----------------
        
            polLength120.append(self.polymerLength[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 120)]) 
            timePol120.append(self.timePolymer[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 120)]) 
     
            Deborah120.append(self.DeborahNum[i][(self.time[i] > 40)  &  (self.time[i] <90)])
            time120.append(self.time[i][(self.time[i] > 40)  &  (self.time[i] <90)])

            ######----------------
  
            polLength160.append(self.polymerLength[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] <200)]) 
            timePol160.append(self.timePolymer[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] <200)])

            Deborah160.append(self.DeborahNum[i][(self.time[i] > 40)  &  (self.time[i] <200)])
            time160.append(self.time[i][(self.time[i] > 40)  &  (self.time[i] <200)])
        
            ######----------------
        
        mean_pol_40, mean_pol_80, mean_pol_120, mean_pol_160 = [],[],[],[] 
        mean_De_40, mean_De_80, mean_De_120, mean_De_160 = [],[],[],[]  
        mean_pol_tot, mean_De_tot = [],[] 
        
        for i in range(len(listDeborah)):   
            mean_pol_40.append(np.mean(polLength40[i]))
            mean_pol_80.append(np.mean(polLength80[i]))
            mean_pol_120.append(np.mean(polLength120[i]))
            mean_pol_160.append(np.mean(polLength160[i]))

            mean_De_40.append(np.mean(Deborah40[i]))
            mean_De_80.append(np.mean(Deborah80[i]))
            mean_De_120.append(np.mean(Deborah120[i]))
            mean_De_160.append(np.mean(Deborah160[i]))


            mean_pol_tot.append(np.mean([mean_pol_40[i], mean_pol_80[i] , mean_pol_120[i] , mean_pol_160[i]]))
            mean_De_tot.append(np.mean([mean_De_40[i], mean_De_80[i] , mean_De_120[i] , mean_De_160[i]])) 

        ######-----------------
       
        
        mean_pol_tot.sort()
        self.Deborah, self.DePolLength  =  mean_De_tot, mean_pol_tot                       # Main list 
        
        #plt.plot(mean_De_tot,mean_pol_tot)
        #plt.plot(self.Deborah, self.polLength) 
        #plt.show()
        
        for j,i in zip(listDeborah,range(len(listDeborah))):   
           with open(self.local_paths[j] + 'DeborahNumber.dat','a+') as f:
               f.write('%12.6f  %12.6f\n' %(mean_pol_tot[i],mean_De_tot[i]) )
 

 #------------------------------------------------------------------------------------------------


    def DeStatistics(self):
        '''
            perform several averages of De Number and of Polimer scaled length 
            then compute a mean of this average in order to validate steady state of simulation
            and plot the correct mean value

        '''
        polLength40, polLength80, polLength120,polLength160 = [],[],[],[]
        timePol40, timePol80, timePol120, timePol160   = [],[],[],[]
        Deborah40, Deborah80, Deborah120, Deborah160   = [],[],[],[]
        time40, time80, time120, time160               = [],[],[],[]
        
        
        listDeborah = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        for i in listDeborah:   


            ####-----------------
            #Deborah0         = self.DeborahNum[(self.time > 0)  &  (self.time < 40)]
            #time0            =       self.time[(self.time > 0)  &  (self.time < 40)]
        
            ######----------------
        
            polLength40.append(self.polymerLength[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 80)]) 
            timePol40.append(self.timePolymer[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 80)]) 

            Deborah40.append(self.DeborahNum[i][(self.time[i] > 40)  &  (self.time[i] < 80)])
            time40.append(self.time[i][(self.time[i] > 40)  &  (self.time[i] < 80)])
        
            ######---------------- 
        
            polLength80.append(self.polymerLength[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 120)]) 
            timePol80.append(self.timePolymer[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 120)]) 
   
            Deborah80.append(self.DeborahNum[i][(self.time[i] > 40)  &  (self.time[i] < 120)])
            time80.append(self.time[i][(self.time[i] > 40)  &  (self.time[i] < 120)])

            ######----------------
        
            polLength120.append(self.polymerLength[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 120)]) 
            timePol120.append(self.timePolymer[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] < 120)]) 
     
            Deborah120.append(self.DeborahNum[i][(self.time[i] > 40)  &  (self.time[i] <90)])
            time120.append(self.time[i][(self.time[i] > 40)  &  (self.time[i] <90)])

            ######----------------
  
            polLength160.append(self.polymerLength[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] <200)]) 
            timePol160.append(self.timePolymer[i][(self.timePolymer[i] > 40) & (self.timePolymer[i] <200)])

            Deborah160.append(self.DeborahNum[i][(self.time[i] > 40)  &  (self.time[i] <200)])
            time160.append(self.time[i][(self.time[i] > 40)  &  (self.time[i] <200)])
        
            ######----------------
        
        mean_pol_40, mean_pol_80, mean_pol_120, mean_pol_160 = [],[],[],[] 
        mean_De_40, mean_De_80, mean_De_120, mean_De_160 = [],[],[],[]  
        mean_pol_tot, mean_De_tot = [],[] 
        
        for i in range(len(listDeborah)):   
            mean_pol_40.append(np.mean(polLength40[i]))
            mean_pol_80.append(np.mean(polLength80[i]))
            mean_pol_120.append(np.mean(polLength120[i]))
            mean_pol_160.append(np.mean(polLength160[i]))

            mean_De_40.append(np.mean(Deborah40[i]))
            mean_De_80.append(np.mean(Deborah80[i]))
            mean_De_120.append(np.mean(Deborah120[i]))
            mean_De_160.append(np.mean(Deborah160[i]))


            mean_pol_tot.append(np.mean([mean_pol_40[i], mean_pol_80[i] , mean_pol_120[i] , mean_pol_160[i]]))
            mean_De_tot.append(np.mean([mean_De_40[i], mean_De_80[i] , mean_De_120[i] , mean_De_160[i]])) 

        ######-----------------
       
        
        #mean_pol_tot.sort()
        self.deborah, self.DePolLen  =  mean_De_tot, mean_pol_tot                       # Main list 
        
#------------------------------------------------------------------------------------
    
    def plotDeborahNumber(self,show = True):
        '''
            Prepare the graph "Polymer scaled Length VS Debora Number" using 
            a plot style created by the script's author.

            The graph are presented in both linear and log-log scale 

        '''
        x,y = self.Deborah, self.DePolLength

        fig,axs = SansItalic_paper(**{'lines.linewidth': 4,'fontsize':20} )(1,2,(15.0,5.0))
        axs[0].set(xlim=[180,2450]) 
        axs[0].set(ylim=[-0.02,0.53]) 
        axs[0].plot(x,y,label='')
        axs[0].plot(x,y,linestyle='', marker='o',markersize=12,color='C0',alpha=0.7)
        axs[0].set_title('Scaled Length VS Deborah Num.',y=1.02,color='black')
        axs[0].set(xlabel='Deborah Number')
        axs[0].set(ylabel='l/l$_0$')
    
        ##
        axs[1].set(ylim=[0.0045,0.55]) 
        axs[1].set_xscale("log", nonposx='clip')
        axs[1].set_yscale("log", nonposy='clip')
        axs[1].yaxis.set_ticks_position('both')
        axs[1].xaxis.set_ticks_position('both')
        axs[1].tick_params(direction='out', length=6)
        axs[1].plot(x,y,color='C0',label='')
        axs[1].plot(x,y,linestyle='', marker='o',markersize=12,color='C0',alpha=0.7)
        axs[1].grid(linestyle=(0, (8, 8)))
     #  axs[1].yaxis.set_major_locator(MultipleLocator(0.01))
        axs[1].get_xaxis().set_tick_params(direction='out')
        ticks = [200,300,400,500,600,700,800,900,1000,1300,1600,1900,2200,2400,2600]
        axs[1].set_xticks(ticks)

        #dic = { 200 : '2 10$^2$', 300 : '3 10$^2$' , 400 : '4 10$^2$', 500 : '' ,600 : '6 10$^2$', 
        #        700 : '' , 800: '' , 900 : '' ,1000 : '10$^3$', 1300: '1.3 10$^3$', 1600:'', 1900: ''}

        #dic = { 200 : '2e2', 300 : '3e2' , 400 : '4e2', 500 : '5e2' ,600 : '6e2', 
        #        700 : '' , 800: '' , 900 : '' ,1000 : '1e3', 1300: '1.3e3', 1600:'', 1900: '1.9e3',
        #        2200:'', 2400: '',2600:'2.6e3'}
        dic = { 200 : '2 10$^2$', 300 : '3 $10^2$' , 400 : '4 10$^2$', 500 : '5 10$^2$' ,600 : '', 
            700 : '' , 800: '8 10$^2$' , 900 : '' ,1000 : '10$^3$', 1300: '1.3 10$^3$', 1600:'', 1900: '1.9 10$^3$',
            2200:'', 2400: '',2600:'2.6 10$^3$'}
        labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
        axs[1].set_xticklabels(labels)
    

        #axs[1].xaxis.set_major_locator(MultipleLocator(1000))
        #axs[1].xaxis.set_minor_locator(MultipleLocator(100))
        plt.yticks([ 0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,0.004,0.005,0.006,0.007,0.008,0.009])

        axs[1].set_title('Scaled Length VS Deborah Num. (log scale)',y=1.02)
        axs[1].set(xlabel='Deborah Number')
        axs[1].set(ylabel='l/l$_0$')
    
        plt.tight_layout()
        plt.savefig('DeborahNumber.pdf',bbox_inches='tight' )  
        if show :
            plt.show()

 #---------------------------------------------------------------------------------------------------------

    def ReynoldsStatistics(self):
        '''
            perform several averages of Re Number and of Polimer scaled length 
            then compute a mean of this average in order to validate steady state of simulation
            and plot the correct mean value

        '''
       
        polLength50, polLength100, polLength150,polLength200 = [],[],[],[]
        timePol50, timePol100, timePol150, timePol200   = [],[],[],[]
        Reynolds50, Reynolds100, Reynolds150, Reynolds200   = [],[],[],[]
        time50, time100, time150, time200               = [],[],[],[]
        
        self.Reynolds, self.RePolLength  = [],[]                           # Main list 

        
        listReynolds = [4,5,6,7,8,9,11,12,14,15]
        for i in listReynolds:   
            

            ####-----------------
        
            polLength50.append(self.polymerLength[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 100)]) 
            timePol50.append(self.timePolymer[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 100)]) 

            Reynolds50.append(self.ReynoldsLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] < 100)])
            time50.append(self.timeLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] < 100)])
        
            ######---------------- 
        
            polLength100.append(self.polymerLength[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 150)]) 
            timePol100.append(self.timePolymer[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 150)]) 
   
            Reynolds100.append(self.ReynoldsLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] < 150)])
            time100.append(self.timeLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] < 150)])

            ######----------------
        
            polLength150.append(self.polymerLength[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 200)]) 
            timePol150.append(self.timePolymer[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 200)]) 
     
            Reynolds150.append(self.ReynoldsLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] <200)])
            time150.append(self.timeLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] <200)])

            ######----------------
  
            polLength200.append(self.polymerLength[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] <250)]) 
            time200.append(self.timePolymer[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] <250)])

            Reynolds200.append(self.ReynoldsLambda[i][(self.timeLambda[i] > 40)  &  (self.timeLambda[i] <200)])
            time200.append(self.timeLambda[i][(self.timeLambda[i] > 40)  &  (self.timeLambda[i] <200)])
        
            ######----------------
        
        mean_pol_50, mean_pol_100, mean_pol_150, mean_pol_200 = [],[],[],[] 
        mean_Re_50, mean_Re_100, mean_Re_150, mean_Re_200 = [],[],[],[]  
        mean_pol_tot, mean_Re_tot = [],[] 
        
        for i in range(len(listReynolds)):   
            mean_pol_50.append(np.mean(polLength50[i]))
            mean_pol_100.append(np.mean(polLength100[i]))
            mean_pol_150.append(np.mean(polLength150[i]))
            mean_pol_200.append(np.mean(polLength200[i]))

            mean_Re_50.append(np.mean(Reynolds50[i]))
            mean_Re_100.append(np.mean(Reynolds100[i]))
            mean_Re_150.append(np.mean(Reynolds150[i]))
            mean_Re_200.append(np.mean(Reynolds200[i]))


            mean_pol_tot.append(np.mean([mean_pol_50[i], mean_pol_100[i] , mean_pol_150[i]]))
            #, mean_pol_160[i]]))
            mean_Re_tot.append(np.mean([mean_Re_50[i], mean_Re_100[i] , mean_Re_150[i]]))
            #, mean_Re_200[i]])) 

        ######-----------------
       
        
        mean_pol_tot.sort(key=float)
        mean_Re_tot.sort(key=float)
        
        self.Reynolds, self.RePolLength  =  mean_Re_tot, mean_pol_tot                       # Main list 

        #plt.plot(self.Reynolds, self.polLength) 
        #plt.plot(self.Reynolds, self.polLength,'s') 
        #plt.show()
        
        for j,i in zip(listReynolds,range(len(listReynolds))):   
           with open(self.local_paths[j] + 'ReynoldsNumber.dat','a+') as f:
               f.write('%12.6f  %12.6f\n' %(mean_pol_tot[i],mean_Re_tot[i]) )
     
############################################################################################################

    def ReStatistics(self):
        '''
            perform several averages of Re Number and of Polimer scaled length 
            then compute a mean of this average in order to validate steady state of simulation
            and plot the correct mean value

        '''
       
        polLength50, polLength100, polLength150,polLength200 = [],[],[],[]
        timePol50, timePol100, timePol150, timePol200   = [],[],[],[]
        Reynolds50, Reynolds100, Reynolds150, Reynolds200   = [],[],[],[]
        time50, time100, time150, time200               = [],[],[],[]
        
        self.rey, self.rePolLen  = [],[]                           # Main list 

        
        listRey = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        for i in listRey:   
            

            ####-----------------
        
            polLength50.append(self.polymerLength[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 100)]) 
            timePol50.append(self.timePolymer[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 100)]) 

            Reynolds50.append(self.ReynoldsLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] < 100)])
            time50.append(self.timeLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] < 100)])
        
            ######---------------- 
        
            polLength100.append(self.polymerLength[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 150)]) 
            timePol100.append(self.timePolymer[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 150)]) 
   
            Reynolds100.append(self.ReynoldsLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] < 150)])
            time100.append(self.timeLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] < 150)])

            ######----------------
        
            polLength150.append(self.polymerLength[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 200)]) 
            timePol150.append(self.timePolymer[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] < 200)]) 
     
            Reynolds150.append(self.ReynoldsLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] <200)])
            time150.append(self.timeLambda[i][(self.timeLambda[i] > 50)  &  (self.timeLambda[i] <200)])

            ######----------------
  
            polLength200.append(self.polymerLength[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] <250)]) 
            time200.append(self.timePolymer[i][(self.timePolymer[i] > 50) & (self.timePolymer[i] <250)])

            Reynolds200.append(self.ReynoldsLambda[i][(self.timeLambda[i] > 40)  &  (self.timeLambda[i] <200)])
            time200.append(self.timeLambda[i][(self.timeLambda[i] > 40)  &  (self.timeLambda[i] <200)])
        
            ######----------------
        
        mean_pol_50, mean_pol_100, mean_pol_150, mean_pol_200 = [],[],[],[] 
        mean_Re_50, mean_Re_100, mean_Re_150, mean_Re_200 = [],[],[],[]  
        mean_pol_tot, mean_Re_tot = [],[] 
        
        for i in range(len(listRey)):   
            mean_pol_50.append(np.mean(polLength50[i]))
            mean_pol_100.append(np.mean(polLength100[i]))
            mean_pol_150.append(np.mean(polLength150[i]))
            mean_pol_200.append(np.mean(polLength200[i]))

            mean_Re_50.append(np.mean(Reynolds50[i]))
            mean_Re_100.append(np.mean(Reynolds100[i]))
            mean_Re_150.append(np.mean(Reynolds150[i]))
            mean_Re_200.append(np.mean(Reynolds200[i]))


            mean_pol_tot.append(np.mean([mean_pol_50[i], mean_pol_100[i] , mean_pol_150[i]]))
            #, mean_pol_160[i]]))
            mean_Re_tot.append(np.mean([mean_Re_50[i], mean_Re_100[i] , mean_Re_150[i]]))
            #, mean_Re_200[i]])) 

        ######-----------------
       
        
#        mean_pol_tot.sort(key=float)
#        mean_Re_tot.sort(key=float)
        self.rey, self.rePolLen  =  mean_Re_tot, mean_pol_tot                       # Main list 

    
#------------------------------------------------------------------------------------

    def plotReynoldsNumber(self,show=True):
        '''
            Prepare the graph "Polymer scaled Length VS Reynolds Number" using 
            a plot style created by the script's author.

            The graph are presented in both linear and log-log scale 

        '''

        x, y = self.Reynolds, self.RePolLength

        fig,axs = SansItalic_paper(**{'lines.linewidth': 6} )(1,2,(15.0,5.0))
        #plt.ylim([0,0.55]) 
        axs[0].plot(x,y,color='C1',linewidth=4)
        axs[0].plot(x,y,linestyle='', marker='o',markersize=12,color='C0',alpha=0.7)
        axs[0].set(ylim=[0.015,0.48])
        axs[0].set(xlim=[67,120])
        axs[0].set_title('Scaled Length VS Re$_\lambda$',y=1.02,color='black')
        axs[0].set(xlabel='Reynolds Number')
        axs[0].set(ylabel='l/l$_0$')
        #axs[0].legend()
    
        axs[1].set_xscale("log", nonposx='clip')
        axs[1].set_yscale("log", nonposy='clip')
        axs[1].yaxis.set_ticks_position('both')
        axs[1].xaxis.set_ticks_position('both')
        axs[1].tick_params(direction='out', length=6)
        axs[1].plot(x,y,color='C1',label='',linewidth=4)
        axs[1].plot(x,y,linestyle='', marker='o',markersize=12,color='C0',alpha=0.7)
        axs[1].grid(linestyle=(0, (7, 7)))
        axs[1].xaxis.set_major_locator(MultipleLocator(10))
        axs[1].yaxis.set_major_locator(MultipleLocator(0.01))
        axs[1].get_xaxis().set_tick_params(direction='out')
        plt.xticks([70,80,90,100,110,120])               #get only ticks we want
        plt.yticks([0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5])               #get only ticks we want
        axs[1].set_title('Scaled Length VS Re$_\lambda$ (log scale)',y=1.02)
        axs[1].set(xlabel='Reynolds Number')
        axs[1].set(ylabel='l/l$_0$')
        axs[1].set(xlim=[67,120])
        axs[1].set(ylim=[0.0255,0.5])
        #plt.legend()
        plt.tight_layout()
        plt.savefig('ReynoldsTaylorNumber.pdf',bbox_inches='tight' )  
        if show: 
            plt.show()

#--------------------------------------------------------------------------------------------
     

    def WeissenbergStatistics(self):
        '''
            perform several averages of We Number and of Polimer scaled length 
            then compute a mean of this average in order to validate steady state of simulation
            and plot the correct mean value

        '''
        self.WePolLength, self.Weissenberg = [],[]
        
        
        listWeissenberg = [0,1,2,4,7,8,9,10,11,14,15]
        for i in listWeissenberg:   
            self.WePolLength.append(np.mean(self.polymerLength[i]))
            self.Weissenberg.append(self.WeissenbergNum[i])
        
        #plt.plot(self.Weissenberg,self.WePolLength)
        #plt.show()
#--------------------------------------------------------------------------------------------
     

    def WeStatistics(self):
        '''
            perform several averages of We Number and of Polimer scaled length 
            then compute a mean of this average in order to validate steady state of simulation
            and plot the correct mean value

        '''
        self.WePolLen, self.We = [],[]
        
        
        listWeissenberg = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        for i in listWeissenberg:   
            self.WePolLen.append(np.mean(self.polymerLength[i]))
            self.We.append(self.WeissenbergNum[i])
        
        #plt.plot(self.Weissenberg,self.WePolLength)
        #plt.show()


   
#--------------------------------------------------------------------------------------------
    def plotWeissenbergNumber(self,show=True):
        '''
            Prepare the graph "Polymer scaled Length VS Weissenberg Number" using 
            a plot style created by the script's author.

            The graph are presented in both linear and log-log scale 

        '''

        x,y= self.Weissenberg, self.WePolLength

        #fig,axs = SansItalic(**{'lines.linewidth': 4} )(1,2,(18.0,6.0))
        fig,axs = SansItalic_paper(**{'lines.linewidth': 6} )(1,2,(15.0,5.0))
        axs[0].set(ylim=[0,0.55]) 
        axs[0].plot(x,y,color='C2',linewidth=4,label='')
        axs[0].plot(x,y,linestyle='', marker='o',markersize=12,color='C0',alpha=0.7)
        axs[0].set_title('Scaled Length VS Weissenberg Num.',y=1.02,color='black')
        axs[0].set(xlabel='Weissenberg Number')
        axs[0].set(ylabel='l/l$_0$')
        axs[0].set(xlim=[45,440])
        axs[0].set(ylim=[-0.01,0.5]) 
        axs[1].set(ylim=[0.004,0.5]) 
        axs[1].set(xlim=[52,448])
        axs[1].set_xscale("log", nonposx='clip')
        axs[1].set_yscale("log", nonposy='clip')
        axs[1].yaxis.set_ticks_position('both')
        axs[1].xaxis.set_ticks_position('both')
        axs[1].tick_params(direction='out', length=6)
        axs[1].plot(x,y,color='C2',label='',linewidth=4)
        axs[1].plot(x,y,linestyle='', marker='o',markersize=12,color='C0',alpha=0.7)
        axs[1].grid(linestyle=(0, (7, 7)))
        axs[1].xaxis.set_major_locator(MultipleLocator(100))
        axs[1].yaxis.set_major_locator(MultipleLocator(0.1))
        axs[1].get_xaxis().set_tick_params(direction='out')
        plt.yticks([0.004,0.005,0.006,0.007,0.008,0.009,0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,0.5,0.4,0.3,0.2])      
        axs[1].set_title('Scaled Length VS Weissenberg Num.(log scale)',y=1.02)
        axs[1].set(xlabel='Weissenberg Number')
        axs[1].set(ylabel='l/l$_0$')
        plt.tight_layout()
        plt.savefig('WeissenbergNumber.pdf',bbox_inches='tight' )  
        if show :
            plt.show()
    
 
#------------------------------------------------------------------------------------
    def chceckLocalFolders(self): 
        ''' 
            verify if the directory already exist:
            if not create it , else do nothing ..
        '''
        print('Verify if local Folders exist...\n')
        for dirs,source_path in zip(self.local_paths,self.sources_paths):
            if not exists(dirs):
                print(dirs + "  doesn't exist I will create it!\n")
                call('mkdir ' + dirs , shell=True)
                print('...copy inside the simulation''s file:')        
                call( 'scp ' + source_path + 'snf_taylorscal.dat ' + dirs , shell=True )
                call( 'scp ' + source_path + 'polength.dat ' + dirs       , shell=True )
                call( 'scp ' + source_path + 'snf_quads.dat ' + dirs      , shell=True )
                call( 'scp ' + source_path + 'memo.dat ' + dirs           , shell=True )
            
        print('.. Done')
        
#------------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------
    
    def update(self):
        print('... update files:\n')
        
        for path,source_path in zip(self.local_paths,self.sources_paths):
            call( 'scp ' + source_path + 'snf_taylorscal.dat ' + path , shell=True )
            call( 'scp ' + source_path + 'polength.dat ' + path       , shell=True )
            call( 'scp ' + source_path + 'snf_quads.dat ' + path      , shell=True )
            call( 'scp ' + source_path + 'memo.dat ' + path           , shell=True )

#---------------------------------------------------------------------------------------------



##############################################################################################

def main():
    

    result1 = TurbulentPolyData(sys.argv)
    
    result1.DeborahStatistics()
    result1.plotDeborahNumber()
    result1.ReynoldsStatistics()
    result1.plotReynoldsNumber()
    result1.WeissenbergStatistics()
    result1.plotWeissenbergNumber()
##########################################################################################

if __name__ == '__main__' :
    main()





    










