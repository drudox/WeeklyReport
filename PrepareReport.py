#!/usr/bin/python3 
# -*- coding : utf-8 -*- 


from TurbulentPolymerSimulation import *
from PlotTimeVSLength import *
from PlotTimeVSReynolds import *
import sys
from subprocess import call


class Report: 
    '''
        Weekly Report 
        @author: Marco Ghiani 
        
        @brief : class that:
                - generate a LaTeX file 
                - create Plot of simulations
                - compile the LaTeX sources File
    '''
    def __init__(self,date,period):
        self.date = date
        self.period = period
        self.Result = TurbulentPolyData(sys.argv[1:])
        """ 
            instance of TurbulentPolyData class catch and compute all the necessary 
            quantities from a Fluid Flow - Polymer 
        """
        self.PlotTvsL = PlotPolLength(sys.argv[1:])
        '''
            define a class that produce a plot of scaled Polymer length versus t/T_{LET}
        '''
        self.PlotTvsRe = PlotReLambda(sys.argv[1:])
        '''
            define a class that produce a plot of Taylor reynolds Number versus t/T_{LET}
        '''
        
        self.Result.getValue()
        
        self.Result.DeborahStatistics()           
        self.Result.DeStatistics()

        self.Result.ReynoldsStatistics()           
        self.Result.ReStatistics()           
        
        self.Result.WeissenbergStatistics()           
        self.Result.WeStatistics()           
        
        self.Result.plotWeissenbergNumber(show=False)
        self.Result.plotDeborahNumber(show=False)
        self.Result.plotReynoldsNumber(show=False)
        
        #call('python3 Plot.py' , shell=True)
        #call('python3 Plot_Re.py' , shell=True)


#--------------------------------------------------------------------------
#
#--------------------------------------------------------------------------
    def prepareLaTeX(self,
                    target,
                    tableTitle,
                    ):
        '''
            this is the method that create the LaTeX source file 
            given 
            a list of target ...
            and prepare a table with the results that comes from the 
            object istantiated from the TurbulentPolyData() class 
        '''
        with open(self.date + '.tex','w') as f:
            f.write('\\documentclass[]{article}\n'+
                    '\\usepackage{amsmath}\n'+
                    '\\usepackage{verbatim}\n'+
                    '\\usepackage{geometry}\n'+
                    '\\usepackage{rotating}\n'+
                    '\\usepackage{subfigure}\n'+
                    '\\usepackage{hyperref}\n'+
                    '\\usepackage{graphicx}'
                    '\\geometry{\n'+
                    '    paper=a4paper, % Change to letterpaper for US letter\n'+
                    '    inner=1.5cm, % Inner margin\n'+
                    '    outer=2.8cm, % Outer margin\n'+
                    '    bindingoffset=.5cm, % Binding offset\n'+
                    '    top=3.5cm, % Top margin\n'+
                    '    bottom=1.5cm, % Bottom margin\n'+
                    '    }\n'+
                    '    \\usepackage{empheq}\n'+
                    '    \\usepackage{xcolor}\n'+
                    '    \\newcommand{\\boxedeq}[2]{\\begin{empheq}[box={\\fboxsep=6pt\\fbox}]{align}\\label{#1}#2\\end{empheq}}\n'+
                    '\\title{Weekly report \\textbf{23-04/30-04}}\n'+
                    '\\author{Marco Ghiani}\n'+
                    '                        '+
                    '%  START THE DOCUMENTS\n'+
                    '\\begin{document}\n'+              
                    '\\maketitle\n\n\n'+
                    '\\section{Targets of this week}\n'+
                    'During this week my aims are: \n'         
                    )
            f.write('\\begin{itemize}\n')
            
            for t in target:
                f.write('   \\item '+ t + '\n')

            f.write('\\end{itemize}\n')
            f.write('')
            
            f.write('\\section{Principal quantities plot}\n')
            

            f.write('\\begin{figure}\n')
            f.write('   \\centering\n')
            f.write('       \\includegraphics[width=1\\textwidth]{DeborahNumber.pdf}\n')
            f.write('   \\label{DeNum}\n')
            f.write('   \\caption{Deborah number VS Pol. Scaled length}\n')
            f.write('\\end{figure}\n')
 
            f.write('')
            f.write('')
            f.write('')
 
            
            f.write('\\begin{figure}\n')
            f.write('   \\centering\n')
            f.write('       \\includegraphics[width=1\\textwidth]{WeissenbergNumber.pdf}\n')
            f.write('   \\label{WeNum}\n')
            f.write('   \\caption{Weissenberg number VS Pol. Scaled length}\n')
            f.write('\\end{figure}\n')
            
            f.write('')
            f.write('')
            f.write('')
 
 
            f.write('\\begin{figure}\n')
            f.write('   \\centering\n')
            f.write('       \\includegraphics[width=1\\textwidth]{ReynoldsTaylorNumber.pdf}\n')
            f.write('   \\label{ReNum}\n')
            f.write('   \\caption{Reynolds Taylor number VS Pol length}\n')
            f.write('\\end{figure}\n')
            
            f.write('')
            f.write('')
            f.write('')
 

            f.write('\\begin{figure}[hb]\n')
            f.write('   \\begin{center}\n')
            f.write('       \\subfigure[][Deborah Number VS pol. scaled length]{\n')
            f.write('       \\includegraphics[width=1\\textwidth]{DeborahNumber.pdf}\n ')
            f.write('   \\label{fig:1}\n ')
            f.write('}\n ')
            
            f.write('       \\subfigure[][Weissenberg number VS Pol. Scaled length]{\n ')
            f.write('       \\includegraphics[width=1\\textwidth]{WeissenbergNumber.pdf}\n ')
            f.write('   \\label{fig:2}\n ')
            f.write('}\n ')

            f.write('       \\subfigure[][Reynolds Taylor number VS Pol Scaled length]{\n')
            f.write('       \\includegraphics[width=1\\textwidth]{ReynoldsTaylorNumber.pdf}\n')
            f.write('   \\label{fig:3}\n')
            f.write('}\n')
            f.write('   \\end{center}\n')
            f.write('   \\caption{Principal polymer-flow quantities}\n')
            f.write('   \\label{fig:quantities}\n')
            f.write('\\end{figure}\n')
            

            f.write('')
            f.write('')
            f.write('')
 
#  -----    TimeVSchainLength.pdf
            f.write('\\begin{figure}\n')
            f.write('   \\centering\n')
            f.write('       \\includegraphics[width=1\\textwidth]{timeVSlength.pdf}\n')
            f.write('   \\label{ReNum}\n')
            f.write('   \\caption{Time/$\\tau_{LET}$  VS Pol chain Scaled length}\n')
            f.write('\\end{figure}\n')
            
            f.write('')
            f.write('')
            f.write('')
 
#  -----    TimeVSReynolds.pdf
            f.write('\\begin{figure}\n')
            f.write('   \\centering\n')
            f.write('      \\includegraphics[width=1\\textwidth]{timeVSreynolds.pdf}\n')
            f.write('   \\label{ReNum}\n')
            f.write('   \\caption{Time/$\\tau_{LET}$  VS Reynolds Taylor Number}\n')
            f.write('\\end{figure}\n')
 


            #---------------------------------------------------------------------------------------------------------
            #  ----            TABLE OF SIMULATIONS RESULT (features of flow and pol chain)             ------
            #---------------------------------------------------------------------------------------------------------
            #f.write('\\newpage')
            f.write('')
            f.write('')
            f.write('')
 
            f.write('\\section{Table of quantities}\n')
            N = 8
            if len(self.Result.simulation) <= N-3:
                  f.write('\\begin{table}\n \\centering\n \\caption{' + 
                           tableTitle + '}\n'+ '\\label{table1}\n' +
                           '\\begin{tabular}{\n')
                  for i in range(len(self.Result.simulation)+1):    
                      f.write('c')
                  f.write('}\n')
                  f.write('\\hline\n'+ '\\textbf{Quantities}')
                  
                  for i in self.Result.simulation:
                      f.write('& '+ '\\textbf{'+i+'}')
                  f.write('\\\\ \n\\hline\n')
                  f.write('\nPol Mol Mass')
                  for i in self.Result.molMass:
                      f.write('& ')
                      f.write('%1.1e' % float(i))
                      f.write(' $g/mol$')
                  f.write('\\\\ \n')
                  
                  f.write('\nKuhn monomer mm')
                      
                  for i in self.Result.kuhnMolMass:
                      f.write(' &  ')
                      f.write('%6.4f' %float(i))
                  f.write('\\\\ \n\\hline\n')
                  
                  f.write('\nZimm relax. t $\\tau_0$') 
                  for i in self.Result.Zimm :
                      f.write(' &  ')
                      f.write(' %8.6f' %float(i))
                  f.write('\\\\ \n\\hline\n')
      
                  f.write('\nPol.L Max \\textcolor{red}{dimless}') 
                  for i in self.Result.undim_lchain_mx:
                      f.write(' &  ') 
                      f.write(' %8.1f' %float(i))
                      f.write(' [-]')
                  f.write('\\\\ \n')
      
                  f.write('\nPol.L Max \\textcolor{red}{dimens}') 
                  for i in self.Result.lchain_mx:
                      f.write(' &  ') 
                      f.write(' %1.6E' %float(i))
                      f.write(' $cm$')
                  f.write('\\\\ \n')
      
                  f.write('\nEq. chain sz ')
                  for i in range(len(self.Result.lchain_mx)):
                      f.write(' &  ')
                      eq_c =  float(self.Result.chainSize[i])/float(self.Result.lchain_mx[i])
                      f.write('%f' %eq_c )
                      f.write('$l/l_{\\text{Max}}$')
                  f.write('\\\\ \n\\hline\n')
                      
                  f.write('\\\\ \n')
                  f.write('$\\lambda$') 
                  for i in self.Result.TaylorNumMean:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$cm$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\eta_{k}$') 
                  for i in self.Result.eta:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$cm$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\tau_\\eta$') 
                  for i in self.Result.tau:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$s$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\tau_{LET}$') 
                  for i in self.Result.LET:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$s$')
                  f.write('\\\\ \n')
      
                  f.write('$\\Delta x$')
                  for i in range(len(self.Result.local_paths)):
                      f.write(' &  ')
                      f.write(' %9.7f' %self.Result.cell_size)
                      f.write('$cm$')
                  f.write('\\\\ \n')
      
                  f.write('$Re_\\lambda$') 
                  for i in self.Result.ReyLambdaMean:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n\\hline\n\\\\\n')
                  

                  f.write('Deborah Number') 
                  for i in self.Result.deborah:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('Reynolds$_\\lambda$ Number') 
                  for i in self.Result.rey:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('Weissenberg Number') 
                  for i in self.Result.We:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
                  f.write('\\end{tabular}')
                  f.write('\\end{table}')

################################################################################


            if len(self.Result.simulation) <= N : 
                  '''  
                        table with less or with 8 simulations
                  '''
                  f.write('\\begin{sidewaystable}\n \\centering\n \\caption{' + 
                           tableTitle + '}\n'+ '\\label{table1}\n' +
                           '\\begin{tabular}{')
                  for i in range(len(self.Result.simulation)+1):    
                      f.write('c')
                  f.write('}\n')
                  f.write('\\hline\n'+ '\\textbf{Quantities}')
                  
                  for i in self.Result.simulation:
                      f.write('& '+ '\\textbf{'+i+'}')
                  f.write('\\\\ \n\\hline\n')
                  f.write('\nPol Mol Mass')
                  for i in self.Result.molMass:
                      f.write('& ')
                      f.write('%1.1e' % float(i))
                      f.write(' $g/mol$')
                  f.write('\\\\ \n')
                  
                  f.write('\nKuhn monomer mm')
                      
                  for i in self.Result.kuhnMolMass:
                      f.write(' &  ')
                      f.write('%6.4f' %float(i))
                  f.write('\\\\ \n\\hline\n')
                  
                  f.write('\nZimm relax. t $\\tau_0$') 
                  for i in self.Result.Zimm :
                      f.write(' &  ')
                      f.write(' %8.6f' %float(i))
                  f.write('\\\\ \n\\hline\n')
      
                  f.write('\nPol.L Max \\textcolor{red}{dimless}') 
                  for i in self.Result.undim_lchain_mx:
                      f.write(' &  ') 
                      f.write(' %8.1f' %float(i))
                      f.write(' [-]')
                  f.write('\\\\ \n')
      
                  f.write('\nPol.L Max \\textcolor{red}{dimens}') 
                  for i in self.Result.lchain_mx:
                      f.write(' &  ') 
                      f.write(' %1.6E' %float(i))
                      f.write(' $cm$')
                  f.write('\\\\ \n')
      
                  f.write('\nEq. chain sz ')
                  for i in range(len(self.Result.lchain_mx)):
                      f.write(' &  ')
                      eq_c =  float(self.Result.chainSize[i])/float(self.Result.lchain_mx[i])
                      f.write('%f' %eq_c )
                      f.write('$l/l_{\\text{Max}}$')
                  f.write('\\\\ \n\\hline\n')
                      
                  f.write('\\\\ \n')
                  f.write('$\\lambda$') 
                  for i in self.Result.TaylorNumMean:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$cm$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\eta_{k}$') 
                  for i in self.Result.eta:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$cm$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\tau_\\eta$') 
                  for i in self.Result.tau:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$s$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\tau_{LET}$') 
                  for i in self.Result.LET:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$s$')
                  f.write('\\\\ \n')
      
                  f.write('$\\Delta x$')
                  for i in range(len(self.Result.local_paths)):
                      f.write(' &  ')
                      f.write(' %9.7f' %self.Result.cell_size)
                      f.write('$cm$')
                  f.write('\\\\ \n')
      
                  f.write('$Re_\\lambda$') 
                  for i in self.Result.ReyLambdaMean:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n\\hline\n\\\\\n')
                  

                  f.write('Deborah Number') 
                  for i in self.Result.deborah:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('Reynolds$_\\lambda$ Number') 
                  for i in self.Result.rey:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('Weissenberg Number') 
                  for i in self.Result.We:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
                  f.write('\\end{tabular}')
                  f.write('\\end{sidewaystable}')


#########################################################################    
                  

            elif len(self.Result.simulation) >= N and  len(self.Result.simulation) <= 2*N : 
                  '''
                        table with more than 8 simulations and less than 2x8 = 16 simulations
                  '''
                  f.write('\\begin{sidewaystable}\n \\centering\n \\caption{' + 
                           tableTitle + '}\n'+ '\\label{table1}\n' +
                           '\\begin{tabular}{')
                  for i in range(len(self.Result.simulation[:N])+1):    
                      f.write('c')
                  f.write('}\n')
                  f.write('\\hline\n'+ '\\textbf{Quantities}')
                  
                  for i in self.Result.simulation[:N]:
                      f.write('& '+ '\\textbf{'+i+'}')
                  f.write('\\\\ \n\\hline\n')
                  f.write('\nPol Mol Mass')
                  for i in self.Result.molMass[:N]:
                      f.write('& ')
                      f.write('%1.1e' % float(i))
                      f.write(' $g/mol$')
                  f.write('\\\\ \n')
                  
                  f.write('\nKuhn monomer mm')
                      
                  for i in self.Result.kuhnMolMass[:N]:
                      f.write(' &  ')
                      f.write('%6.4f' %float(i))
                  f.write('\\\\ \n\\hline\n')
                  
                  f.write('\nZimm relax. t $\\tau_0$') 
                  for i in self.Result.Zimm[:N] :
                      f.write(' &  ')
                      f.write(' %8.6f' %float(i))
                  f.write('\\\\ \n\\hline\n')
      
                  f.write('\nPol.L Max \\textcolor{red}{dimless}') 
                  for i in self.Result.undim_lchain_mx[:N]:
                      f.write(' &  ') 
                      f.write(' %8.1f' %float(i))
                      f.write(' [-]')
                  f.write('\\\\ \n')
      
                  f.write('\nPol.L Max \\textcolor{red}{dimens}') 
                  for i in self.Result.lchain_mx[:N]:
                      f.write(' &  ') 
                      f.write(' %1.6E' %float(i))
                      f.write(' $cm$')
                  f.write('\\\\ \n')
      
                  f.write('\nEq. chain sz ')
                  for i in range(len(self.Result.lchain_mx[:N])):
                      f.write(' &  ')
                      eq_c =  float(self.Result.chainSize[i])/float(self.Result.lchain_mx[i])
                      f.write('%f' %eq_c )
                      f.write('$l/l_{\\text{Max}}$')
                  f.write('\\\\ \n\\hline\n')
                      
                  f.write('\\\\ \n')
                  f.write('$\\lambda$') 
                  for i in self.Result.TaylorNumMean[:N]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$cm$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\eta_{k}$') 
                  for i in self.Result.eta[:N]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$cm$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\tau_\\eta$') 
                  for i in self.Result.tau[:N]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$s$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\tau_{LET}$') 
                  for i in self.Result.LET[:N]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$s$')
                  f.write('\\\\ \n')
      
                  f.write('$\\Delta x$')
                  for i in range(len(self.Result.local_paths[:N])):
                      f.write(' &  ')
                      f.write(' %9.7f' %self.Result.cell_size)
                      f.write('$cm$')
                  f.write('\\\\ \n')
      
                  f.write('$Re_\\lambda$') 
                  for i in self.Result.ReyLambdaMean[:N]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n\\hline\n\\\\\n')
                  

                  f.write('Deborah Number') 
                  for i in self.Result.deborah[:N]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('Reynolds$_\\lambda$ Number') 
                  for i in self.Result.rey[:N]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('Weissenberg Number') 
                  for i in self.Result.We[:N]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('\\end{tabular}')
                  f.write('\\end{sidewaystable}')
            
                #pass
                  f.write('\\begin{sidewaystable}\n \\centering\n \\caption{' + 
                           tableTitle + '}\n'+ '\\label{table1}\n' +
                           '\\begin{tabular}{')
                  for i in range(len(self.Result.simulation[N:])+1):    
                      f.write('c')
                  f.write('}\n')
                  f.write('\\hline\n'+ '\\textbf{Quantities}')
                  
                  for i in self.Result.simulation[N:]:
                      f.write('& '+ '\\textbf{'+i+'}')
                  f.write('\\\\ \n\\hline\n')
                  f.write('\nPol Mol Mass')
                  for i in self.Result.molMass[N:]:
                      f.write('& ')
                      f.write('%1.1e' % float(i))
                      f.write(' $g/mol$')
                  f.write('\\\\ \n')
                  
                  f.write('\nKuhn monomer mm')
                      
                  for i in self.Result.kuhnMolMass[N:]:
                      f.write(' &  ')
                      f.write('%6.4f' %float(i))
                  f.write('\\\\ \n\\hline\n')
                  
                  f.write('\nZimm relax. t $\\tau_0$') 
                  for i in self.Result.Zimm[N:] :
                      f.write(' &  ')
                      f.write(' %8.6f' %float(i))
                  f.write('\\\\ \n\\hline\n')
      
                  f.write('\nPol.L Max \\textcolor{red}{dimless}') 
                  for i in self.Result.undim_lchain_mx[N:]:
                      f.write(' &  ') 
                      f.write(' %8.1f' %float(i))
                      f.write(' [-]')
                  f.write('\\\\ \n')
      
                  f.write('\nPol.L Max \\textcolor{red}{dimens}') 
                  for i in self.Result.lchain_mx[N:]:
                      f.write(' &  ') 
                      f.write(' %1.6E' %float(i))
                      f.write(' $cm$')
                  f.write('\\\\ \n')
      
                  f.write('\nEq. chain sz ')
                  for i in range(len(self.Result.lchain_mx[N:])):
                      f.write(' &  ')
                      eq_c =  float(self.Result.chainSize[i])/float(self.Result.lchain_mx[i])
                      f.write('%f' %eq_c )
                      f.write('$l/l_{\\text{Max}}$')
                  f.write('\\\\ \n\\hline\n')
                      
                  f.write('\\\\ \n')
                  f.write('$\\lambda$') 
                  for i in self.Result.TaylorNumMean[N:]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$cm$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\eta_{k}$') 
                  for i in self.Result.eta[N:]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$cm$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\tau_\\eta$') 
                  for i in self.Result.tau[N:]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$s$')
                  f.write('\\\\ \n')
                  
                  f.write('$\\tau_{LET}$') 
                  for i in self.Result.LET[N:]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                      f.write('$s$')
                  f.write('\\\\ \n')
      
                  f.write('$\\Delta x$')
                  for i in range(len(self.Result.local_paths[N:])):
                      f.write(' &  ')
                      f.write(' %9.7f' %self.Result.cell_size)
                      f.write('$cm$')
                  f.write('\\\\ \n')
      
                  f.write('$Re_\\lambda$') 
                  for i in self.Result.ReyLambdaMean[N:]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n\\hline\n\\\\\n')
                  

                  f.write('Deborah Number') 
                  for i in self.Result.deborah[N:]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('Reynolds$_\\lambda$ Number') 
                  for i in self.Result.rey[N:]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('Weissenberg Number') 
                  for i in self.Result.We[N:]:
                      f.write(' &  ')
                      f.write(' %9.7f' %float(i))
                  f.write('\\\\ \n')
 
                  f.write('\\end{tabular}')
                  f.write('\\end{sidewaystable}')
            

            f.write('\n\\end{document}')
        

        call('pdflatex '+ self.date+'.tex' ,shell=True) 
        call('open '+ self.date+'.pdf',shell=True )


def main():
    
   doc = Report('5Aug','1/08/2019-5/08/2009')

   #doc.DeborahStatistics()
   doc.prepareLaTeX(['first','second','third'], 'Polymer molar mass = 3e7 g/mol')


if __name__ == '__main__':
    main()
