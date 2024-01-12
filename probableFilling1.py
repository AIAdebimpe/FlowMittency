import numpy as np
from time import time
from scipy.optimize import bisect, newton_krylov, root_scalar
import os
from utilities import Computations


class ProbableFilling:
    global ins
    def __init__(self, obj, createFile=False):
        global do
        self.obj = obj
        self.createFile = createFile
        do = Computations(self)
        self.__writeHeaders__()
        

    def __getattr__(self, name):
       return getattr(self.obj, name)
    
    def __computePcstar1__(self):
        global guess1, guess2, err
        N = 2000
        self._fillQuasi = 1-self.fluid        

        '''for c in [1e-6, 1e-3, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                   200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000,
                   4000, 5000, 6000, 7000, 8000, 9000, 10000]:'''
        self.randNum = self.rand([N, self.totElements])
        '''for c in [1e-6, 1e-3, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
                 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000]:
            #for c in [900, 1000]:'''
        cList =  (c for c in [1e-6, 1e-3, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 
                              400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000])
        for c in cList:
            self.c = c
            guess1, guess2, err = -1000, 1000, 1e6

            while guess2 <= 10000:
                #err = 1e6
                try:

                    #print('*******************************', c)
                    self.Pcstar = bisect(
                        self.funPcstar, guess2, guess1, args=(c, self.fw), xtol=1e-5)
                    self.prob = probfinal
                    self.pfw = self.computeFw()
                    self.Pcstar = finalPcstar
                    if (abs(self.pfw-self.fw) > 0.005):
                        print('check the convergence')
                        from IPython import embed; embed()

                    if (abs(self.pfw-self.fw) <= 0.02):
                        self.psatW = do.Saturation(self.probareaW, self.AreaSPhase)
                        self.writeProbResult(self.pfQ)
                    break
                except ValueError:
                    guess1 += 1000
                    guess2 += 1000

            if guess2 == 10000:
                print('ive got to 10000')
                from IPython import embed; embed()

    def __computePcstar__(self):
        global guess1, guess2, err
        N = 2000
        self._fillQuasi = 1-self.fluid
        self.randNum = self.rand([N, self.totElements])
        cList = [1e-6, 1e-3, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 
                    400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000]
        n = 41
        xList = np.linspace(*[0, 10000], n)

        for c in cList:
            err = 1e6
            s = 1
            #from IPython import embed; embed()
            y0 = self.funPcstar(xList[s-1], c, self.fw)
            y1 = self.funPcstar(xList[s], c, self.fw)
            while s < n-1:
                while (np.sign(y0) == np.sign(y1)) and (s < n-1):
                    y0 = y1
                    s += 1
                    y1 = self.funPcstar(xList[s], c, self.fw)
                print(y0, y1, xList[s], c, s)
                try:
                    self.Pcstar = bisect(self.funPcstar, xList[s-1], xList[s], 
                                    args=(c, self.fw), xtol=1e-5)
                    self.prob = self.computeProbability(self.Pcstar, c)
                    self.pfw = self.computeFw()
                    self.psatW = do.Saturation(self.probareaW, self.AreaSPhase)
                    self.c = c
                    self.writeProbResult(self.pfQ)
                    assert abs(self.pfw-self.fw) < 0.02
                    break
                except AssertionError:
                    y0 = y1
                    s += 1
                    y1 = self.funPcstar(xList[s], c, self.fw)
                except ValueError:
                    break
        
        del self.randNum


            
            
        
        
        '''for c in cList:
            err = 1e6
            y = [self.funPcstar(x, c, self.ins.fw) for x in xList]
            # Find where adjacent signs are not equal
            sign_changes = np.where(np.sign(y[:-1]) != np.sign(y[1:]))[0]
            ''''''root_finders = (root_scalar(f=self.funPcstar,
                                        args=(c, self.ins.fw),
                                        bracket=(xList[s], xList[s+1]))
                            for s in sign_changes)''''''
            roots = np.array([bisect(self.funPcstar, xList[s], xList[s+1], 
                            args=(c, self.ins.fw), xtol=1e-5) for s in sign_changes])
            #roots = np.array([r.root if r.converged else np.nan for r in root_finders])
            roots = np.unique(roots[~np.isnan(roots)])
            print(c, roots)
            try:
                err = np.array([self.funPcstar(x, c, self.ins.fw) for x in roots])
                self.Pcstar = roots[abs(err) == abs(err).min()]
                self.prob = self.computeProbability(self.Pcstar, c)
                self.pfw = self.computeFw()
                self.psatW = do.Saturation(self.probareaW, self.AreaSPhase)
                self.c = c
                self.writeProbResult(self.pfQ)
            except ValueError:
                pass'''
                
                
        #from IPython import embed; embed()
            
            
            
            
        #roots = multi_root(lambda x: self.funPcstar(x, c, self.ins.fw), [-5000, 10000])

        '''while guess2 <= 10000:
            #err = 1e6
            try:
                
                
                #print('*******************************', c)
                self.Pcstar = bisect(
                    self.funPcstar, guess2, guess1, args=(c, self.ins.fw), xtol=1e-5)
                self.prob = probfinal
                self.pfw = self.computeFw()
                self.Pcstar = finalPcstar
                if (abs(self.pfw-self.ins.fw) > 0.005):
                    print('check the convergence')
                    sol = newton_krylov(lambda x: self.funPcstar(x, c, self.ins.fw), 1000, 
                            line_search='wolfe')
                    #from IPython import embed; embed()
                    from IPython import embed; embed()
                if (abs(self.pfw-self.ins.fw) <= 0.02):
                    self.psatW = do.Saturation(self.probareaW, self.AreaSPhase)
                    self.c = c
                    self.writeProbResult(self.pfQ)
                break
            except ValueError:
                guess1 += 1000
                guess2 += 1000
        if guess2 == 10000:
            print('ive got to 10000')
            from IPython import embed; embed()'''
            
    
    def funPcstar(self, x, c, fw):
        global probfinal, finalPcstar, err
        self.prob = self.computeProbability(x, c)
        self.pfw = self.computeFw()
        if abs((self.pfw-fw)/fw) < err:
            probfinal = self.prob.copy()
            finalPcstar = x
            err = abs((self.pfw-fw)/fw)
        #print(x, c, (self.pfw-fw)/fw, self.prob.sum(), self.pfw, fw)
        return (self.pfw-fw)/fw


    def computeProbability(self, x, c):
        arrI = (self.PcI - x)/c
        arrD = (x - self.PcD)/c
        fillProb = 1.0/(1.0 + np.exp(arrD-arrI))
        prob = (fillProb > self.randNum).mean(axis=0)
        cond = (np.abs(prob-self._fillQuasi) > 0.5)
        prob[cond] = self._fillQuasi[cond]
        return prob


    def computeFw(self):
        cond = self.trappedNW & (self.prob == self._fillQuasi)

        self.probareaW = self.cornerArea+self.prob*self.centerArea
        probgw = np.zeros(self.totElements)
        probgw[1:-1] = self.cornerCond[1:-1]+self.prob[1:-1]*(
            self.centerArea[1:-1]/self.AreaSPhase[1:-1]*self.gwSPhase[1:-1])
        probgnw = (1-self.prob)*self.centerCond
        
        self.probareaW[cond] = self.AreaWPhase[cond]
        probgw[cond] = self.gWPhase[cond]
        probgnw[cond] = self.gNWPhase[cond]
        probgwL = do.computegL(probgw)
        probgnwL = do.computegL(probgnw)
        self.pqW = do.computeFlowrate(probgwL)
        self.pqNW = do.computeFlowrate(probgnwL)
        self.pkrw = self.pqW/self.qwSPhase
        self.pkrnw = self.pqNW/self.qnwSPhase
        return self.pqW/(self.pqW+self.pqNW)
    
    
    def writeProbResult(self, fQ):
        print('pSw: %7.6g  \tpqW: %8.6g  \tpkrw: %8.6g  \tpqNW: %8.6g  \tpkrnw:\
              %6.6g  \tPcstar: %8.6g \tpfw: %6.3g \tc: %6.3g \tqt0: %8.6g' % (
              self.psatW, self.pqW, self.pkrw, self.pqNW, self.pkrnw,
              self.Pcstar, self.pfw, self.c, (
                self.pqW+self.pqNW)/self.Area_,))
            
        fQ.write('\n%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g' % (
            self.psatW, self.pqW, self.pkrw, self.pqNW, self.pkrnw,
            self.Pcstar, self.pfw, self.capPresMin, self.c, (
                self.pqW+self.pqNW)/self.Area_, ))
        
    
    def __writeHeaders__(self):
        _num = 1
        prev = os.path.join(self.dirname, "Results/FlowmodelOOP_"+
                            self.title+"_Probable_"+str(_num)+".csv")
        while True:
            file_name = os.path.join(self.dirname, "Results/FlowmodelOOP_"+
                                 self.title+"_Probable_"+str(_num)+".csv")
            if os.path.isfile(file_name): _num += 1
            else:
                if self.createFile:
                    self.pfQ = open(file_name, 'a+')
                    self.pfQ.write('\n======================================================================')
                    ''''self.pfQ.write('\n'+'%'+'Fluid properties:\nsigma (mN/m)  \tmu_w (cP)  \tmu_nw (cP)')
                    self.pfQ.write('\n'+'%'+ '\t%.6g\t\t%.6g\t\t%.6g' % (
                        self.sigma, self.muw, self.munw, ))
                    self.pfQ.write('\n'+'%'+' calcBox: \t %.6g \t %.6g' % (
                        self.calcBox[0], self.calcBox[1], ))
                    self.pfQ.write('\n'+'%'+'Wettability:')
                    self.pfQ.write('\n'+'%'+'model \tmintheta \tmaxtheta \tdelta \teta \tdistmodel')
                    self.pfQ.write('\n'+'%'+'%.6g\t\t%.6g\t\t%.6g\t\t%.6g\t\t%.6g\t' % (
                        self.wettClass, round(self.minthetai*180/np.pi,3), round(self.maxthetai*180/np.pi,3), self.delta, self.eta ,) + str(self.distModel),)
                    self.pfQ.write('\nmintheta \tmaxtheta \tmean  \tstd')
                    self.pfQ.write('\n'+'%'+'%3.6g\t\t%3.6g\t\t%3.6g\t\t%3.6g' % (
                        round(self.contactAng.min()*180/np.pi,3), round(self.contactAng.max()*180/np.pi,3), round(self.contactAng.mean()*180/np.pi,3), round(self.contactAng.std()*180/np.pi,3)))
                    
                    self.pfQ.write('\nPorosity:  '+'%3.6g' % (self.porosity))
                    self.pfQ.write('\nMaximum pore connection:  '+'%3.6g' % (self.maxPoreCon))
                    self.pfQ.write('\nAverage pore-to-pore distance:  '+'%3.6g' % (self.avgP2Pdist))
                    self.pfQ.write('\nMean pore radius:  '+'%3.6g' % (self.Rarray[self.poreList].mean()))
                    self.pfQ.write('\nAbsolute permeability:  '+'%3.6g' % (self.absPerm))'''
                    
                    self.pfQ.write('\n======================================================================')
                    #self.pfQ.write("\n"+"%"+"Sw\t qW(m3/s)\t krw\t qNW(m3/s)\t krnw\t Pc\t Invasions")
                else: self.pfQ = open(prev, 'a+')
                break
            prev = file_name
                    



