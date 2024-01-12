import pandas as pd
import numpy as np
from matplotlib import markers, axes, pyplot as plt
import os
pd.options.mode.chained_assignment = None

class makePlot():
    muw, munw, sigma = 0.821e-3, 0.838e-3, 57e-3
    lnetwork, absK, absK0 = 3.0035e-3, 2.7203e-12, 1.85e-12
    area = lnetwork**2
    lmean, rmean = 0.0001973, 2.2274e-05
    por = 0.2190

    def __init__(self, num, title,  results,
                 compWithLitData=False, compWithPrevData=False, drain=False, imbibe=False, exclude=None, include=None):

        self.colorlist = ['g', 'c', 'y', 'm', 'k', 'b', 'lightcoral', 'lime',    
                          'navy', 'tomato', 'khaki', 'olive', 'gold', 'teal', 'darkcyan', 'tan', 'limegreen']
        self.markerlist = ['v', '^', '<', '>', 'p', 'P','d', 'D', 'h', 'H', 's', 'o', 'v', '^', 
                           's', 'v', 'o', '^',  'd', 'D', 'h', 'H']
        self.linelist = ['--', ':', '-.', '--', (0, (1, 1)), (0, (5, 10)), (0, (5, 1)), 
                         (0, (3, 1, 1, 1)), (0, (3, 10, 1, 10, 1, 10)), (0, (3, 1, 1, 1, 1, 1))]
        self.num = num
        self.title = title
        self.compWithLitData = compWithLitData
        self.compWithPrevData = compWithPrevData
        self.drain = drain
        self.imbibe = imbibe
        self.exclude = exclude
        self.include = include
        self.results = results
        self.img_dir = "./result_images/"
        os.makedirs(os.path.dirname(self.img_dir), exist_ok=True)

        if self.drain:
            drainageBank(self)
        elif self.imbibe:
            imbibitionBank(self)
            
    def pcSw(self):
        if self.drain:
            filename = self.img_dir+'Pc_vs_Sw_Drainage_{}_{}.jpg'.format(
                self.title, self.num)
        else:
            filename = self.img_dir+'Pc_vs_Sw_Imbibition_{}_{}.jpg'.format(
                self.title, self.num)
        
        leg = []
        ind = 0        
                    
        for val in self.results.keys():
            res = self.results[val]
            if val == 'Literature data':
                res = res['pcSw']
                res1 = res.loc[res['source'] == 'MICP']
                res2 = res.loc[res['source'] != 'MICP']
                if not res1.empty:
                    plt.scatter(res1['satW'], res1['capPres']/1000, s=30, marker='o',
                                facecolors='none', edgecolors='b')
                    leg.append('MICP')
                if not res2.empty:
                    plt.scatter(res2['satW'], res2['capPres']/1000, s=30, marker='s',
                                facecolors='none', edgecolors='k')
                    leg.append(val)
            elif val == 'model':
                plt.plot(res['satW'], res['capPres']/1000, '-r', linewidth=2)
                leg.append(val)
            else:
                plt.plot(res['satW'], res['capPres']/1000, linestyle=self.linelist[ind], 
                         color=self.colorlist[ind], linewidth=2)
                leg.append(val)
                ind += 1
            
        plt.ylabel('Capillary Pressure(kPa)')
        plt.legend(leg)
        plt.ylim(0, 25)
        plt.xlim(0, 1)
        plt.xlabel('Sw')
        plt.savefig(filename, dpi=500)
        plt.close()

    def krSw(self):
        if self.drain:
            filename = self.img_dir+'kr_vs_Sw_Drainage_{}_{}.jpg'.format(
                self.title, self.num)
        else:
            filename = self.img_dir+'kr_vs_Sw_Imbibition_{}_{}.jpg'.format(
                self.title, self.num)
        
        leg = []
        j = 0
        for val in self.results.keys():
            res = self.results[val]
            if val == 'Literature data':
                res = res['krSw']
                print(res)
                plt.scatter(res['satW'], res['krw'], s=30, marker='s',
                    facecolors='none', edgecolors='k')
                leg.append('Literature data (krw)')
                plt.scatter(res['satW'], res['krnw'], s=30, marker='o',
                            facecolors='none', edgecolors='b')
                leg.append('Literature data (krnw)')
                
            elif val == 'model':
                plt.plot(res['satW'], res['krw'], linestyle='-',
                        color='r', linewidth=2)
                plt.plot(res['satW'], res['krnw'], linestyle='-',
                        color='r', linewidth=2, label = '_nolegend_')
                leg.append(val)
            
            else:
                plt.plot(res['satW'], res['krw'], linestyle=self.linelist[j], linewidth=2,
                        color=self.colorlist[j])
                plt.plot(res['satW'], res['krnw'], linestyle=self.linelist[j], linewidth=2,
                        color=self.colorlist[j], label = '_nolegend_')
                j += 1
                leg.append(val)
            
        plt.ylabel('Relative Permeability')
        plt.legend(labels=leg)
        plt.xlabel('Sw')
        plt.xlim(0, 1.0)
        plt.ylim(0, 1.0)
        plt.savefig(filename, dpi=500)
        plt.close()

    def krSw_c(self):
        filename = self.img_dir+'kr_vs_Sw_Probable_{}_{}.jpg'.format(
                self.title, self.num)
        leg = []
        resImb = self.results['imbibition']
        resProb = self.results['model_c']
        resProb1 = self.results['model_c_1']
        print(resProb1)
        
        plt.plot(resImb['satW'], resImb['krw'], linestyle='-', color='r', linewidth=2)
        leg.append('imbibition')
        plt.plot(resImb['satW'], resImb['krnw'], linestyle='-', color='r', linewidth=2,
                 label = '_nolegend_')
        
        #from IPython import embed; embed()
        #for ind, c in enumerate(np.unique(resProb['c'])):
        #'''
        for ind, c in enumerate([1,100,200,300,500,1000,1200]):
            resM = resProb.loc[(resProb['c'] == c)]
            resM.sort_values(by='fw', inplace=True)
            resM = self.removeAbnormalities2(resM)
            plt.plot(resM['psatW'], resM['pkrw'], '--v', color=self.colorlist[ind], 
                    markersize=2)
            leg.append('c={}'.format(c))
            plt.plot(resM['psatW'], resM['pkrnw'], '--v', color=self.colorlist[ind],       
                    markersize=2, label='_nolegend_')
            
        '''
        for ind, c in enumerate([5,100,200,300,500,1000,1100]):
            resM = resProb1[(resProb1['c'] == c)]
            resM.sort_values(by='fw', inplace=True)
            resM = self.removeAbnormalities2(resM)
            plt.plot(resM['psatW'], resM['pkrw'], '--*', color=self.colorlist[ind], 
                    markersize=4)
            leg.append('c={}'.format(c))
            plt.plot(resM['psatW'], resM['pkrnw'], '--*', color=self.colorlist[ind],
                    markersize=4, label='_nolegend_') '''
            
        plt.legend(leg)
        plt.xlabel('Sw')
        plt.xlim(0.0, 1.0)
        plt.ylim(0.0, 1.0)
        plt.ylabel('Relative permeability')
        plt.tight_layout()
        plt.savefig(filename, dpi=500)
        print('im done')
        plt.close()

    def pressgrad_c_fw_single(self):
        filename = self.img_dir+'gradP_vs_Ca_Probable_{}_{}.jpg'.format(
                        self.title, self.num)
        
        xlim0, xlim1 = -7, -4
        ylim0, ylim1 = 4, 7

        atol, beta = 2e-2, 1e-7
        resProb = self.results['model_c_fw']
        resProb['vt'] = (resProb['pqW']+resProb['pqNW'])/self.area
        nablaP0 = 1/self.lnetwork
        resProb['lambda_t'] = resProb['pkrw']/self.muw + resProb['pkrnw']/self.munw
        resProb['alpha'] = 1/(self.sigma*((1-resProb['pfw'])/self.munw+resProb['pfw']/self.muw))
        resProb['loggradP'] = np.log10(beta/resProb['alpha']*resProb['c']*nablaP0/resProb['vt'])
        resProb['logCa'] = np.log10(beta*resProb['c'])
        leg = []
        for ind, fw in enumerate([0.2, 0.4, 0.5, 0.6, 0.7, 0.8]):
            resM = resProb[np.isclose(resProb['pfw'], fw, atol=atol)]
            resM.sort_values(by='c', inplace=True)
            #resM = self.removeAbnormalities(resM)
            print(resM)
            plt.plot(resM['logCa'], resM['loggradP'], '--v', color=self.colorlist[ind], 
                    markersize=3)
            leg.append('fw={}'.format(fw))

        plt.legend(leg)
        #plt.xlim(xlim0, xlim1)
        #plt.ylim(ylim0, ylim1)
        plt.xlabel('Log10(Ca)')
        plt.ylabel('Log10(Pressure Gradient Pa/m)')
        plt.tight_layout()
        plt.savefig(filename, dpi=500)
        print('im done')
        plt.close()
   

    def pressgrad_c_fw(self):        
        xlim0, xlim1 = -7, -4
        ylim0, ylim1 = 4, 7
        atol = 2e-2
        resProb = self.results['model_c_fw']
        resProb['vt'] = (resProb['pqW']+resProb['pqNW'])/self.area
        
        #from IPython import embed; embed()
        resProb['lambda_t'] = self.absK*(resProb['pkrw']/self.muw + resProb['pkrnw']/self.munw)
        resProb['nablaP0'] = resProb['vt']/resProb['lambda_t']
        resProb['alpha'] = 1/(self.sigma*((1-resProb['pfw'])/self.munw+resProb['pfw']/self.muw))

        resYihuai = self.results['Yihuai']
        resG = self.results['Gao']
        beta=[1e-7,1e-7,1e-7,1e-7,1e-7,1e-7]
        for ind, fw in enumerate([0.2, 0.4, 0.5, 0.6, 0.7, 0.8]):
            filename = self.img_dir+'gradP_vs_Ca_Probable_{}_{}_{}.jpg'.format(
                self.title, self.num, ind+1)
            leg = []
            resM = resProb[np.isclose(resProb['pfw'], fw, atol=atol)]
            resM.sort_values(by='c', inplace=True)

            resM['loggradP'] = np.log10(beta[ind]/resM['alpha']*resM['c']*
                                            resM['nablaP0']/resM['vt'])
            resM['logCa'] = np.log10(beta[ind]*resM['c'])
            #from IPython import embed; embed()
            #print(resM)
            #resM = self.removeAbnormalities_f(resM)
            #print(resM)
            #from IPython import embed; embed()
            resY = resYihuai[resYihuai['fw']==fw]
            plt.plot(resM['logCa'], resM['loggradP'], '--v', color='r', 
                  markersize=5)
            #plt.scatter(resM['logCa'], resM['loggradP'], s=30, marker='v', color='r')
            leg.append('model, fw={}'.format(fw))

            plt.plot(resY['logCa'], resY['loggradP'],'--^', color='b', 
                       markersize=5)
            #plt.scatter(resY['logCa'], resY['loggradP'], s=30, marker='s', 
             #           facecolors='none', edgecolors='b')
            leg.append('Zhang et. al. fw={}'.format(fw))
            if fw == 0.5:
                #plt.scatter(resG['logCa'], resG['loggradP'], s=30, marker='s',
                 #           facecolors='none', edgecolors='g')
                plt.plot(resG['logCa'], resG['loggradP'],'-->', color='g',
                        markersize=5)
                leg.append('Gao et. al. fw={}'.format(fw))

            plt.legend(leg)
            #plt.xlim(xlim0, xlim1)
            #plt.ylim(ylim0, ylim1)
            plt.xlabel('Log10(Ca)')
            plt.ylabel('Log10(Pressure Gradient Pa/m)')
            plt.tight_layout()
            plt.savefig(filename, dpi=500)
            print('im done')
            plt.close()

    def pressgrad_c_fw_1(self):        
        xlim0, xlim1 = -7, -4
        ylim0, ylim1 = 4, 7
        atol = 2e-2
        resProb = self.results['model_c_fw']
        # from IPython import embed; embed()

        resProb['lambda_t'] = self.absK*(resProb['pkrw']/self.muw + resProb['pkrnw']/self.munw)
        resProb['nablaP0'] = resProb['vt']/resProb['lambda_t']
        resProb['alpha'] = 1/(self.sigma*((1-resProb['pfw'])/self.munw+resProb['pfw']/self.muw))
        
        resYihuai = self.results['Yihuai']
        resG = self.results['Gao']
        print(resProb)
        #from IPython import embed; embed() 
        for ind, fw in enumerate([0.2, 0.4, 0.5, 0.6, 0.7, 0.8]):
            filename = self.img_dir+'gradP_vs_Ca_Probable_{}_{}_{}.jpg'.format(
                self.title, self.num, ind+1)
            leg = []
            resM = resProb[np.isclose(resProb['pfw'], fw, atol=atol)]
            beta=resM['vt0'].iloc[0]
            alpha = resM['alpha'].mean()
            nablaP0 = resM['nablaP0'].mean()
            resM.sort_values(by='c', inplace=True)
            #resM['loggradP'] = np.log10(beta/resM['alpha']*resM['c']*
             #                           resProb['nablaP0']/resM['vt'])
            resM['loggradP'] = np.log10(beta/alpha*resM['c']*
                            nablaP0/resM['vt'])
            resM['logCa'] = np.log10(beta*resM['c'])
            #from IPython import embed; embed()
            #print(resM)
            #resM = self.removeAbnormalities_f(resM)
            #print(resM)
            #from IPython import embed; embed()
            resY = resYihuai[resYihuai['fw']==fw]
            plt.plot(resM['logCa'], resM['loggradP'], '--v', color='r', 
                markersize=5)
            #plt.scatter(resM['logCa'], resM['loggradP'], s=30, marker='v', color='r')
            leg.append('model, fw={}'.format(fw))
            plt.plot(resY['logCa'], resY['loggradP'],'--^', color='b', 
                    markersize=5)
            #plt.scatter(resY['logCa'], resY['loggradP'], s=30, marker='s', 
            #           facecolors='none', edgecolors='b')
            leg.append('Zhang et. al. fw={}'.format(fw))
            if fw == 0.5:
                #plt.scatter(resG['logCa'], resG['loggradP'], s=30, marker='s',
                #           facecolors='none', edgecolors='g')
                plt.plot(resG['logCa'], resG['loggradP'],'-->', color='g',
                        markersize=5)
                leg.append('Gao et. al. fw={}'.format(fw))
            plt.legend(leg)
            #plt.xlim(xlim0, xlim1)
            #plt.ylim(ylim0, ylim1)
            plt.xlabel('Log10(Ca)')
            plt.ylabel('Log10(Pressure Gradient Pa/m)')
            plt.tight_layout()
            plt.savefig(filename, dpi=500)
            print('im done')
            plt.close()

    def removeAbnormalities(self, data):
        data1 = self.removeAbnormalities_f(data)
        data2 = self.removeAbnormalities_b(data)
        if len(data1) >= len(data2):
            return data1
        else:
            return data2
            
    def removeAbnormalities_1(self, data):
        cond = True
        data['keep'] = True
        while cond:
            data['keep'].iloc[:-1] = (
                (data['lambda_t'].values[1:] >= data['lambda_t'].values[:-1]))
            cond = ~data['keep'].all()
            data = data[data['keep']]
        return data

    def removeAbnormalities_f(self, data):
        cond = True
        data['keep'] = True
        while cond:
            data['keep'].iloc[:-1] = (
                (data['loggradP'].values[1:] >= data['loggradP'].values[:-1]) &
                (data['lambda_t'].values[1:] >= data['lambda_t'].values[:-1]))
            cond = ~data['keep'].all()
            data = data[data['keep']]            
        return data
    
    def removeAbnormalities_b(self, data):
        cond = True
        data['keep'] = True
        while cond:
            data['keep'].iloc[1:] = (
                (data['loggradP'].values[:-1] <= data['loggradP'].values[1:]) &
                (data['lambda_t'].values[:-1] <= data['lambda_t'].values[1:]))
            cond = ~data['keep'].all()
            data = data[data['keep']]            
        return data

    def removeAbnormalities2(self, data):
        cond = True
        data['keep'] = True
        while cond:
            data['keep'].iloc[1:] = (data['psatW'].values[1:] >= data['psatW'].values[:-1]) & (
                               (data['pkrw'].values[1:] >= data['pkrw'].values[:-1]) &
                              (data['pkrnw'].values[1:] <= data['pkrnw'].values[:-1]))
            cond = ~data['keep'].all()
            data = data[data['keep']]            

        return data
    
    def removeAbnormalities1(self, data):
        cond = True
        data['keep'] = True
        while cond:
            data['keep'].iloc[1:] = (data['psatW'].values[1:] >= data['psatW'].values[:-1]) & (
                    (data['loggradP'].values[1:] >= data['loggradP'].values[:-1]))
            cond = ~data['keep'].all()
            data = data[data['keep']]

        return data


class drainageBank:
    def __init__(self, obj):
        self.obj = obj
        if self.compWithLitData:
            self.__compWithLitData__()
        if self.compWithPrevData:
            self.__compWithPrevData__()

    def __getattr__(self, name):
        return getattr(self.obj, name)

    def __compWithLitData__(self):
        self.results['Literature data'] = {}
        self.results['Literature data']['pcSw'] = pd.read_csv(
            './results_csv/Exp_Results_Bentheimer_Drainage_Pc_Sw.csv',
            names=['source', 'satW', 'Pc', 'capPres'], sep=',',
            skiprows=1, index_col=False)
        self.results['Literature data']['krSw'] = pd.read_csv(
            './results_csv/Exp_Results_Bentheimer_Drainage_kr_Sw.csv',
            names=['satW', 'krw', 'krnw'], sep=',',
            skiprows=1, index_col=False)
        
        
        self.results['Valvatne et al.'] = pd.read_csv(
            './results_csv/pnflow_Bentheimer_Drainage_010725.csv',
            names=['satW', 'capPres', 'krw', 'krnw', 'RI'], sep=',',
            skiprows=1, index_col=False)
        
    def __compWithPrevData__(self):
        if self.include:
            todo = list(self.include)
        else:
            todo = np.arange(1, self.num).tolist()
            if self.exclude:
                todo = np.setdiff1d(todo, self.exclude).tolist()

        while True:
            try:
                n = todo.pop(0)
                self.results['model_'+str(n)] = pd.read_csv(
                    "./results_csv/FlowmodelOOP_{}_Drainage_{}.csv".format(self.title, n),
                    names=['satW', 'qWout', 'krw', 'qNWout', 'krnw', 'capPres', 'invasions'],
                    sep=',', skiprows=18, index_col=False)
            except FileNotFoundError:
                pass
            except IndexError:
                break


class imbibitionBank():
    def __init__(self, obj):
        self.obj = obj
        if self.compWithLitData:
            self.__compWithLitData__()
        if self.compWithPrevData:
            self.__compWithPrevData__()

    def __getattr__(self, name):
        return getattr(self.obj, name)
        
    def __compWithLitData__(self):
        self.results['Literature data'] = {}
        self.results['Literature data']['pcSw'] = pd.read_csv(
            './results_csv/Exp_Results_Bentheimer_Imbibition_Pc_Sw.csv',
            names=['source', 'satW', 'Pc', 'capPres'], sep=',',
            skiprows=1, index_col=False)
        self.results['Literature data']['krSw'] = pd.read_csv(
            './results_csv/Exp_Results_Bentheimer_Imbibition_kr_Sw.csv',
            names=['satW', 'krw', 'krnw'], sep=',',
            skiprows=1, index_col=False)
        
        self.results['Valvatne et al.'] = pd.read_csv(
            './results_csv/pnflow_Bentheimer_Imbibition_010725.csv', names=[
                'satW', 'capPres', 'krw', 'krnw', 'RI'], sep=',', skiprows=1,
            index_col=False)
            
    def __compWithPrevData__(self):
        if self.include:
            todo = list(self.include)
        else:
            todo = np.arange(1, self.num).tolist()
            if self.exclude:
                todo = np.setdiff1d(todo, self.exclude).tolist()
        while True:
            try:
                n = todo.pop(0)
                self.results['model_'+str(n)] = pd.read_csv(
                    "./results_csv/FlowmodelOOP_{}_Imbibition_{}.csv".format(self.title, n),
                    names=['satW', 'qWout', 'krw', 'qNWout', 'krnw', 'capPres', 'invasions'],
                    sep=',', skiprows=18, index_col=False)
            except FileNotFoundError:
                pass
            except IndexError:
                break