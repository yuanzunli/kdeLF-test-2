# version: 2021_08_24_11_20

import math
import numpy as np
import time
import os.path
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import minimize
from scipy.integrate import dblquad
import kdeLF.kde_fortran_t

Mpc=3.08567758E+022 
c=299792.458
logfilename='log.txt'
from astropy.cosmology import FlatLambdaCDM
from multiprocessing import cpu_count
Number_of_threads=cpu_count()
Number_of_threads_used=min(Number_of_threads-1,20)


class KdeLF():
    """A Flexible Method for Estimating Luminosity Functions via Kernel Density Estimation"""
    # solid_angle: The solid angle subtended by the survey, in the unit of sr.
    
    
    Lmin,Lmax,z1,z2,xz,yL,red,lum,weight,ndata,f_pilot,hi=0,0,0,0,0,0,0,0,0,0,0,0   
    MLE_done,MLE_result=False,0
    
    
    def __init__(self, sample_file, zbin, f_lim, solid_angle, H0=71, Om0=0.27, Lbin=None,automatic=False,small_sample=False,
                 adaptive=False,pilot=None,unbounded_lum=True,epson=1e-9):        
        self.sample_file=sample_file
        self.zbin=zbin
        self.Omega=solid_angle
        self.H0=H0
        self.Om0=Om0
        self.Lbin=Lbin
        self.automatic=automatic
        self.adaptive=adaptive
        self.unbounded_lum=unbounded_lum      
        self.small_sample=small_sample
        KdeLF.z1=min(zbin)
        KdeLF.z2=max(zbin)        
        self.f_lim=f_lim
        self.cosmo = FlatLambdaCDM(H0=self.H0, Om0=self.Om0)
        self.h1l,self.h1h,self.h2l,self.h2h,self.bh=0,0,0,0,0
        self.set_beta_fixed=False
        self.fixed_beta=0.4
        self.Neff=0
        self.weighting=False
        self.absolute_magnitude=False
        self.first=True
        self.data_loaded=False        
        self.epson=epson
        self.pilot=pilot
        self.mcmc_fit_done=False 

    def __call__(self, x):   
        return self.f_lim(x)
            
    #************************************************************************  

    def load_sample(self):   

        kdeLF.kde_fortran_t.check()
        with open(self.sample_file, 'r') as f: 
            try:
                red,lum,weight = np.loadtxt(f, usecols=(0,1,2), unpack=True)
                self.weighting=True  
            except:
                red,lum = np.loadtxt(f, usecols=(0,1), unpack=True)            
                print('z & L data loaded')
                
                #self.weighting=True
                weight=np.zeros(len(red))+0.5
                
                
        if self.Lbin is None:
            Lmin,Lmax=np.min(lum),np.max(lum)
            b=Lmin % 1
            if b>0.5:
                Lmin=np.floor(Lmin)+0.5
            else:
                Lmin=np.floor(Lmin)
        
            b=Lmax % 1
            if b>0.5:
                Lmax=np.ceil(Lmax)
            else:
                Lmax=np.floor(Lmax)+0.5
        else:
            Lmin=min(Lbin)
            Lmax=max(Lbin)
        KdeLF.Lmin=Lmin
        KdeLF.Lmax=Lmax
        #************************************* 
        select=((red>KdeLF.z1) & (red<KdeLF.z2))
        KdeLF.red=red[select]
        KdeLF.lum=lum[select]
        KdeLF.ndata=len(KdeLF.red)        
        if np.max(KdeLF.lum)<0:
            self.absolute_magnitude=True
        
        kdeLF.kde_fortran_t.params.weighting, kdeLF.kde_fortran_t.params.absolute_magnitude = self.weighting, self.absolute_magnitude             
        if self.weighting is True:
            KdeLF.weight=weight[select]
            KdeLF.weight=1/KdeLF.weight
            #print('min(weight)=',min(KdeLF.weight))
            if min(KdeLF.weight)<1:
                print('please check the weight values')
                return
            self.Neff=np.sum(KdeLF.weight)        
            #print('self.Neff',self.Neff)                                   
            kdeLF.kde_fortran_t.params.weight,kdeLF.kde_fortran_t.params.nw=KdeLF.weight,self.Neff                    
        else:
            self.Neff=KdeLF.ndata       
        
        KdeLF.yL=KdeLF.lum-self.f_lim(KdeLF.red)
        KdeLF.xz=np.log( (KdeLF.red-KdeLF.z1)/(KdeLF.z2-KdeLF.red) )
        kdeLF.kde_fortran_t.params.xz,kdeLF.kde_fortran_t.params.ylum,kdeLF.kde_fortran_t.params.red=KdeLF.xz,KdeLF.yL,KdeLF.red
        # Pass values to arrays in the fortran type object 'kde' 
        input_params=np.array([self.H0,self.Om0,KdeLF.z1,KdeLF.z2,KdeLF.Lmin,KdeLF.Lmax])
        kdeLF.kde_fortran_t.f2py_value(KdeLF.ndata,input_params,Number_of_threads_used)         
        kdeLF.kde_fortran_t.params.small_sample,kdeLF.kde_fortran_t.params.adaptive = self.small_sample,self.adaptive 
       
        rato=kdeLF.kde_fortran_t.initialize(self.epson)
        rato=rato*100
        print('epson',self.epson,'rato','%.6f' %rato,'%')
        
        limx=np.linspace(np.min(KdeLF.red),KdeLF.z2,50)
        limy=self.f_lim(limx)             # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        kdeLF.kde_fortran_t.params.limx,kdeLF.kde_fortran_t.params.limy=limx,limy        
        
        
#        fig, ax = plt.subplots(1, 1)
#        plt.plot(KdeLF.xz,KdeLF.yL,'.',ms=1,color=(0.0,0.0,0.0),alpha=0.5) 
#        plt.plot(KdeLF.xz,-KdeLF.yL,'.',ms=1,color=(1.0,0.0,0.0),alpha=0.5)        
#        plt.show()
        
        kdeLF.kde_fortran_t.params.unbounded_lum=self.unbounded_lum
        self.data_loaded=True
        #kdeLF.kde_fortran_t.time_lim()
        
    #************************************************************************
    
    def K(self,x,y):
        return 1/(2*np.pi)*np.exp(-0.5*(x*x+y*y))
    
    def f_ref(self,x,y):
        global h1,h2    
        x_m=(x-KdeLF.xz)/h1
        y_m=(y-KdeLF.yL)/h2
        y_p=(y+KdeLF.yL)/h2
        if self.weighting is True:
            f_kde=(self.K(x_m,y_m)+self.K(x_m,y_p))*KdeLF.weight
            return np.sum(f_kde) / (h1*h2*self.Neff)
        else:
            f_kde=self.K(x_m,y_m)+self.K(x_m,y_p)
            return np.sum(f_kde) / (h1*h2*KdeLF.ndata)

    def set_f_pilot_value(self,h1_pilot,h2_pilot):
        global h1,h2
        h1=h1_pilot
        h2=h2_pilot
        KdeLF.f_pilot=np.empty(KdeLF.ndata)
        KdeLF.hi=np.empty(KdeLF.ndata)
        f = open('f_pilot.txt', "w")
        for i in range(KdeLF.ndata):
            xi=KdeLF.xz[i]
            yi=KdeLF.yL[i]
            KdeLF.f_pilot[i]=self.f_ref(xi,yi)
            print('%.6f' %xi,'%.6f' %yi,'%.8f' % KdeLF.f_pilot[i],file=f)
        f.close()
        if self.set_beta_fixed:
            KdeLF.hi=KdeLF.f_pilot**(-self.fixed_beta)       
        return    
    
    
    def f_ada(self,x,y):
        global h10,h20
        ha1=h10*KdeLF.hi
        ha2=h20*KdeLF.hi
        x_m=(x-KdeLF.xz)/ha1
        y_m=(y-KdeLF.yL)/ha2
        y_p=(y+KdeLF.yL)/ha2
        f_kde=( self.K(x_m,y_m)+self.K(x_m,y_p) ) / (ha1*ha2)
        return np.sum(f_kde) / KdeLF.ndata
        
 
#    def fhi1(pdelta1):
#        global h1,h2    
#        xi=xz.reshape(ndata,1)
#        yi=yL.reshape(ndata,1)    
#        x_m=(xi-xz)/h1
#        y_m=(yi-yL)/h2
#        y_p=(yi+yL)/h2    
#        temp=( K(x_m,y_m) + K(x_m,y_p) ) /(h1*h2)
#        fh_i=np.log( (temp.sum(1) - 1/(2*np.pi*h1*h2)) * 2/(2*ndata-1)/(red+delta1) )
#        return np.sum(fh_i)
#    
#    def lnlike1(para):
#        global xz,h1,h2
#        h1=para[0]
#        h2=para[1]
#        delta1=para[2]    
#        xz=np.log( (red-z1)/(z2-red) )
#        return fhi1(delta1) 
#        
#    def fhi2(xi,yi):
#        global xz,h1,h2
#        x_m=(xi-xz)/h1
#        y_m=(yi-yL)/h2
#        y_p=(yi+yL)/h2
#        temp=( K(x_m,y_m) + K(x_m,y_p) ) /(h1*h2)
#        return np.sum(temp) - 1/(2*np.pi*h1*h2)    
#        
#    def lnlike2(para):
#        global xz,h1,h2
#        h1=para[0]
#        h2=para[1]
#        delta1=para[2]    
#        xz=np.log( (red-z1)/(z2-red) + delta1 )         
#        f_hi=np.zeros(ndata)
#        for i in range(ndata):
#            f_hi[i]=np.log( fhi2(xz[i],yL[i]) * 2/(2*ndata-1)/(red[i]+delta1) )    
#        lnlik=np.sum(f_hi)    
#        print(lnlik)
#        return lnlik
#    
#    import kde
#    from scipy.integrate import dblquad

    def f1d(self,y):
        global h
        y_m=(y-KdeLF.yL)/h
        y_p=(y+KdeLF.yL)/h
        f_kde=( np.exp(-y_m**2/2) + np.exp(-y_p**2/2) )/np.sqrt(2*np.pi)
        return np.sum(f_kde) / (h*KdeLF.ndata)

    def f1da(self,y):
        global h_0
        ha=h_0*KdeLF.hi
        y_m=(y-KdeLF.yL)/ha
        y_p=(y+KdeLF.yL)/ha
        f_kde=( np.exp(-y_m**2/2) + np.exp(-y_p**2/2) )/np.sqrt(2*np.pi) / ha
        return np.sum(f_kde) / KdeLF.ndata
        

    def set_f_pilot_1d(self,h_pilot):
        global h
        h=h_pilot
        KdeLF.f_pilot=np.empty(KdeLF.ndata)
        KdeLF.hi=np.empty(KdeLF.ndata)        
        for i in range(KdeLF.ndata):
            yi=KdeLF.yL[i]
            KdeLF.f_pilot[i]=self.f1d(yi)  
    
  
    def prepare_KDE(self,theta):                     #$$$$$$$$$$$$$$$ small_sample not finished
        global h10,h20
        global h1,h2
        if self.data_loaded is False:
            self.load_sample()
        
        if self.adaptive is True:
            h10,h20,beta=theta
            if self.pilot is None:
                self.adaptive=False
                self.get_optimal_h()
                self.adaptive=True
                h1_pilot,h2_pilot=KdeLF.MLE_result
            else:    
                h1_pilot,h2_pilot=self.pilot            
            self.set_f_pilot_value(h1_pilot,h2_pilot)
            KdeLF.hi=KdeLF.f_pilot**(-beta)
            kdeLF.kde_fortran_t.params.hi, kdeLF.kde_fortran_t.params.h10, kdeLF.kde_fortran_t.params.h20, kdeLF.kde_fortran_t.params.beta = KdeLF.hi, h10, h20, beta            
        else:
            h1,h2=theta
            kdeLF.kde_fortran_t.params.h1, kdeLF.kde_fortran_t.params.h2 = h1, h2           
        return 
        
        
    def p(self,z,L):      
        y=L-self.f_lim(z)
        x=np.log( (z-KdeLF.z1)/(KdeLF.z2-z) )
        result=self.f_ref(x,y) * ( 1/(z-KdeLF.z1) + 1/(KdeLF.z2-z) )      
        return result
        
    def p_ada(self,z,L):
        y=L-self.f_lim(z)
        x=np.log( (z-KdeLF.z1)/(KdeLF.z2-z) )
        result=self.f_ada(x,y) * ( 1/(z-KdeLF.z1) + 1/(KdeLF.z2-z) )
        return result        
    
    def p1d(self,z,L):
        y=L-self.f_lim(z)
        result=self.f1d(y)/(KdeLF.z2-KdeLF.z1)
        return result

    def p1da(self,z,L):
        y=L-self.f_lim(z)
        result=self.f1da(y)/(KdeLF.z2-KdeLF.z1)
        return result                
        
    def phi_kde(self,z,L,theta):        
        global h1,h2
        global h10,h20
                
        dvdz=self.cosmo.differential_comoving_volume(z).value      
        if self.small_sample is True:
            if self.adaptive:
                result=self.p1da(z,L)*self.Neff/(dvdz*self.Omega) 
            else:
                result=self.p1d(z,L)*self.Neff/(dvdz*self.Omega)        
        else:        
            if self.adaptive:
                h10,h20,beta=theta
                result=self.p_ada(z,L)*self.Neff/(dvdz*self.Omega)
            else:
                h1,h2=theta                
                result=self.p(z,L)*self.Neff/(dvdz*self.Omega) 
        return result    


    def log10phi(self,z,L,theta):      #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$           
        dvdz=self.cosmo.differential_comoving_volume(z).value    
        if self.adaptive:
            h10,h20,beta=theta
            kdeLF.kde_fortran_t.params.h10,kdeLF.kde_fortran_t.params.h20,kdeLF.kde_fortran_t.params.beta=h10,h20,beta    
        else:
            h1,h2=theta
            kdeLF.kde_fortran_t.params.h1,kdeLF.kde_fortran_t.params.h2=h1,h2            
        result=kdeLF.kde_fortran_t.prob(z,L)*self.Neff/(dvdz*self.Omega)
        return np.log10(result) 
 
#    
#    def g(z):
#        return np.fmax(Lmin,f_lim(z))
    
#    def p_see():
#        z=1.0
#        L=np.arange(f_lim(z),Lmax,0.1)
#        for i in range(len(L)):
#            print(z,L[i],p(z,L[i]))
#        return
#    
#    
#    def pp(L,z):
#        x=np.log( (z-z1)/(z2-z) )
#        y=L-f_lim(z)    
#        return kdeLF.kde_fortran_t.f_ref_f2py(x,y) * ( 1/(z-z1) + 1/(z2-z) )
    
    ##################################################################################################################################
    def _2lnlike(self,para):
        h1,h2=para     
        #print(para)                                          
        result= -2 * kdeLF.kde_fortran_t.lnlike(h1,h2)                
        return result
        

    def _2lnlike_1d(self,para):
        result=-2*kdeLF.kde_fortran_t.lnlike_1d(para)
        return result
        

    def _2lnlike_1da(self,para):
        h_0,beta=para        
        result=-2*kdeLF.kde_fortran_t.lnlike_1da(h_0,beta)
        return result  
                     
        
    def _2lnlike_adaptive(self,para):
        global h10,h20
        #print('para=',para)    
        if self.set_beta_fixed:
            h10,h20=para
            beta=self.fixed_beta
        else:
            h10,h20,beta=para    
        result= -2 * kdeLF.kde_fortran_t.lnlike_ada(h10,h20,beta)
        return result 
        

    def lnprior(self,theta):                       
        if self.adaptive:        
            if self.set_beta_fixed:
                h10, h20 = theta
                if self.h1l < h10 < self.h1h and self.h2l < h20 < self.h2h:
                    return 0.0                
            else:
                h10, h20, beta = theta
                if self.h1l < h10 < self.h1h and self.h2l < h20 < self.h2h and 0.0 < beta < self.bh:
                    return 0.0                
        else:        
            h1, h2 = theta
            if self.h1l < h1 < self.h1h and self.h2l < h2 < self.h2h:
                return 0.0
        return -np.inf
        
    
    def lnprior1(self,theta):
        if self.adaptive:        
            if self.set_beta_fixed:
                h10, h20 = theta
                result=-np.log(1+h10*h10) - np.log(1+h20*h20)                
                if h10>0 and h20>0:
                    return result              
            else:
                h10, h20, beta = theta
                result=-np.log(1+h10*h10) - np.log(1+h20*h20)
                if h10>0 and h20>0 and 0.0 < beta < self.bh:
                    return result                
        else:            
            h1, h2 = theta
            result=-np.log(1+h1*h1) - np.log(1+h2*h2)
            if h1>0 and h2>0:
                return result        
        return -np.inf             
    
    
    def lnprob(self,theta):
        if self.uniform_priors is True:
            lp = self.lnprior(theta)
        else:
            lp = self.lnprior1(theta)   
        if not np.isfinite(lp):
            return -np.inf
        if self.adaptive:
            result = lp - self._2lnlike_adaptive(theta)/2  
        else:
            result = lp - self._2lnlike(theta)/2         
        return result 
    

    def time_test(self):
        n=min(Number_of_threads,16)        
        Xn=np.zeros(n)
        Tn=np.zeros(n)
        x0=[ [0.15,0.15],[0.3,0.2] ]
        for j in range(0,n):
            kdeLF.kde_fortran_t.params.process=j+1        
            time_start=time.time()        
            for i in range(2):        
                yy=self.lnlike(x0[i])
            time_end=time.time()
            Xn[j]=j+1
            Tn[j]=time_end-time_start
            print(Xn[j],Tn[j])
            if j>0:
                if Tn[j]>Tn[j-1]:
                    N_optimal=j
        
        plt.plot(Xn,Tn,'ko--')
        plt.ylabel(r'Total time (s)',fontsize=13)
        plt.xlabel(r'Number of threads')
        #plt.show()                     
        return N_optimal
        
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    def get_h_1d(self,bound_h):
        global h
        from scipy.optimize import minimize_scalar
        KdeLF.xz=KdeLF.red
        kdeLF.kde_fortran_t.params.xz=KdeLF.red
        print('bandwidth for 1d estimator,')                  
        print('    bound for h:  ', bound_h)
        res = minimize_scalar(self._2lnlike_1d, bounds=bound_h, method='bounded')          
        h=res.x
        print('    Optimal h:    ', '%.4f' % h)              
        return np.array([h]) 
   
    def get_h_1da(self,h_pilot,bounds_h0_beta):
        global h_0
        self.set_f_pilot_1d(h_pilot)
        kdeLF.kde_fortran_t.params.hi, kdeLF.kde_fortran_t.params.f_pilot, kdeLF.kde_fortran_t.params.set_beta_fixed=KdeLF.hi, KdeLF.f_pilot,self.set_beta_fixed 
        x0=np.array([h,self.fixed_beta])        
        print('global bandwidth and beta for adaptive 1d estimator,')
        print('    Initial h0 & beta:     ',x0)        
        print('    bounds for h0 & beta:  ',bounds_h0_beta)
        res = minimize(self._2lnlike_1da, x0, method='Powell', bounds=bounds_h0_beta,options={'xtol': 1e-4, 'ftol': 1e-4, 'disp': False})       
        h_0,beta=res.x
        KdeLF.hi=KdeLF.f_pilot**(-beta)
        print('    Optimal h0 & beta:')
        print('        ','%.4f' % h_0, '%.4f' % beta)
        return np.array([h_0,beta]) 
    
    def get_h_2d(self,x0,bounds):
        global h1,h2           

        print('bandwidths for 2d estimator,')
        print('    Initial h1 & h2:     ',x0)
        print('    bounds for h1 & h2:  ',bounds)
        res = minimize(self._2lnlike, x0, method='Powell', bounds=bounds,options={'xtol': 1e-4, 'ftol': 1e-4, 'disp': False})            
        h1,h2=res.x
        print('    Optimal h1 & h2:')
        print('        ','%.4f' % h1, '%.4f' % h2)
        return np.array([h1,h2])
    
   
    def get_h_2da(self,h1_pilot,h2_pilot,x0,bounds_h0,bound_beta):
        global h10,h20
        self.set_f_pilot_value(h1_pilot,h2_pilot)
        kdeLF.kde_fortran_t.params.hi, kdeLF.kde_fortran_t.params.f_pilot, kdeLF.kde_fortran_t.params.set_beta_fixed=KdeLF.hi, KdeLF.f_pilot,self.set_beta_fixed            

        print('global bandwidths and beta for adaptive 2d estimator,')
        bounds=bounds_h0[:2]        
        if self.set_beta_fixed is False:            
            bounds.append(bound_beta)                
            x0=np.append(x0,self.fixed_beta)            
            print('    Initial h10, h20 & beta:     ',x0)
            print('    bounds for h10, h20 & beta:  ',bounds)
        else:            
            print('beta is fixed at value',self.fixed_beta)
            print('Initial bandwidths for adaptive:',x0)
            print('bounds for h10,h20:',bounds)            
        res = minimize(self._2lnlike_adaptive, x0, method='Powell', bounds=bounds,options={'xtol': 1e-4, 'ftol': 1e-4, 'disp': False})            
        if self.set_beta_fixed:
            h10,h20=res.x
            beta=self.fixed_beta
            print('    Optimal h10, h20 & beta (fixed):')
            
        else:
            h10,h20,beta=res.x
            KdeLF.hi=KdeLF.f_pilot**(-beta)
            #KdeLF.hi=(1+KdeLF.f_pilot)**(-beta)
            #KdeLF.hi=np.exp(-beta*KdeLF.f_pilot)
            print('    Optimal h10, h20 & beta:')
        print('        ','%.4f' % h10, '%.4f' % h20, '%.4f' % beta,'\n')
        
        return np.array([h10,h20,beta])    
    
    
    def get_optimal_h(self,initial_bandwidths=[0.15,0.15],bandwidth_bound=[(0.001, 1.0),(0.001, 1.0)],beta_bound=(0.01,0.6),set_beta_fixed=False,beta_value=0.3,quiet=False):
        """
        parameter minimize_method: scipy.optimize.minimize 'Powell'.                                  
         
        """        

        self.set_beta_fixed, self.fixed_beta = set_beta_fixed, beta_value         
        if self.data_loaded is False:
            self.load_sample()        

        if not quiet:
            print("Maximum likelihood estimation by scipy.optimize 'Powell' method,")
            print('redshift bin:','(',KdeLF.z1,',',KdeLF.z2,')')
            print('sample size of this bin:',KdeLF.ndata)
        x0=initial_bandwidths
        bounds=bandwidth_bound[:2]
        bound=[]
        bound.append(bandwidth_bound[1])
        time_start=time.time() 
        if self.small_sample is True:         
            KdeLF.MLE_result=self.get_h_1d(bounds[1])
            h_pilot=KdeLF.MLE_result[0]
            if self.adaptive is True:            
                print('pilot bandwidth:')
                print('    hp:', '%.4f' % h_pilot,'\n')                
                bound.append(beta_bound)
                KdeLF.MLE_result=self.get_h_1da(h_pilot,bound)
            KdeLF.MLE_done=True
            return KdeLF.MLE_result 
        
        #if small_sample is False, the ^phi or ^phi_a will work:                
  
        if self.adaptive is False:
            KdeLF.MLE_result=self.get_h_2d(x0,bounds)           
        else:  # adaptive is True
            if self.pilot is None:
                KdeLF.MLE_result=self.get_h_2d(x0,bounds)
                h1_pilot,h2_pilot=KdeLF.MLE_result[0],KdeLF.MLE_result[1]
            else:
                h1_pilot,h2_pilot=self.pilot[0],self.pilot[1]           
            print('pilot bandwidths:')                 
            print('    h1p & h2p:', '%.4f' % h1_pilot, '%.4f' % h2_pilot,'\n') 
            KdeLF.MLE_result=self.get_h_2da(h1_pilot,h2_pilot,x0,bounds,beta_bound)
        KdeLF.MLE_done=True     
        time_end=time.time()
        print("Cost total time:",time_end-time_start)
        return KdeLF.MLE_result
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f(self,z,L2):
        return self.f_lim(z)-L2
    
    def rho(self,L,z):
        dvdz=self.cosmo.differential_comoving_volume(z).value        
        return dvdz*self.Omega
    
    def get_binLF(self):        
        if self.data_loaded is False:
            self.load_sample()
        z1=KdeLF.z1
        z2=KdeLF.z2
        Lmin=KdeLF.Lmin 
        Lmax=KdeLF.Lmax
        Lbin=0.3
        nLbin=int((Lmax-Lmin)/Lbin)        
        L=np.empty(nLbin)
        phi=np.empty(nLbin)
        zz1=max(1e-6,z1)
        L0=max(Lmin,self.f_lim(zz1))
        f_lim_z2=self.f_lim(z2)
        L1=L0
        L2=L1+Lbin        
        k=0
        for i in range(nLbin):
            num=0
            for j in range(KdeLF.ndata): 
                #if KdeLF.red[j]>=z1 and KdeLF.red[j]<z2 and KdeLF.lum[j]>=L1 and KdeLF.lum[j]<L2:
                if KdeLF.lum[j]>=L1 and KdeLF.lum[j]<L2:
                    num=num+1
            #print('L1,L2,num',L1,L2,num)
            Lo=self.f_lim((z1+z2)/2)
            Lc=L1+Lbin/2
            if num>0 and Lc>=Lo:
                if L2<f_lim_z2:
                    x0=optimize.brentq(self.f,zz1,z2,args=(L2))
                    #print('x0',x0)
                else:
                    x0=z2
                area = dblquad(self.rho, z1, x0, lambda z: max(self.f_lim(z),L1), lambda z: L2)
                #print('area',area)
                phi_page=num/area[0]
                L[k]=Lc
                phi[k]=np.log10(phi_page)
                #print(k,L[k],phi[k])
                k=k+1
            L1=L1+Lbin
            L2=L2+Lbin
            if L1>=Lmax:
                break            
        result=np.empty((2,k)) 
        result[0]=L[0:k];result[1]=phi[0:k]             #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    
        return result
    #############################################################################
            

    def cdf(self,x):
        global z2s,clf
        return clf    

    def ECDF(self,data):
    # Compute ECDF 
        x1 = np.sort(data)        
        x=np.unique(data)
        if self.absolute_magnitude is True:
            x1=x1[::-1]
            x=x[::-1]        
        n = x1.size
        m=x.size
        y=np.zeros(m)
        num=0
        s=0
        for i in range(m):
            ni=0
            for j in range(s,n):
                if x[i]==x1[j]:
                    ni+=1
            s+=ni
            num=num+ni
            y[i]=num/n
        return np.array([x,y])
        
    
    def KStest(self,CDF_redshift=False):        
        from scipy import stats
        global z2s,clf     
      
        if self.mcmc_fit_done is True:
            self.prepare_KDE(self.mcmc_fit)
                
        lum=np.unique(KdeLF.lum) 
        if self.absolute_magnitude is True:
            lum=lum[::-1]        
        
        red=np.sort(KdeLF.red)
        n=len(lum)
        clf=np.zeros(n)
        czf=np.zeros(n)
        z2s=np.zeros(n)
        f_lim_z2=self.f_lim(KdeLF.z2)
        
        z1=max(1e-6,KdeLF.z1)
        for i in range(n):
            L2=lum[i]
            if abs(L2)<abs(f_lim_z2):
                z2s[i]=optimize.brentq(self.f,z1,KdeLF.z2,args=(L2))
                #print('z2s',z2s[i])
            else:
                z2s[i]=KdeLF.z2
        
        #time_start=time.time()
        #clf=kdeLF.kde_fortran_t.clf(lum,z2s)
        clf=kdeLF.kde_fortran_t.clf_new(lum,z2s)
        #time_end=time.time()
        #print("Cost total time CDF_luminosity:",time_end-time_start)      
        
        ecdf=self.ECDF(KdeLF.lum)        
        dl=abs(clf-ecdf[1])
        dl=max(dl)
      
        if CDF_redshift is True:
            czf=kdeLF.kde_fortran_t.czf(red)
            dz=abs(czf-ecdf)
            dz=max(dz)
                      
        #ks=stats.ks_1samp(lum, self.cdf)
        #print('d',dl,dz)
        #print(ks)
        plt.figure(figsize=(8,6)) 
        ax=plt.axes([0.1,0.1, 0.85, 0.85])
        ax.tick_params(direction='in', top=True, right=True, labelsize=12)
        
        for i in range(n-1):
            if i==0:
                plt.plot([ecdf[0][i], ecdf[0][i+1]], [ecdf[1][i], ecdf[1][i]], 'k', linewidth=0.8,label='empirical',alpha=0.7)
            else:
                plt.plot([ecdf[0][i], ecdf[0][i+1]], [ecdf[1][i], ecdf[1][i]], 'k', linewidth=0.8,alpha=0.7)            
            plt.plot([ecdf[0][i+1], ecdf[0][i+1]], [ecdf[1][i], ecdf[1][i+1]], 'k', linewidth=0.8,alpha=0.7)        
        
        if self.small_sample is False:
            if self.adaptive is False:        
                plt.plot(lum,clf,'r-',linewidth=1.5, label=r'$\mathrm{expected~by~}\hat{\phi}$',alpha=0.7)
            else:
                plt.plot(lum,clf,'r-',linewidth=1.5, label=r'$\mathrm{expected~by~}\hat{\phi}_{\mathrm{a}}$',alpha=0.7)       
        else:
            if self.adaptive is False:        
                plt.plot(lum,clf,'r-',linewidth=1.5, label=r'$\mathrm{expected~by~}\hat{\phi}_{\mathrm{1}}$',alpha=0.7)
            else:
                plt.plot(lum,clf,'r-',linewidth=1.5, label=r'$\mathrm{expected~by~}\hat{\phi}_{\mathrm{1a}}$',alpha=0.7)         
        plt.legend(fontsize=12)
        if self.absolute_magnitude is True:              
            plt.xlim(lum[0]*0.99, lum[n-1]*1.01)
            plt.xlabel(r'$M$',fontsize=18)
        else:
            plt.xlabel(r'$\log_{10} ~ L$',fontsize=18)    
        plt.ylabel('CDFs',fontsize=18)        
        #plt.show()
        plt.savefig('kstest.png')
        return dl          #ks.statistic,ks.pvalue
    
       
         
    def get_lgLF(self,z=None, plot_binLF=False, extrapolation=True):                 
        """        
        This function get the logarithm (log10) of the luminosity function (LF) estimated via kernel density estimation.        . 
        Input parameter z: The redshift at which the LF is estimated. By default, z is the center of the 
                           redshift bin (z1,z2), i.e., z=(z1+z2)/2.
        return: 2d array includes the luminosity and log LF vectors.        
        """
        if z is None:
            z=np.mean(KdeLF.red)        
        if (z<=KdeLF.z1) or (z>=KdeLF.z2):
            print('INPUT ERROR: The input redshift should satisfy',str(KdeLF.z1),"< z <",str(KdeLF.z2))
            return        
        if self.absolute_magnitude is True:
            LL=np.linspace(self.f_lim(z),KdeLF.Lmin,100)
        else:
            LL=np.linspace(self.f_lim(z),KdeLF.Lmax,100)
        
        rho=np.zeros(len(LL))
        for i in range(len(LL)):
            rho[i]=self.phi_kde(z,LL[i],KdeLF.MLE_result)                
        result=np.array([LL,np.log10(rho)])         
        
        lmax=max(KdeLF.lum)
        lmin=min(KdeLF.lum)
        LFmin=np.log10(rho[0])         
        print('lmax',lmax,lmin)
        
        if self.absolute_magnitude is True:
            LFmax=np.log10(self.phi_kde(z,lmin,KdeLF.MLE_result))
        else:
            LFmax=np.log10(self.phi_kde(z,lmax,KdeLF.MLE_result))     
        y1=LFmax-0.25*(LFmin-LFmax)
        y2=LFmin+0.1*(LFmin-LFmax)
        y1=max(y1,np.log10(rho[len(rho)-1]))
        print("log10 LF is given at z=",z)
        if extrapolation:
            print("The vertical dot dash line marks the high luminosity limit, L=",lmax,',','beyond which the LF is extrapolated.')        

        plt.figure(figsize=(8,6)) 
        ax=plt.axes([0.13,0.1, 0.82, 0.85])
        ax.tick_params(direction='in', top=True, right=True, labelsize=12) 

        if plot_binLF is True: 
            binLF=self.get_binLF()   
            plt.plot(binLF[0],binLF[1],'o',mfc='white',mec='black',ms=9,mew=1.0,alpha=0.7,label=r'$\hat{\phi}_{\mathrm{bin}}$') # plot_bin
        if self.small_sample is False:
            if self.adaptive is False:        
                plt.plot(result[0],result[1],color=(1.0, 0.0, 0.0),linewidth=2.0, label=r'$\hat{\phi}$')
            else:
                plt.plot(result[0],result[1],color=(1.0, 0.0, 0.0),linewidth=2.0, label=r'$\hat{\phi}_{\mathrm{a}}$')        
        else:
            if self.adaptive is False:        
                plt.plot(result[0],result[1],color=(1.0, 0.0, 0.0),linewidth=2.0, label=r'$\hat{\phi}_{\mathrm{1}}$')
            else:
                plt.plot(result[0],result[1],color=(1.0, 0.0, 0.0),linewidth=2.0, label=r'$\hat{\phi}_{\mathrm{1a}}$')                
        
        if self.absolute_magnitude is True:
            tx=lmax-(lmax-KdeLF.Lmin)*0.618
            ty=y1+(y2-y1)*0.88
            ax.text(tx,ty,'z='+'%.3f' %z,fontsize=13,bbox=dict(boxstyle='square,pad=0.3', fc='yellow', ec='k',lw=0.5 ,alpha=0.4)) 
            plt.plot([lmin,lmin],[-50,-0.1],linewidth=0.6, linestyle=(0, (5, 5, 1, 5)))
            plt.ylabel(r'$\log_{10}( \phi(z,M) ~/~ {\rm Mpc}^{-3} ~ \mathrm{mag}^{-1} )$',fontsize=18)
            plt.xlabel(r'$M$',fontsize=18)            
            plt.xlim(lmax*0.96,KdeLF.Lmin)
            
        else:
            tx=lmin+(KdeLF.Lmax-lmin)*0.618
            ty=y1+(y2-y1)*0.88
            ax.text(tx,ty,'z='+'%.3f' %z,fontsize=13,bbox=dict(boxstyle='square,pad=0.3', fc='yellow', ec='k',lw=0.5 ,alpha=0.4))
            plt.plot([lmax,lmax],[-20,-0.1],linewidth=0.6, linestyle=(0, (5, 5, 1, 5)))
            plt.ylabel(r'$\log_{10}( \phi(z,L) ~/~ {\rm Mpc}^{-3} ~ \Delta L^{-1} )$',fontsize=18)
            plt.xlabel(r'$\log_{10} ~ L$',fontsize=18)
            plt.xlim(lmin, KdeLF.Lmax)
        plt.ylim(y1,y2)
        plt.legend(loc='lower left', fontsize=12)
        if self.adaptive is True:
            figname='z'+'%.3f' %z+'a.png'
        else:
            figname='z'+'%.3f' %z+'.png'
        plt.savefig(figname)
        plt.show()
        return result


    def get_lgLF_given_h(self,bandwidth,z=None,global_bandwidth=None,beta_value=None):
        global h1,h2
        global h10,h20
        if self.data_loaded is False:
            self.load_sample()
        if self.adaptive is False:        
            h1,h2=bandwidth            
            KdeLF.MLE_result=np.array([h1,h2])
        else:    # self.adaptive is True
            if global_bandwidth is None:
                print("Please set values for the argument 'global_bandwidth'! This should be a list or array with 2 elements")
                return            
            h10,h20=global_bandwidth
            h1,h2=bandwidth
            h1_pilot,h2_pilot=bandwidth           
            self.set_f_pilot_value(h1_pilot,h2_pilot)
            if self.set_beta_fixed:
                beta=self.fixed_beta                
            else:
                if beta_value is None:
                    print("Please set a value for the argument 'beta_value'!")
                    return
                beta=beta_value
                KdeLF.hi=KdeLF.f_pilot**(-beta)
                
            KdeLF.MLE_result=np.array([h10,h20,beta])    
                                    
        if z is None:
            result=self.get_lgLF()
        else:
            z0=z
            result=self.get_lgLF(z=z0)        
        return result

    
    def run_mcmc(self, max_n=3000,Ntau=50, initial_point=None, parallel=False,uniform_priors=True):
        #import os.path
        #import os
        import emcee
        
        self.uniform_priors=uniform_priors               #$$$$$$$$$$$$$$$$$$$$$$$$        
        # Set up the backend
        # Don't forget to clear it in case the file already exists
        if self.adaptive is True:
            chain_name = "chain_a" + str(KdeLF.z1) + "_" + str(KdeLF.z2) + ".h5"
        else:
            chain_name = "chain_" + str(KdeLF.z1) + "_" + str(KdeLF.z2) + ".h5"    
        if os.path.isfile(chain_name):
            print("The file '",chain_name,"' already exists, please clear or rename it before run a new MCMC!" )
            return
                
        if KdeLF.MLE_done is False and initial_point is None:
            print('Maximum likelihood estimation not done. It starts now...')             
            self.get_optimal_h()
        if initial_point is not None:
            if (self.adaptive and self.set_beta_fixed) or (self.adaptive is False):
                if len(initial_point) != 2:
                    print("initial_point should be a array or list contains two elements. e.g., [0.15,0.15]")
                    return
            else:
                if len(initial_point) != 3:
                    print("initial_point should be a array or list contains three elements. e.g., [0.15,0.15,0.1]")        
                    return            
            KdeLF.MLE_result=initial_point
            if self.data_loaded is False:
                self.load_sample()
        
        self.h1l, self.h1h =KdeLF.MLE_result[0]/5, 3*KdeLF.MLE_result[0]          #$$$$$$$$$$$$$$$$$$$$$$$$$$
        self.h2l, self.h2h =KdeLF.MLE_result[1]/5, 3*KdeLF.MLE_result[1]
        if len(KdeLF.MLE_result) == 3:
            if self.adaptive:
                self.bh=3*KdeLF.MLE_result[2]        
        
#        if self.adaptive is False:        
#            print('use uniform priors on h1,h2')        
        
        # Initialize the walkers
        npara=len(KdeLF.MLE_result)         
        pos = KdeLF.MLE_result + 1e-4 * np.random.randn(32, npara)  
        nwalkers, ndim = pos.shape
        backend = emcee.backends.HDFBackend(chain_name)
        backend.reset(nwalkers, ndim)
    
        # Initialize the sampler
        if parallel is True:                   
            from multiprocessing import Pool
            os.environ["OMP_NUM_THREADS"] = "1"
            kdeLF.kde_fortran_t.params.process=1
            with Pool() as pool:                
                sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob, backend=backend, pool=pool) 
        
                index = 0
                autocorr = np.empty(max_n)
    
                # This will be useful to testing convergence
                old_tau = np.inf     
                # Now we'll sample for up to max_n steps
                for sample in sampler.sample(pos, iterations=max_n, progress=True):
                    # Only check convergence every 50 steps
        
                    its=sampler.iteration
                    if its % 50 == 0:
                        f = open(logfilename, "a")
                        print('completed:', its,"/",max_n, ',', 'progress:','%.2f' % (its/max_n*100),"%",file=f)
                        f.close()   
                    if sampler.iteration % 50:
                        continue
    
                    # Compute the autocorrelation time so far
                    # Using tol=0 means that we'll always get an estimate even
                    # if it isn't trustworthy
                    tau = sampler.get_autocorr_time(tol=0)
                    autocorr[index] = np.mean(tau)
                    index += 1
    
                    # Check convergence
                    converged = np.all(tau * Ntau < sampler.iteration)
                    converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                    if converged:
                        break
                    old_tau = tau        
        else:        
            sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob, backend=backend)    
            # We'll track how the average autocorrelation time estimate changes
            index = 0
            autocorr = np.empty(max_n)
    
            # This will be useful to testing convergence
            old_tau = np.inf     
            # Now we'll sample for up to max_n steps
            for sample in sampler.sample(pos, iterations=max_n, progress=True):
                # Only check convergence every 50 steps
        
                its=sampler.iteration
                if its % 50 == 0:
                    f = open(logfilename, "a")
                    print('completed:', its,"/",max_n, ',', 'progress:','%.2f' % (its/max_n*100),"%",file=f)
                    f.close()   
                if sampler.iteration % 50:
                    continue
    
                # Compute the autocorrelation time so far
                # Using tol=0 means that we'll always get an estimate even
                # if it isn't trustworthy
                tau = sampler.get_autocorr_time(tol=0)
                autocorr[index] = np.mean(tau)
                index += 1
    
                # Check convergence
                converged = np.all(tau * Ntau < sampler.iteration)
                converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                if converged:
                    break
                old_tau = tau         
    
        n = 50 * np.arange(1, index + 1)
        y = autocorr[:index]
        plt.plot(n, n / float(Ntau), "--k", label=r"$\tau = N/$"+str(Ntau))
        plt.plot(n, y, "o-", label=r"$\tau = \hat{\tau}(N)$")
        plt.xlim(0, n.max())
        plt.ylim(0, y.max() + 0.1 * (y.max() - y.min()))
        plt.xlabel("number of samples, $N$")
        plt.ylabel(r"mean $\hat{\tau}$");        
        plt.legend(fontsize=12);        

        if self.adaptive is True:
            fig_name = "tau_" + str(KdeLF.z1) + "_" + str(KdeLF.z2) + "a.png"
        else:
            fig_name = "tau_" + str(KdeLF.z1) + "_" + str(KdeLF.z2) + ".png"
        plt.savefig(fig_name)
        
        self.chain_analysis(chain_name)
        return
        
    ##############################################################################
    def chain_analysis(self,chain_name):               
        import emcee        

        if os.path.isfile(chain_name) is False:
            print("The file '",chain_name,"' does not exist, please check the file name and path!" )
            return
            
        reader = emcee.backends.HDFBackend(chain_name)

        if self.adaptive:
            if self.set_beta_fixed:
                ndim=2
                labels = [r"$h_{10}$", r"$h_{20}$"]
                labs=[' h10',' h20']
            else:
                ndim=3
                labels = [r"$h_{10}$", r"$h_{20}$", r"$\beta$"]
                labs=[' h10',' h20','beta']
        else:
            ndim=2
            labels = [r"$h_{1}$", r"$h_{2}$"]
            labs=['h1','h2']


        fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
        samples = reader.get_chain()         
        for i in range(ndim):
            ax = axes[i]            
            ax.plot(samples[:, :, i], "k", alpha=0.3)   
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number");

        if self.adaptive is True:
            fig_name = "walkers_" + str(KdeLF.z1) + "_" + str(KdeLF.z2) + "a.png"
        else:
            fig_name = "walkers_" + str(KdeLF.z1) + "_" + str(KdeLF.z2) + ".png"
        fig.savefig(fig_name)

        try:
            tau = reader.get_autocorr_time()
            burnin = int(2 * np.max(tau))
            thin = int(0.5 * np.min(tau))
            flat_samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
            log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)        
            print("burn-in: {0}".format(burnin))
            print("thin: {0}".format(thin))
            print("flat chain shape: {0}".format(flat_samples.shape))
            print("flat log prob shape: {0}".format(log_prob_samples.shape))       
            dis=max(200,burnin*1.5)
            goodchain=True
        except:
            print("WARNING: The chain is shorter than 50 times the integrated autocorrelation time for", ndim, "parameter(s). Use this estimate with caution and run a longer chain!")
            dis=200
            goodchain=False            
        self.samples=samples[dis:,:,:].reshape((-1, ndim))            
        
        import corner
        figure = corner.corner(self.samples, labels=labels,
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_fmt='.3f', title_kwargs={"fontsize": 12})                       
         

        if self.adaptive is True:
            fig_name = "triangle_" + str(KdeLF.z1) + "_" + str(KdeLF.z2) + "a.png"
        else:
            fig_name = "triangle_" + str(KdeLF.z1) + "_" + str(KdeLF.z2) + ".png"
        figure.savefig(fig_name)
        
        #self.plot_posterior_LF()       
        
        #from IPython.display import display, Math                                
        self.mcmc_fit=np.zeros(ndim)
        print(' ')
        print('********** MCMC bestfit and 1 sigma errors **********')        
        for i in range(ndim):
            mcmc = np.percentile(self.samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
            txt = txt.format(mcmc[1], q[0], q[1], labels[i])            
            print(labs[i], '=','%.4f' %mcmc[1], '-', '%.4f' %q[0], '+', '%.4f' %q[1])
            #display(Math(txt))
            self.mcmc_fit[i]=mcmc[1]
        print(' ')    
        self.mcmc_fit_done=True        
        
        if goodchain is True:
            self.samples=flat_samples                                
        #print('self.samples',self.samples.shape)
        return     
        

    def phi_post(self,z, mag, theta):              
        phi=[]
        i=0
        for m in mag:
            lfi = self.log10phi(z,m,theta)
            i+=1
            phi.append(lfi) 
        return phi
        

    def plot_posterior_LF(self,z=None,sigma=1,bestfit_colour='red',error_band_colour='orange',pilot_bandwidths=None):
        from tqdm import tqdm
        if self.data_loaded is False:
            self.load_sample()        
        if z is None:
            z=np.mean(KdeLF.red)        
        if (z<=KdeLF.z1) or (z>=KdeLF.z2):
            print('INPUT ERROR: The input redshift should satisfy',str(KdeLF.z1),"< z <",str(KdeLF.z2))
            return        
        if self.adaptive is True:            
            if KdeLF.MLE_done is False and pilot_bandwidths is None:
                #print('Maximum likelihood estimation not done. It starts now...')             
                self.get_optimal_h()        
        
        nmags = 100
        if self.absolute_magnitude is True:
            mags=np.linspace(self.f_lim(z),KdeLF.Lmin,nmags)
        else:
            mags=np.linspace(self.f_lim(z),KdeLF.Lmax,nmags)        
        
        nsample = min(5000,len(self.samples))
        
        rsample = self.samples[np.random.randint(len(self.samples), size=nsample)]
        phi = np.zeros((nsample, nmags))
    
        for i, theta in enumerate(tqdm(rsample)):       
            phi[i] = self.phi_post(z,mags,theta)       
            
        if sigma==1:        
            up = np.percentile(phi, 15.87, axis=0)       
            down = np.percentile(phi, 84.13, axis=0)   
        if sigma==3:
            up = np.percentile(phi, 0.135, axis=0)         
            down = np.percentile(phi, 99.865, axis=0)     
        
        fig = plt.figure(figsize=(7, 7), dpi=100)
        ax = fig.add_subplot(1, 1, 1)
        ax.tick_params('both', which='major', length=7, width=1)
        ax.tick_params('both', which='minor', length=3, width=1)        
        ax.fill_between(mags, down, y2=up, color='orange', alpha=0.4) 
        
        bf = np.median(self.samples, axis=0)
        print('bestfit MCMC',bf)
        phi_fit = self.phi_post(z,mags,bf)  
        ax.plot(mags, phi_fit, lw=1, c='red', zorder=3)        

        lmax=max(KdeLF.lum)
        lmin=min(KdeLF.lum)
        LFmin=phi_fit[0]        
        print('lmax',lmax,lmin)        
        if self.absolute_magnitude is True:
            LFmax=np.log10(self.phi_kde(z,lmin,bf))           
        else:
            LFmax=np.log10(self.phi_kde(z,lmax,bf))
        y1=LFmax-0.25*(LFmin-LFmax)
        y2=LFmin+0.1*(LFmin-LFmax)
        y1=max(y1,phi_fit[len(phi_fit)-1])        
        
        if self.absolute_magnitude is True:
            tx=lmax-(lmax-KdeLF.Lmin)*0.618
            ty=y1+(y2-y1)*0.88
            ax.text(tx,ty,'z='+'%.3f' %z,fontsize=13,bbox=dict(boxstyle='square,pad=0.3', fc='yellow', ec='k',lw=0.5 ,alpha=0.4)) 
            #plt.plot([lmin,lmin],[-50,-0.1],linewidth=0.6, linestyle=(0, (5, 5, 1, 5)))
            plt.ylabel(r'$\log_{10}( \phi(z,M) ~/~ {\rm Mpc}^{-3} ~ \mathrm{mag}^{-1} )$',fontsize=18)
            plt.xlabel(r'$M$',fontsize=18)            
            plt.xlim(lmax,KdeLF.Lmin)
        else:
            tx=lmin+(KdeLF.Lmax-lmin)*0.618
            ty=y1+(y2-y1)*0.88
            ax.text(tx,ty,'z='+'%.3f' %z,fontsize=13,bbox=dict(boxstyle='square,pad=0.3', fc='yellow', ec='k',lw=0.5 ,alpha=0.4))
            plt.plot([lmax,lmax],[-20,-0.1],linewidth=0.6, linestyle=(0, (5, 5, 1, 5)))
            plt.ylabel(r'$\log_{10}( \phi(z,L) ~/~ {\rm Mpc}^{-3} ~ \Delta L^{-1} )$',fontsize=18)
            plt.xlabel(r'$\log_{10} ~ L$',fontsize=18)
            plt.xlim(lmin, KdeLF.Lmax)
        plt.ylim(y1,y2)        
        if self.adaptive is True:
            figname='lfa_z'+'%.3f' %z+'.png'
        else:
            figname='lf_z'+'%.3f' %z+'.png'
        plt.savefig(figname)
        #plt.show()
        result=np.array([mags,phi_fit,down,up])
        return result  



##############################################################################################################################################

#class LFsampler():
#    """sample a luminosity function"""    
#    def __init__(self, phi, z_bound, L_bound, f_lim, solid_angle,H0=70,Om0=0.3):        
#        self.phi=phi
#        self.z_bound=z_bound
#        self.omega=solid_angle
#        self.L_bound=L_bound
#        self.f_lim=f_lim 
#        self.H0=H0
#        self.Om0=Om0
#        self.cosmo = FlatLambdaCDM(H0=self.H0, Om0=self.Om0)

#    def __call__(self, x):   
#        return self.f_lim(x)
#        
#    def __call__(self, z,L):   
#        return self.phi(z,L)  
#   
#         
#    def rho(self,L,z):
#        dvdz=self.cosmo.differential_comoving_volume(z).value
#        result=self.phi(z,L)*omega*dvdz  
#        return result


#ans=dblquad(rho, z1, z2, lambda z: L1, lambda z: L2)
#print(ans)
#num=ans[0]
#ndata=math.ceil(num)
#print(num,ndata)

#def p(L,z):
#    dvdz=cosmo.differential_comoving_volume(z).value
#    return phi_real(z,L)*dvdz*omega/num

#ans=dblquad(p, z1, z2, lambda z: L1, lambda z: L2)
#print(ans)        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
              



