import numpy as np
import sys

class globalSLRModel:
    '''This is a model to calculate global SLR based on the contributions from five components.
    Input to this model are only the global mean surface air temperature anomaly and ocean heat content changes.

    This is based on a variety of different models (see description of individual components).

    References:
    Perette et al. (2013): "A scaling approach to project regional sea level rise and its uncertainties"
    Marti et al. (2022)
    Levitus et al. (2012)
    Church et al. (2011)
    Wong et al. (2017) (BRICKv0.3 / mimiBRICK)
    Fox-Kemper et al. (2021) (IPCC report, Chapter 9 (?))
    Nauels et al. (2017) (MAGICC6)
    Li et al. (2020)
    Wigley and Raper (2005)
    Shaffer et al. (2014)

    ##########################################################
    #### Thermosteric SLR as a function of ocean energy input
    ##########################################################
    Thermosteric SLR follows ocean heat content changes
    approximately linearly with rate factor r (m/YJ)
    Use here r=0.11 because of large possible range of values
    and because it fits best with IPCC estimates.
    Reference values are:
    0.145 (Marti et al.2022)
    0.12 (Levitus et al., 2012)
    0.15 (Church et al., 2011)
    0.105 (BRICK v0.3)
    0.1199 from IPCC estimates
    0.115 from MPIESM fit
    0.0975 from fitting OBS data of OHC to thermoSLR

    ##########################################################
    ##### SLR from Mountain glaciers #########################
    ##########################################################
    Calculating the mountain glacier component as in Perrette et al. (2013)

    ##########################################################
    ##### SLR from land water storage change #################
    ##########################################################
    Very little data exist. This here is just a very simple
    first implementation taken over from BRICKv0.3 (Wong, 2017).
    -> BRICK just describes a linear increase with some 
       random fluctuations (which only starts in year 2000)

    ##########################################################
    ##### SLR from Greenland Ice Sheet #######################
    ##########################################################
    Calculation as in MAGICC6 (Nauels, 2017).
    Implementation separates the two components and adds them together.
    In the discharge component we have a choice between a low and
    a high discharge scenario, default is high, because it
    seems to fit better with Observational data.

    ##########################################################
    ##### SLR from Antarctic Ice Sheet #######################
    ##########################################################
    Calculation of Antarctic ice sheet volume based on the implementation in BRICK (Wong et al., 2017).
    This employs the DAIS model (Shaffer, 2014)
    '''

    def __init__(self, T, OHC_change, sy=1750, ey=2300, dt=1.0, dbg=0):

        ########## Set some global variables ################
        self.T_anomaly = T
        self.OHC_change = OHC_change

        self.dbg = dbg
        
        # Timestep:
        self.sy = sy
        self.ey = ey
        self.dt = dt
        self.nyears = int((self.ey - self.sy) / self.dt) + 1
        self.time = np.linspace(self.sy+self.dt*0.5, self.ey-self.dt*0.5, self.nyears)

        ########## Parameter settings ###################

        ####################################################################################
        ###########          Tune these parameters          ################################
        ####################################################################################
        ### Tuning for thermosteric component
        self.thermo_m_per_YJ = 0.11

        ### Tuning for Land water storage component
        self.LWS_rate = 0.0003

        ### Tuning for Mountain glacier component


        ### Tuning for Greenland Ice Sheet component
        self.GIS_rho = 7.933e-4
        self.GIS_eps = 0.4722
        self.GIS_outlet_max = 0.05363 


        ### Tuning for Antarctic Ice Sheet component
        self.AIS_gamma = 2.
        self.AIS_alpha = 0.35
        self.AIS_mu = 8.7
        self.AIS_nu = 0.012
        self.AIS_P0 = 0.35
        self.AIS_kappa = 0.04
        self.AIS_f0 = 1.2
        self.AIS_h0 = 1700.
        self.AIS_c = 95.
        self.AIS_b0 = 775.
        self.AIS_slope = 0.0006
        


        ####################################################################################
        ####################################################################################

        ## Other parameters
        # thermosteric component
        self.thermo_startyear = self.sy

        # Land water storage component
        self.LWS_shape = 0.00018
        self.LWS_startyear = 2000

        # Mountain glacier component
        self.MG_startyear = 1850
        self.MG_T_eq = 0.0
        self.MG_beta_0 = 0.0008
        self.MG_n = 1.646
        self.MG_V_0 = 0.41

        # GrIS component
        self.GIS_startyear = self.sy
        self.GIS_s = 5
        self.GIS_v = 0.0001148
        self.GIS_X = 0.0
        self.GIS_p = 2.0169
        self.GIS_SMB_max = 7.36    # m

        # AntIS component
        self.AIS_startyear = self.sy+1
        self.AIS_rho_w = 1030.
        self.AIS_rho_i = 917.
        self.AIS_rho_m = 4000.
        self.AIS_Rad0 = 1.864e6
        self.AIS_dell = self.AIS_rho_w / self.AIS_rho_i
        self.AIS_eps1 = self.AIS_rho_i /(self.AIS_rho_m - self.AIS_rho_i)
        self.AIS_eps2 = self.AIS_rho_w /(self.AIS_rho_m - self.AIS_rho_i)

        ########## Initialize arrays ################
        self.SLR_thermo = np.zeros((self.nyears))
        self.SLR_LWS    = np.zeros((self.nyears))
        self.SLR_MG     = np.zeros((self.nyears))
        self.SLR_GIS    = np.zeros((self.nyears))
        self.SLR_AIS    = np.zeros((self.nyears))
        self.SLR_total  = np.zeros((self.nyears))

        self.SLR_GIS_SMB = np.zeros((self.nyears))
        self.SLR_GIS_DIS = np.zeros((self.nyears))
        self.GIS_outlet_vdis = np.zeros((self.nyears))


        self.AIS_Volume = np.zeros((self.nyears))
        self.AIS_Radius = np.zeros((self.nyears))


    def reset_SLR(self):
        ########## Initialize arrays again ################
        self.SLR_thermo = np.zeros((self.nyears))
        self.SLR_LWS    = np.zeros((self.nyears))
        self.SLR_MG     = np.zeros((self.nyears))
        self.SLR_GIS    = np.zeros((self.nyears))
        self.SLR_AIS    = np.zeros((self.nyears))
        self.SLR_total  = np.zeros((self.nyears))

        self.SLR_GIS_SMB = np.zeros((self.nyears))
        self.SLR_GIS_DIS = np.zeros((self.nyears))
        self.GIS_outlet_vdis = np.zeros((self.nyears))

        self.AIS_Volume = np.zeros((self.nyears))
        self.AIS_Radius = np.zeros((self.nyears))


    def align(self, year):
        index = int(year - self.sy)
        print('Aligning SLR to be 0 in year '+str(year)+' (index: '+str(index)+')')
        self.SLR_thermo = self.SLR_thermo - self.SLR_thermo[index]
        self.SLR_LWS = self.SLR_LWS - self.SLR_LWS[index]
        self.SLR_MG = self.SLR_MG - self.SLR_MG[index]
        self.SLR_GIS = self.SLR_GIS - self.SLR_GIS[index]
        self.SLR_AIS = self.SLR_AIS - self.SLR_AIS[index]
        self.SLR_total = self.SLR_total - self.SLR_total[index]


    def integrate(self, silent=True):
        if not silent: print('Start integrating...')

        self.__AIS_init()
        self.__GIS_init()
        for i in range(0, self.nyears-1):
            if self.dbg==1: print('Year: ', i)
            if self.time[i] >= self.thermo_startyear: self.__update_SLR_thermo(i)
            if self.time[i] >= self.LWS_startyear:    self.__update_SLR_LWS(i)
            if self.time[i] >= self.MG_startyear:     self.__update_SLR_MG(i)
            if self.time[i] >= self.GIS_startyear:    self.__update_SLR_GIS(i)
            if self.time[i] >= self.AIS_startyear:    self.__update_SLR_AIS(i)

            self.SLR_total[i+1] = self.SLR_thermo[i+1] + self.SLR_LWS[i+1] + self.SLR_MG[i+1] + self.SLR_GIS[i+1] + self.SLR_AIS[i+1]

        if not silent: print('   ...finished')



    def __update_SLR_thermo(self,i):
        if self.dbg==1: print('   SLR thermo: ', i)
        self.SLR_thermo[i+1] = self.SLR_thermo[i] + self.thermo_m_per_YJ * 1.e-24 * self.OHC_change[i]
        return

    def __update_SLR_LWS(self, i):
        if self.dbg==1: print('   SLR LWS: ', i)        
        self.SLR_LWS[i+1] = self.SLR_LWS[i] + np.random.normal(loc=self.LWS_rate, scale=self.LWS_shape)
        return

    def __update_SLR_MG(self, i):
        if self.dbg==1: print('   SLR MG: ', i)
        change = self.MG_beta_0 * (self.T_anomaly[i] - self.MG_T_eq) * (1.0 - self.SLR_MG[i]/self.MG_V_0)**self.MG_n
        self.SLR_MG[i+1] = self.SLR_MG[i] + change
        return


    def __GIS_init(self):
        self.GIS_outlet_vdis[:] = self.GIS_outlet_max
        return

    def __update_SLR_GIS(self, i):
        if self.dbg==1: print('   SLR GIS: ', i)
        T = self.T_anomaly[i]

        if T >= 0.0: term1 = (self.GIS_X*T + (1.-self.GIS_X)*T**self.GIS_p)
        else: term1 = 0.0
        term2 = np.sqrt(1.0 - self.SLR_GIS_SMB[i]/self.GIS_SMB_max)
    

        if self.dbg == 'GIS': print(i, 'Tanomaly:', T)
        if self.dbg == 'GIS': print(i, 'term1:', term1)
        if self.dbg == 'GIS': print(i, 'term2:', term2)
        

        self.SLR_GIS_SMB[i+1] = self.SLR_GIS_SMB[i] + self.GIS_v * term1 * term2
        if self.dbg == 'GIS': print(i, 'SLR_GIS_SMB:',  self.SLR_GIS_SMB[i+1])

        tmp = self.GIS_rho * self.GIS_outlet_vdis[i] * np.exp(self.GIS_eps*T)
        if tmp < 0: tmp = 0.0

        self.GIS_outlet_vdis[i+1] = self.GIS_outlet_vdis[i] - tmp
        self.SLR_GIS_DIS[i+1] = self.GIS_s*(self.GIS_outlet_max - self.GIS_outlet_vdis[i+1])

        self.SLR_GIS[i+1] = self.SLR_GIS_DIS[i+1] + self.SLR_GIS_SMB[i+1]

        if self.dbg == 'GIS':
            print(i, 'tmp:', tmp)
            print(i, 'GIS_outlet_vdis:', self.GIS_outlet_vdis[i+1] )
            print(i, 'SLR_GIS_DIS:', self.SLR_GIS_DIS[i+1])
            print(i, 'SLR_GIS:', self.SLR_GIS[i+1])


        return

    def __AIS_init(self):
        R = self.AIS_Rad0
        rc = self.AIS_b0/self.AIS_slope
        V = np.pi * (1+self.AIS_eps1) * ( (8./15.) * self.AIS_mu**0.5 * R**2.5 - (1./3.)*self.AIS_slope*R**3)
        if R>rc: V = V - np.pi*self.AIS_eps2 * ( (2./3.)  * self.AIS_slope*(R**3-rc**3)-self.AIS_b0*(R**2-rc**2) )
        self.AIS_Volume[:] = V
        self.AIS_Radius[:] = R

        if self.dbg == 'AIS':
            print('rc init: ', rc)
            print('V init: ', V)
            print('V init: ', V)
            print('R init: ', R)

        return


    def __update_SLR_AIS(self, i):
        if self.dbg==1: print('   SLR AIS: ', i)
        Tg = self.T_anomaly[i]
        R = self.AIS_Radius[i]
        V = self.AIS_Volume[i]
        SL = self.SLR_total[i]
        dSL = SL - self.SLR_total[i-1]
        AIS_includes_dSLais = 1.0

        # define some constants
        Pi = np.pi
        Tf = -1.8
        Toc_0 = 0.72
        Aoc = 3.619e14
        lf    = -1.18

        # From fitting Antarctic air temperature to global mean temperature anomaly
        c1 = 14.27863951
        c2 = 0.77385984
        Ta = (Tg - c1) / c2
    
        # Connecting Antarctic ocean temperature to Antarctic air temperature
        a_anto = 0.3
        b_anto = 0.5
        Toc = Tf + (a_anto*Tg + b_anto-Tf) / (1. + np.exp(-Tg+(Tf-b_anto)/a_anto))

        if self.dbg == 'AIS': print(i, 'SL:', SL)
        if self.dbg == 'AIS': print(i, 'dSL:', dSL)
        if self.dbg == 'AIS': print(i, 'Ta:', Ta)
        if self.dbg == 'AIS': print(i, 'Toc:', Toc)

        # Start model
        hr   = self.AIS_h0 + self.AIS_c * Ta        # equation 5
        rc   = (self.AIS_b0 - SL)/self.AIS_slope    # application of equation 1 (paragraph after eq3)
        P    = self.AIS_P0 * np.exp(self.AIS_kappa*Ta) # equation 6
        beta = self.AIS_nu * P**(0.5)      # equation 7 (corrected with respect to text)

             
        rR = R - ((hr - self.AIS_b0 + self.AIS_slope*R)**2) / self.AIS_mu
        Btot_1 = P * Pi * R**2 - \
        Pi * beta * (hr - self.AIS_b0 + self.AIS_slope*R) * (R*R - rR*rR) - \
        (4. * Pi * beta * self.AIS_mu**0.5 *   (R-rR)**2.5) / 5.  + \
        (4. * Pi * beta * self.AIS_mu**0.5 * R*(R-rR)**1.5) / 3.
        Btot_2 = P * Pi*R**2
             
        if hr > 0: Btot = Btot_1
        else: Btot = Btot_2

        if self.dbg == 'AIS': print(i, 'Btot:', Btot)
                                  
        F_1   = 0.   # no ice flux
        ISO_1 = 0.   # (third term equation 14) NAME?
        fac_1 = Pi * (1.+self.AIS_eps1) * (4./3. * self.AIS_mu**0.5 * R**1.5 - self.AIS_slope*R**2) # ratio dV/dR (eq 14) 

        ## if R > rc
        # In case there is a marine ice sheet / grounding line
        fac_2   = fac_1 - 2.*Pi*self.AIS_eps2 * (self.AIS_slope*R**2 - self.AIS_b0*R) # correction fac (eq 14)
        Hw = self.AIS_slope*R - self.AIS_b0 + SL  # equation 10
        
        # Ice speed at grounding line (equation 11)
        Speed = self.AIS_f0 * ((1.-self.AIS_alpha) + self.AIS_alpha * ((Toc - Tf)/(Toc_0 - Tf))**2) * \
                (Hw**self.AIS_gamma) / ( (self.AIS_slope*self.AIS_Rad0 - self.AIS_b0)**(self.AIS_gamma-1.) )
        
        F_2   = 2.*Pi*R * self.AIS_dell * Hw * Speed   # equation 9

        # ISO term depends on dSL_tot (third term equation 14 !! NAME)
        c_iso = 2.*Pi*self.AIS_eps2* (self.AIS_slope*rc**2 - (self.AIS_b0/self.AIS_slope)*rc)

        # first term is zero if dSL represents only non-AIS components (dSLais=0)
        # second term is zero if dSL represents all components (dSLais=1)
        ISO_2 = AIS_includes_dSLais * c_iso *  dSL + (1.-AIS_includes_dSLais) * \
                    ((1.-c_iso)/c_iso) * (dSL - lf * (Btot - F_2) / Aoc)
                 
        if R>rc:
            F = F_2
            ISO = ISO_2
            fac = fac_2
        else:
            F = F_1
            ISO = ISO_1
            fac = fac_1

        self.AIS_Radius[i+1] = R + (Btot-F+ISO)/fac
        self.AIS_Volume[i+1] = V + (Btot-F+ISO)

        self.SLR_AIS[i+1] = 57. * (1. - self.AIS_Volume[i+1]/self.AIS_Volume[0])

        if self.dbg == 'AIS':
            print(i, 'F:', F)
            print(i, 'ISO:', ISO)
            print(i, 'fac:', fac)

            print(i, 'R: ', self.AIS_Radius[i+1])
            print(i, 'V: ', self.AIS_Volume[i+1])

            print(i, 'SLR_AIS: ', self.SLR_AIS[i+1])

            print()



        return


    def getSLRTotal(self):    return self.SLR_total
    def getSLRThermo(self):    return self.SLR_thermo
    def getSLRLWS(self):    return self.SLR_LWS
    def getSLRMG(self):    return self.SLR_MG
    def getSLRGIS(self):    return self.SLR_GIS
    def getSLRAIS(self):    return self.SLR_AIS
    def getTime(self):         return self.time
    def getTanomaly(self):     return self.T_anomaly
    def getOHCchange(self):     return self.OHC_change
