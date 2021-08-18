import numpy as np
import nucrat
import os
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
import Opacity

def properties(stellar_m, f_d, R_g):

        # constants in SI units
        kB  = 1.3806*10**(-23.0) #m^2*kg*s^(-2)*K^(-1)
        GN  = 6.67*10**(-11.0) #m^3*kg^(-1)*s^(-2)
        pi4 = 4.0 * np.pi
        eps = 10.0**(-10)
        rhoconvert = 5.62*10**20 #1kg/m^3 = 5.62*10^20 GeV/cm^3

        #Initialize values for radius finder
        #factor = 1.0*10**(11)
        prec   = 0.05
        firstpass = True
        printing  = False

        solarmass = 1.988*10**30.0 # kg solar mass
        GeV       = 1.602*10**(-10.0) #J
        parsec    = 3.086*10**(16) # m
        solarluminosity = 3.828*10**26.0 #W solar luminosity
        m_u       = 1.6605*10**(-27.0) #kg atomic unit
        m_H       = 1.007825*m_u
        c         = 3.0*10**8.0 #m/s
        sigma     = 5.6704*10**(-8.0) #W*m^(-2)*K^(-4) stefan-boltzmann constant
        yearinsec = 31536000.0 #s one year in seconds
        AU = 1.496*10**11 #m one au in meters

        # component constant
        X = 0.76
        Y = 0.24
        mu = (2*X+0.75*Y)**(-1.0)

        starmin = 8.5*10**32.0 #kg below starmin, n=1.5
        starmax = 2.0*10**33.0 #kg above starmax, n=3.0

        ######
        ###
        ###  Quantities to set
        ###
        ######

        # Halo model
        rho_halo = solarmass/parsec**3
        vbar     = 10.0**4 # m/s
        z        = 20.0
        fz       = 1.14*10**7
        rho_halo = fz*rho_halo


        # Star
        # M        = 500.0*solarmass
        # Rguess   = 10.0**7*M/solarmass*200.0
        # fdm      = 0.01
        dotM     = 10.0**(-3)*solarmass
        DeltaM   = 1.0*solarmass #mass increment for mass building up
        DltTm    = DeltaM*yearinsec/dotM

        # Dark matter model
        m_DM     = 100*GeV/c**2 #kg mass of DM particles converted to kg
        m_phi    = 0.01*GeV/c**2 #kg mass of scalar mediator converted to kg
        f_EM     = 0.66666
        sigma_v  = 3.0*10**(-32.0) #m^3/s average annihilation cross section of dark matter particles
        sigma_SD = 10.0**(-44) #m^2 
        sigma_ss = 1.0**(0)   #cm^2/g
        sigma_ss = sigma_ss*0.1*m_DM #conversion to m^2

        ######
        ###
        ### End of the quantities to set
        ###
        ######
        x = np.array([-2.398, -1.398, -0.398, 0.699, 1.699, 2.74])
        y = np.array([-2.523, -2.0, -1.523, -1.0, -0.523, 0.0])
        f = interp1d(x,y, fill_value='extrapolate') # f is g_psy which produces correct relic density, dependent on DM mass
        alpha = (10**f(2))**2/(4*np.pi) #
        M = stellar_m*solarmass
        fdm = f_d
        Rguess   = R_g
        # temperature finder
        def tsub(rho,P,guess):
                t = guess
                Not_Good = True
                while Not_Good:
                        fx = P-rho/(m_u*mu)*kB*t-1.0/3.0*sigma/c*t**4
                        dfx = -rho/(m_u*mu)*kB-4.0/3.0*sigma/c*t**3
                        dx = -fx/dfx
                        t = t+0.5*dx
                        if abs(dx/t) <= 1.0*10**(-5.0):
                                Not_Good = False
                return t

        poly = np.zeros(4)
        ksi = np.zeros(4)
        deri = np.zeros(4)
        ratio = np.zeros(4)
        # polytrope star model parameters
        poly[0] = 1.5 #n
        ksi[0] = 3.65375 #zn
        deri[0] = 2.71406 #theta
        ratio[0] = 5.99071 #rho_ratio

        poly[1] = 2.0
        ksi[1] = 4.35287
        deri[1] = 2.41105
        ratio[1] = 11.40254

        poly[2] = 2.5
        ksi[2] = 5.35528
        deri[2] = 2.18720
        ratio[2] = 23.40646

        poly[3] = 3.0
        ksi[3] = 6.89685
        deri[3] = 2.01824
        ratio[3] = 54.1825
        # parameter determines whether the star reaches maximum mass or Eddington luminosity
        big_enough = False
        # total lifetime of a dark star
        total_time = 0.0

        Rold    = 0.0
        L_GCold = 0.0
        rho0    = 0.0

        R = Rguess
        # timecount = timecount + 1
        if M <= starmin:
                n = poly[0]
                zn = ksi[0]
                theta = deri[0]
                rho_ratio = ratio[0]
        elif M > starmin and M <=starmax:
                n = 1.5 + (3.0-1.5)*(M-starmin)/(starmax-starmin)
                if n <= 2.0:
                        zn = ksi[0]+(ksi[1]-ksi[0])*(n-poly[0])/(poly[1]-poly[0])
                        theta = deri[0]+(deri[1]-deri[0])*(n-poly[0])/(poly[1]-poly[0])
                        rho_ratio = ratio[0]+(ratio[1]-ratio[0])*(n-poly[0])/(poly[1]-poly[0])
                elif n > 2.0 and n <=2.5:
                        zn = ksi[1]+(ksi[2]-ksi[1])*(n-poly[1])/(poly[2]-poly[1])
                        theta = deri[1]+(deri[2]-deri[1])*(n-poly[1])/(poly[2]-poly[1])
                        rho_ratio = ratio[1]+(ratio[2]-ratio[1])*(n-poly[1])/(poly[2]-poly[1])
                else:
                        zn = ksi[2]+(ksi[3]-ksi[2])*(n-poly[2])/(poly[3]-poly[2])
                        theta = deri[2]+(deri[3]-deri[2])*(n-poly[2])/(poly[3]-poly[2])
                        rho_ratio = ratio[2]+(ratio[3]-ratio[2])*(n-poly[2])/(poly[3]-poly[2])
        else:
                n = poly[3]
                zn = ksi[3]
                theta = deri[3]
                rho_ratio = ratio[3]

        power = np.arange(10,1,step=-0.5)
        masses = 10**power
        radius = np.zeros(len(power))
        lum = np.zeros(len(power))
        lum_perc = np.zeros(len(power))
        temp = np.zeros(len(power))
        heating_para = 1.0/masses
        # parameter indicating whether two luminosities are close enough
        Not_Close = True
        # factor    = 1.0*10**(10)
        firstpass = True
        incR      = True
        decR      = True
        # m_DM = masses[k]*GeV/c**2
        resetR = True

        count = 0
        while Not_Close:
                #R = Rguess
                rho_bar = 3.0*M/(pi4*R**3.0)
                rho_cg  = rho_ratio*rho_bar
                rho_cd  = fdm*rho_cg
                # rho_cd  = rho_halo
                # fdm     = rho_cd / rho_cg
                K_g = pi4*GN*(R/zn)**2*rho_cg**(1.0-1.0/n)/(n+1.0)
                K_d = GN*M/(3.0*R) # (m/s)**2
                v_disp = K_d**0.5 / (3.0*10**8) # velocity dispersion of DM in the DS, useful when considering sommerfeld enhancement
                a = v_disp/(2*alpha)
                c_param = 6*alpha*m_DM/(np.pi**2*m_phi)
                S = np.pi/a*np.sinh(2*np.pi*a*c_param)/(np.cosh(2*np.pi*a*c_param)-np.cos(2*np.pi*np.sqrt(c_param-a**2*c_param**2)))*((c_param-1)**2+4*a**2*c_param**2)/(1+4*a**2*c_param**2) # Sommerfeld enhancement term 
                sigma_v = S*3.0*10**(-32.0)*(v_disp/0.1)**2/1.549 # Sommerfeld enhancement goes like 1/v for p-wave annihilations, the velocity at relic density time is 0.1, normalization factor 1.549 is the enhancement at v = 0.1
                alpha_g = (zn/R)**2           # 1/length**2 unit
                alpha_d = pi4*GN*rho_cg/K_d   # 1/length**2 unit
                #if printing:
                #    print('rho_bar=',rho_bar)
                #    print('R=',R)
                #    print('M=',M)
                #    print('alpha_d=',alpha_g)

                # model to solve
                def model(z,r):
                        rho_g = z[0]
                        rho_d = z[1]
                        phi_g = z[2]
                        phi_d = z[3]
                        drho_gdr = phi_g
                        drho_ddr = phi_d
                        dphi_gdr = -alpha_g*(rho_g**n + fdm*rho_d) - 2.0*phi_g/r
                        dphi_ddr = -alpha_d*(rho_g**n + fdm*rho_d)*rho_d - 2.0*phi_d/r + phi_d**2/rho_d
                        dzdr = [drho_gdr,drho_ddr,dphi_gdr,dphi_ddr]
                        return dzdr

                # initial condition
                z0 = [1,1,0,0]

                # number of radius points
                n_p = 15001

                # radius points
                #r = np.linspace(10,R,n_p)
                #delta_r = r[1]-r[0]

                n_p1=4000
                n_p2=4001
                n_p = n_p1+n_p2
                #r = np.logspace(1,np.log10(R),n_p)
                r1 = np.logspace(1,np.log10(0.8*R),n_p1)
                r2 = np.logspace(np.log10(0.81*R),np.log10(R),n_p2)
                r = [*r1, *r2]
                #with open('file_r.txt', 'w') as f:
                #  for item in r:
                #       f.write("%s\n" % item)
                r = np.array(r)
                delta_r = np.zeros(r.size)
                delta_r[0] = r[0]
                for i in range(0,n_p-1):
                        delta_r[i+1] = r[i+1]-r[i]
                #print("r=",r)
                #print("delta_r=",delta_r)

                # store solutions
                rho_g = np.zeros(r.size)
                rho_d = np.zeros(r.size)
                phi_g = np.zeros(r.size)
                phi_d = np.zeros(r.size)

                # initial conditions
                rho_g[0] = z0[0]
                rho_d[0] = z0[1]
                phi_g[0] = z0[2]
                phi_d[0] = z0[3]

                # solve our group of ODE
                for i in range(1,n_p):
                        r_span = [r[i-1],r[i]]
                        # integrate to next step
                        z = odeint(model,z0,r_span)
                        # store solution
                        rho_g[i] = z[1][0]
                        rho_d[i] = z[1][1]
                        phi_g[i] = z[1][2]
                        phi_d[i] = z[1][3]
                        # next initial condition
                        z0 = z[1]

                # eliminate all the nan's in both profiles
                rho_g = rho_cg*rho_g[~np.isnan(rho_g)]**n
                rho_d = rho_cd*rho_d[~np.isnan(rho_d)]

                # find the mass profile (useful in finding photosphere)
                stellar_mass = np.zeros(rho_g.size)
                stellar_mass[0] = pi4/3.0*r[0]**3*rho_g[0]
                for i in range(rho_g.size-1):
                        stellar_mass[i+1] = (stellar_mass[i] + \
                                        pi4/3.0*(r[i+1]**3-r[i]**3)*rho_g[i+1])

                # integrated number of DM particles
                NDM = np.zeros(rho_d.size)
                NDM[0] = pi4/3.0*r[0]**3*rho_d[0]
                for i in range(rho_d.size-1):
                        NDM[i+1] = NDM[i] + pi4/3.0*(r[i+1]**3-r[i]**3)*rho_d[i+1]
                NDM = NDM / m_DM

                # define the integrated luminosity
                totlumin = np.zeros(rho_g.size)

                # find the temperature profile and pressure profile
                T = np.zeros(rho_g.size)
                P = K_g*rho_g**(1.0+1.0/n)
                for i in range(rho_g.size):
                        T[i] = tsub(rho_g[i],P[i],P[i]*mu*m_u/rho_g[i]/kB)

                with open('file_r.txt', 'w') as f:
                        for item in r:
                                f.write("%s\n" % item)
                with open('file_t.txt', 'w') as f:
                        for item in T:
                                f.write("%s\n" % item)
                with open('file_g.txt', 'w') as f:
                        for item in rho_g:
                                f.write("%s\n" % item)
                with open('file_d.txt', 'w') as f:
                        for item in rho_d:
                                f.write("%s\n" % item)

                # find the position of photosphere
                photo = np.zeros(rho_g.size)
                for i in range(rho_g.size):
                        photo[i] = Opacity.opac(np.log10(rho_g[i]),T[i])*P[i]-2.0/3.0*GN*M/R**2

                R_photo = 0.0
                T_eff = 0.0
                for i in range(rho_g.size-1):
                        if photo[i]>=0 and photo[i+1]<=0:
                                R_photo = r[i]
                                T_eff   = T[i]
                                M_tot   = stellar_mass[i]
                                im      = i
                                break
                for i in range(rho_d.size-1):
                        if R_photo <= r[i]:
                                N_DM    = NDM[i]
                                break

                # This part is about capture and self-capture, which is defuted by other researchers like Haibo et al.
                # C_ann = rho_d[0]**2*pi4/3.0*r[0]**3
                # for i in range(rho_d.size-1):
                #      C_ann = C_ann + pi4/3.0*(r[i+1]**3-r[i]**3)*rho_d[i+1]**2
                # C_ann = 2.0*sigma_v*C_ann/(m_DM*N_DM)**2

                # # Compute the capture rate
                # B = 1.5*GN*stellar_mass[0]/(r[0]*vbar**2)
                # coeff = B - 1 + np.exp(-B)
                # C_cap = rho_g[0]*coeff*r[0]**3
                # for i in range(rho_g.size-1):
                #      B     = 1.5*GN*stellar_mass[i+1]/(r[i+1]*vbar**2)
                #      coeff = B - 1 + np.exp(-B)
                #      C_cap = C_cap + rho_g[i+1]*coeff*(r[i+1]**3-r[i]**3)
                # C_cap = C_cap*2.0*rho_halo*X/m_H/m_DM*sigma_SD*vbar/np.sqrt(np.pi)*pi4/3.0

                # # Compute the rate of self-capture
                # B = 1.5*GN*stellar_mass[0]/(r[0]*vbar**2)
                # coeff = B - 1 + np.exp(-B)
                # C_self = rho_d[0]*coeff*r[0]**3
                # for i in range(rho_d.size-1):
                #      B      = 1.5*GN*stellar_mass[i+1]/(r[i+1]*vbar**2)
                #      coeff  = B - 1 + np.exp(-B)
                #      C_self = C_self + rho_d[i+1]*coeff*(r[i+1]**3-r[i]**3)
                # C_self = 2.0*rho_halo*vbar/np.sqrt(np.pi)*pi4/3.0*C_self*sigma_ss/m_DM**2/N_DM

                # tau_EQ    = 1.0/np.sqrt((C_self/2.0)**2 + C_cap*C_ann)
                # tnh       = np.tanh(DltTm/tau_EQ)
                # if DltTm < 20.0*tau_EQ :
                #      N_cap = C_cap * tnh/ (1/tau_EQ - (C_self/2.0)*tnh)
                # else :
                #      N_cap = C_cap * tnh/ (C_cap*C_ann/C_self + C_self/2.0*(1.0-tnh))
                # Gamma_ann = 0.5*C_ann*N_cap**2
                # tau_TH    = m_DM/(rho_d[0]*sigma_ss*vbar)

                # if True :
                #   print('C_cap  = ',C_cap  )
                #   print('C_self = ',C_self )
                #   print('C_ann  = ',C_ann  )
                #   print('Regime = ', C_self/np.sqrt(C_cap*C_ann))
                #   print('tau_EQ = ',tau_EQ )
                #   print('time   = ',DltTm )
                #   print('N_cap  = ',N_cap)
                #   print('N_DM   = ',N_DM)
                #   print('G_ann  = ',Gamma_ann)
                

                # calculate effective luminosity
                L_eff = pi4*R_photo**2.0*sigma*T_eff**4.0

                # include luminosity contributed from different sources:
                # DM annihilation
                L_DM = rho_d[0]**2*pi4/3.0*r[0]**3
                m_dark = rho_d[0]*pi4/3.0*r[0]**3
                m_bar = rho_g[0]*pi4/3.0*r[0]**3
                for i in range(rho_d.size):
                        if r[i] <= R_photo :
                                L_DM = L_DM + rho_d[i]**2*pi4/3.0*(r[i+1]**3-r[i]**3)
                                m_dark += rho_d[i]*pi4/3.0*(r[i+1]**3-r[i]**3)
                                m_bar += rho_g[i]*pi4/3.0*(r[i+1]**3-r[i]**3)
                L_DM = sigma_v*L_DM*c**2/m_DM

                mass_ratio = m_dark/m_bar

                #print('sigma_v = ', sigma_v)
                #print('c = ', c)
                #print('m_DM = ', m_DM)
                #print('rho_DM = ', rho_d[0])
                #print('L_DM = ', L_DM)

                # nuclear reaction
                L_NR = 0.0
                for i in range(rho_d.size):
                        rhog  = 10**(-3)*rho_g[i] #convert the density from kg/m^3 to g/cm^3 
                        num_p = rhog*X/m_u
                        num_d = num_p/7000.0
                        num_He3 = 0.0
                        nucr = nucrat.nucrat(rhog,T[i],num_p,num_d,num_He3,DltTm,0.24,0.0,0.0,0.0)
                        L_NR = L_NR + nucr*10.0**(-4.0)*pi4*r[i]**2.0*delta_r[i]*rho_g[i]

                # gravitational contraction
                L_GC  = 0.0
                #for i in range(rho_g.size):
                #      if r[i] <= R_photo :
                #            L_GC = L_GC + (pi4*rho_g[i]*r[i]**2)**2*delta_r[i]*GN/3.0
                #L_GC = (L_GC-L_GCold)/DltTm
                #L_GC = 0.6*GN * M * dotM/R_photo/yearinsec

                # Dark matter capture heating
                # L_cap = f_EM*Gamma_ann*2.0*m_DM*c**2
                L_cap = 0.0

                L_total = L_DM + L_NR + L_GC + L_cap
                diff    = (L_total-L_eff)/(L_total+L_eff)

                count = count + 1
                if True:
                        print("M=",stellar_m)
                        print("f_d=",f_d)
                        print("mass_ratio=",mass_ratio)
                        print("R=",R)
                        print("R_photo=",R_photo)
                        print("T_eff=",T_eff)
                        print("L_eff=",L_eff)
                        print("L_DM=",L_DM)
                        print("L_cap=",L_cap)
                        print('L_GC = ',L_GC)
                        print('L_NR = ',L_NR)
                        print("L_total=",L_total)
                        print('ratio=',diff)
                        print('count=',count)
                #exit()
                # if count > 20:
                #         prec = 0.1
                aa = 0.05
                diff = -diff
                low_prec = 0.1
                factor = 2*aa*R
                if firstpass:
                        if diff > prec :
                                if decR : 
                                        incR = False
                                        R = R - factor
                                else:
                                        firstpass = False
                                        Rhigh = R
                                        Rlow = R - factor
                                        R = 0.5*(Rhigh+Rlow)
                                # if factor > aa*R :
                                #         factor = aa*R
                                

                                f = open('fileR.txt', 'a')
                                f.write("%s\n" % R)
                        elif diff < -prec :
                                print("R should increase")
                                if incR :
                                        decR= False
                                        R = R + factor
                                else:
                                        firstpass = False
                                        Rlow = R
                                        Rhigh = R + factor
                                        R= 0.5*(Rhigh+Rlow)
                                # if factor > aa*R :
                                #         factor = aa*R
                                
                                f = open('fileR.txt', 'a')
                                f.write("%s\n" % R)
                        else:
                                Not_Close = False
                                #print("got it")
                                R_final = R
                                Rguess  = R
                                Rold    = R
                                R_photo_final = R_photo
                                L_final = L_total
                                L_GCold = L_GC
                else:
                        if diff > prec :
                                Rhigh = R
                                R = 0.5*(Rhigh+Rlow)
                                # if factor > aa*R :
                                #         factor = aa*R
                                # R = R - factor
                                # f = open('fileR.txt', 'a')
                                # f.write("%s\n" % R)
                        elif diff < -prec :
                                factor = 0.5*factor
                                Rlow = R
                                R = 0.5*(Rhigh+Rlow)
                                # if factor > aa*R :
                                #         factor = aa*R
                                # R = R + factor
                                # f = open('fileR.txt', 'a')
                                # f.write("%s\n" % R)
                        else:
                                Not_Close = False
                                #print("got it")
                                R_final = R
                                Rguess  = R
                                Rold    = R
                                R_photo_final = R_photo
                                L_final = L_total
                                L_GCold = L_GC



        if True :
                print("luminosity = ",L_final/solarluminosity, "solar luminosity")
                print("luminosity from DM = ",L_DM/solarluminosity, "solar luminosity")
                # print("luminosity from capture = ",L_cap)
                print("luminosity from fusion = ",L_NR/solarluminosity, "solar luminosity")
                print("T_eff  = ", T_eff)
                print("Radius = ", R_photo_final/AU, "AU")
                print("DM mass fraction = ", mass_ratio)
                print("rhoDMc = ", rho_d[0]*rhoconvert)
                print("rhogc = ", rho_g[0]*rhoconvert)
                # print('C_cap  = ', C_cap  )
                # print('C_self = ', C_self )
                # print('C_ann  = ', C_ann  )
                # print('Regime = ', C_self/np.sqrt(C_cap*C_ann))
                # print('tau_EQ = ', tau_EQ/10**15)
                # print('tau_TH = ', tau_TH)
                # print('G_ann  = ', Gamma_ann)
                #exit()
        return L_DM/solarluminosity, L_final/solarluminosity, T_eff, R_photo_final/AU, rho_d[0]*rhoconvert, rho_g[0]*rhoconvert, mass_ratio

# plt.plot(heating_para,radius)
# plt.xlabel('Heating Parameter')
# plt.ylabel('Radius [AU]')
# plt.xscale('log')
# plt.show()

# plt.plot(heating_para,temp)
# plt.xlabel('Heating Parameter')
# plt.ylabel('Temperature [K]')
# plt.xscale('log')
# plt.show()

# plt.plot(heating_para,lum)
# plt.xlabel('Heating Parameter')
# plt.ylabel(r'$L_{total}~[L_{\odot}]$')
# plt.xscale('log')
# plt.show()

# plt.plot(heating_para,lum_perc)
# plt.xlabel('Heating Parameter')
# plt.ylabel(r'$L_{DM}/L_{total}$')
# plt.xscale('log')
# plt.show()

# heating_para = heating_para.reshape((-1,1))
# radius = radius.reshape((-1,1))
# temp = temp.reshape((-1,1))
# lum = lum.reshape((-1,1))
# lum_perc = lum_perc.reshape((-1,1))
# data = np.hstack((heating_para,radius,temp,lum,lum_perc))

# np.savetxt("DS_properties_500_solar.txt",data,delimiter='\t',header="heating_para\tradius[AU]\ttemp[K]\tlum[solar]\tlum_perc")
# final_data = [R_final, M_tot, T_eff, L_final]
# with open('STELLAR_DATA.txt', 'a') as f:
#      for item in final_data:
#          f.write("%s, " % item)
#      f.write('\n')

# calculate Eddington luminosity and compare it to current luminosity
# if current luminosity is greater than Eddington luminosity, stop the mass building
# L_Edd = 1.216*10**(31.0)*M/solarmass
# print("L_cap    = ", L_cap)
# print("L_fusion = ", L_NR)
# print("L_grav   = ", L_GC)
# print("L_final  = ", L_final)
# print("L_Edd    = ", L_Edd)
# if L_final < L_Edd:
#         total_time = total_time + DltTm/yearinsec
#         if printing :
#           print("total time = ",total_time," years")
#           print("current mass = ",M/solarmass," solar masses")
#           print("corressponding radius = ",R_final)
#           print("effective temperature = ",T_eff)
#           print("luminosity = ",L_final/solarluminosity," solar luminosity")
#           print("luminosity from DM = ",L_DM/solarluminosity," solar luminosity")
#           print("luminosity from capture = ",L_cap/solarluminosity," solar luminosity")
#           print("luminosity from fusion = ",L_NR/solarluminosity," solar luminosity")
#           print("Eddington luminosity = ",L_Edd/solarluminosity," solar luminosity")
#         M = M + DeltaM
# else:
#         print("life time of dark star phase =",total_time," years")
#         print("maximum mass reached, M = ",M/solarmass," solar mass")
#         print("corressponding radius = ",R_final)
#         print("current luminosity = ",L_final/solarluminosity," solar luminosity")
#         print("luminosity from DM = ",L_DM/solarluminosity," solar luminosity")
#         print("luminosity from capture = ",L_cap/solarluminosity," solar luminosity")
#         print("luminosity from fusion = ",L_NR/solarluminosity," solar luminosity")
#         print("Eddington luminosity = ",L_Edd/solarluminosity," solar luminosity")
#         big_enough = True

# exit()

# if os.path.exists("STELLAR_DATA.txt"):
#   os.remove("STELLAR_DATA.txt")
# open("STELLAR_DATA.txt","w+")
# timecount = 0
# while not big_enough:
#         R = Rguess
#         timecount = timecount + 1
#         if M <= starmin:
#                 n = poly[0]
#                 zn = ksi[0]
#                 theta = deri[0]
#                 rho_ratio = ratio[0]
#         elif M > starmin and M <=starmax:
#                 n = 1.5 + (3.0-1.5)*(M-starmin)/(starmax-starmin)
#                 if n <= 2.0:
#                         zn = ksi[0]+(ksi[1]-ksi[0])*(n-poly[0])/(poly[1]-poly[0])
#                         theta = deri[0]+(deri[1]-deri[0])*(n-poly[0])/(poly[1]-poly[0])
#                         rho_ratio = ratio[0]+(ratio[1]-ratio[0])*(n-poly[0])/(poly[1]-poly[0])
#                 elif n > 2.0 and n <=2.5:
#                         zn = ksi[1]+(ksi[2]-ksi[1])*(n-poly[1])/(poly[2]-poly[1])
#                         theta = deri[1]+(deri[2]-deri[1])*(n-poly[1])/(poly[2]-poly[1])
#                         rho_ratio = ratio[1]+(ratio[2]-ratio[1])*(n-poly[1])/(poly[2]-poly[1])
#                 else:
#                         zn = ksi[2]+(ksi[3]-ksi[2])*(n-poly[2])/(poly[3]-poly[2])
#                         theta = deri[2]+(deri[3]-deri[2])*(n-poly[2])/(poly[3]-poly[2])
#                         rho_ratio = ratio[2]+(ratio[3]-ratio[2])*(n-poly[2])/(poly[3]-poly[2])
#         else:
#                 n = poly[3]
#                 zn = ksi[3]
#                 theta = deri[3]
#                 rho_ratio = ratio[3]

#         # parameter indicating whether two luminosities are close enough
#         Not_Close = True
#         factor    = 1.0*10**(10)
#         firstpass = True
#         incR      = True
#         decR      = True

#         count = 0
#         while Not_Close:
#                 #R = Rguess
#                 rho_bar = 3.0*M/(pi4*R**3.0)
#                 rho_cg  = rho_ratio*rho_bar
#                 rho_cd  = fdm*rho_cg
#                 # rho_cd  = rho_halo
#                 # fdm     = rho_cd / rho_cg
#                 K_g = pi4*GN*(R/zn)**2*rho_cg**(1.0-1.0/n)/(n+1.0)
#                 K_d = GN*M/(3.0*R) # (m/s)**2 
#                 alpha_g = (zn/R)**2           # 1/length**2 unit
#                 alpha_d = pi4*GN*rho_cg/K_d   # 1/length**2 unit
#                 #if printing:
#                 #    print('rho_bar=',rho_bar)
#                 #    print('R=',R)
#                 #    print('M=',M)
#                 #    print('alpha_d=',alpha_g)
 
#                 # model to solve
#                 def model(z,r):
#                         rho_g = z[0]
#                         rho_d = z[1]
#                         phi_g = z[2]
#                         phi_d = z[3]
#                         drho_gdr = phi_g
#                         drho_ddr = phi_d
#                         dphi_gdr = -alpha_g*(rho_g**n + fdm*rho_d) - 2.0*phi_g/r
#                         dphi_ddr = -alpha_d*(rho_g**n + fdm*rho_d)*rho_d - 2.0*phi_d/r + phi_d**2/rho_d
#                         dzdr = [drho_gdr,drho_ddr,dphi_gdr,dphi_ddr]
#                         return dzdr

#                 # initial condition
#                 z0 = [1,1,0,0]

#                 # number of radius points
#                 n_p = 15001

#                 # radius points
#                 #r = np.linspace(10,R,n_p)
#                 #delta_r = r[1]-r[0]

#                 n_p1=4000
#                 n_p2=4001
#                 n_p = n_p1+n_p2
#                 #r = np.logspace(1,np.log10(R),n_p)
#                 r1 = np.logspace(1,np.log10(0.9*R),n_p1)
#                 r2 = np.logspace(np.log10(0.91*R),np.log10(R),n_p2)
#                 r = [*r1, *r2]
#                 #with open('file_r.txt', 'w') as f:
#                 #  for item in r:
#                 #       f.write("%s\n" % item)
#                 r = np.array(r)
#                 delta_r = np.zeros(r.size)
#                 delta_r[0] = r[0]
#                 for i in range(0,n_p-1):
#                       delta_r[i+1] = r[i+1]-r[i]
#                 #print("r=",r)
#                 #print("delta_r=",delta_r)

#                 # store solutions
#                 rho_g = np.zeros(r.size)
#                 rho_d = np.zeros(r.size)
#                 phi_g = np.zeros(r.size)
#                 phi_d = np.zeros(r.size)

#                 # initial conditions
#                 rho_g[0] = z0[0]
#                 rho_d[0] = z0[1]
#                 phi_g[0] = z0[2]
#                 phi_d[0] = z0[3]

#                 # solve our group of ODE
#                 for i in range(1,n_p):
#                         r_span = [r[i-1],r[i]]
#                         # integrate to next step
#                         z = odeint(model,z0,r_span)
#                         # store solution
#                         rho_g[i] = z[1][0]
#                         rho_d[i] = z[1][1]
#                         phi_g[i] = z[1][2]
#                         phi_d[i] = z[1][3]
#                         # next initial condition
#                         z0 = z[1]

#                 # eliminate all the nan's in both profiles
#                 rho_g = rho_cg*rho_g[~np.isnan(rho_g)]**n
#                 rho_d = rho_cd*rho_d[~np.isnan(rho_d)]

#                 # find the mass profile (useful in finding photosphere)
#                 stellar_mass = np.zeros(rho_g.size)
#                 stellar_mass[0] = pi4/3.0*r[0]**3*rho_g[0]
#                 for i in range(rho_g.size-1):
#                    stellar_mass[i+1] = (stellar_mass[i] + \
#                                        pi4/3.0*(r[i+1]**3-r[i]**3)*rho_g[i+1])

#                 # integrated number of DM particles
#                 NDM = np.zeros(rho_d.size)
#                 NDM[0] = pi4/3.0*r[0]**3*rho_d[0]
#                 for i in range(rho_d.size-1):
#                     NDM[i+1] = NDM[i] + pi4/3.0*(r[i+1]**3-r[i]**3)*rho_d[i+1]
#                 NDM = NDM / m_DM

#                 # define the integrated luminosity
#                 totlumin = np.zeros(rho_g.size)

#                 # find the temperature profile and pressure profile
#                 T = np.zeros(rho_g.size)
#                 P = K_g*rho_g**(1.0+1.0/n)
#                 for i in range(rho_g.size):
#                         T[i] = tsub(rho_g[i],P[i],P[i]*mu*m_u/rho_g[i]/kB)

#                 with open('file_r.txt', 'w') as f:
#                   for item in r:
#                        f.write("%s\n" % item)
#                 with open('file_t.txt', 'w') as f:
#                   for item in T:
#                        f.write("%s\n" % item)
#                 with open('file_g.txt', 'w') as f:
#                   for item in rho_g:
#                        f.write("%s\n" % item)
#                 with open('file_d.txt', 'w') as f:
#                   for item in rho_d:
#                        f.write("%s\n" % item)
 
#                 # find the position of photosphere
#                 photo = np.zeros(rho_g.size)
#                 for i in range(rho_g.size):
#                         photo[i] = Opacity.opac(np.log10(rho_g[i]),T[i])*P[i]-2.0/3.0*GN*M/R**2

#                 R_photo = 0.0
#                 T_eff = 0.0
#                 for i in range(rho_g.size-1):
#                         if photo[i]>=0 and photo[i+1]<=0:
#                                 R_photo = r[i]
#                                 T_eff   = T[i]
#                                 M_tot   = stellar_mass[i]
#                                 im      = i
#                                 exit
#                 for i in range(rho_d.size-1):
#                         if R_photo <= r[i]:
#                                 N_DM    = NDM[i]
#                                 exit

#                 # This part is about capture and self-capture, which is defuted by other researchers like Haibo et al.
#                 # C_ann = rho_d[0]**2*pi4/3.0*r[0]**3
#                 # for i in range(rho_d.size-1):
#                 #      C_ann = C_ann + pi4/3.0*(r[i+1]**3-r[i]**3)*rho_d[i+1]**2
#                 # C_ann = 2.0*sigma_v*C_ann/(m_DM*N_DM)**2

#                 # # Compute the capture rate
#                 # B = 1.5*GN*stellar_mass[0]/(r[0]*vbar**2)
#                 # coeff = B - 1 + np.exp(-B)
#                 # C_cap = rho_g[0]*coeff*r[0]**3
#                 # for i in range(rho_g.size-1):
#                 #      B     = 1.5*GN*stellar_mass[i+1]/(r[i+1]*vbar**2)
#                 #      coeff = B - 1 + np.exp(-B)
#                 #      C_cap = C_cap + rho_g[i+1]*coeff*(r[i+1]**3-r[i]**3)
#                 # C_cap = C_cap*2.0*rho_halo*X/m_H/m_DM*sigma_SD*vbar/np.sqrt(np.pi)*pi4/3.0

#                 # # Compute the rate of self-capture
#                 # B = 1.5*GN*stellar_mass[0]/(r[0]*vbar**2)
#                 # coeff = B - 1 + np.exp(-B)
#                 # C_self = rho_d[0]*coeff*r[0]**3
#                 # for i in range(rho_d.size-1):
#                 #      B      = 1.5*GN*stellar_mass[i+1]/(r[i+1]*vbar**2)
#                 #      coeff  = B - 1 + np.exp(-B)
#                 #      C_self = C_self + rho_d[i+1]*coeff*(r[i+1]**3-r[i]**3)
#                 # C_self = 2.0*rho_halo*vbar/np.sqrt(np.pi)*pi4/3.0*C_self*sigma_ss/m_DM**2/N_DM

#                 # tau_EQ    = 1.0/np.sqrt((C_self/2.0)**2 + C_cap*C_ann)
#                 # tnh       = np.tanh(DltTm/tau_EQ)
#                 # if DltTm < 20.0*tau_EQ :
#                 #      N_cap = C_cap * tnh/ (1/tau_EQ - (C_self/2.0)*tnh)
#                 # else :
#                 #      N_cap = C_cap * tnh/ (C_cap*C_ann/C_self + C_self/2.0*(1.0-tnh))
#                 # Gamma_ann = 0.5*C_ann*N_cap**2
#                 # tau_TH    = m_DM/(rho_d[0]*sigma_ss*vbar)

#                 # if True :
#                 #   print('C_cap  = ',C_cap  )
#                 #   print('C_self = ',C_self )
#                 #   print('C_ann  = ',C_ann  )
#                 #   print('Regime = ', C_self/np.sqrt(C_cap*C_ann))
#                 #   print('tau_EQ = ',tau_EQ )
#                 #   print('time   = ',DltTm )
#                 #   print('N_cap  = ',N_cap)
#                 #   print('N_DM   = ',N_DM)
#                 #   print('G_ann  = ',Gamma_ann)
                

#                 # calculate effective luminosity
#                 L_eff = pi4*R_photo**2.0*sigma*T_eff**4.0

#                 # include luminosity contributed from different sources:
#                 # DM annihilation
#                 L_DM = rho_d[0]**2*pi4*r[0]**3
#                 for i in range(rho_d.size):
#                       if r[i] <= R_photo :
#                              L_DM = L_DM + rho_d[i]**2*pi4*(r[i+1]**3-r[i]**3)
#                 L_DM = sigma_v*L_DM*c**2/m_DM

#                 #print('sigma_v = ', sigma_v)
#                 #print('c = ', c)
#                 #print('m_DM = ', m_DM)
#                 #print('rho_DM = ', rho_d[0])
#                 #print('L_DM = ', L_DM)

#                 # nuclear reaction
#                 L_NR = 0.0
#                 for i in range(rho_d.size):
#                         rhog  = 10**(-3)*rho_g[i] #convert the density from kg/m^3 to g/cm^3 
#                         num_p = rhog*X/m_u
#                         num_d = num_p/7000.0
#                         num_He3 = 0.0
#                         nucr = nucrat.nucrat(rhog,T[i],num_p,num_d,num_He3,DltTm,0.24,0.0,0.0,0.0)
#                         L_NR = L_NR + nucr*10.0**(-4.0)*pi4*r[i]**2.0*delta_r[i]*rho_g[i]
 
#                 # gravitational contraction
#                 L_GC  = 0.0
#                 #for i in range(rho_g.size):
#                 #      if r[i] <= R_photo :
#                 #            L_GC = L_GC + (pi4*rho_g[i]*r[i]**2)**2*delta_r[i]*GN/3.0
#                 #L_GC = (L_GC-L_GCold)/DltTm
#                 #L_GC = 0.6*GN * M * dotM/R_photo/yearinsec

#                 # Dark matter capture heating
#                 # L_cap = f_EM*Gamma_ann*2.0*m_DM*c**2
#                 L_cap = 0.0

#                 L_total = L_DM + L_NR + L_GC + L_cap
#                 diff    = (L_total-L_eff)/(L_total+L_eff)

#                 count = count + 1
#                 if True:
#                   print("R=",R)
#                   print("R_photo=",R_photo)
#                   print("T_eff=",T_eff)
#                   print("L_eff=",L_eff)
#                   print("L_DM=",L_DM)
#                   print("L_cap=",L_cap)
#                   print('L_GC = ',L_GC)
#                   print('L_NR = ',L_NR)
#                   print("L_total=",L_total)
#                   print('ratio=',diff)
#                   print('count=',count)
#                   #exit()

#                 aa = 0.05
#                 diff = -2*diff
#                 if firstpass:
#                     if diff > prec :
#                       if decR : 
#                         incR = False
#                       else:
#                         firstpass = False
#                         factor = 0.5*factor
#                       if factor > aa*R :
#                         factor = aa*R
#                       R = R - factor

#                       f = open('fileR.txt', 'a')
#                       f.write("%s\n" % R)
#                     elif diff < -prec :
#                       print("R should increase")
#                       if incR :
#                         decR= False
#                       else:
#                         firstpass = False
#                         factor = 0.5*factor
#                       if factor > aa*R :
#                         factor = aa*R
#                       R = R + factor
#                       f = open('fileR.txt', 'a')
#                       f.write("%s\n" % R)
#                     else:
#                       Not_Close = False
#                       #print("got it")
#                       R_final = R
#                       Rguess  = R
#                       Rold    = R
#                       R_photo_final = R_photo
#                       L_final = L_total
#                       L_GCold = L_GC
#                 else:
#                     if diff > prec :
#                       factor = 0.5*factor
#                       if factor > aa*R :
#                         factor = aa*R
#                       R = R - factor
#                       f = open('fileR.txt', 'a')
#                       f.write("%s\n" % R)
#                     elif diff < -prec :
#                       factor = 0.5*factor
#                       if factor > aa*R :
#                         factor = aa*R
#                       R = R + factor
#                       f = open('fileR.txt', 'a')
#                       f.write("%s\n" % R)
#                     else:
#                       Not_Close = False
#                       #print("got it")
#                       R_final = R
#                       Rguess  = R
#                       Rold    = R
#                       R_photo_final = R_photo
#                       L_final = L_total
#                       L_GCold = L_GC

#         if True :
#                   print("luminosity = ",L_final)
#                   print("luminosity from DM = ",L_DM)
#                   # print("luminosity from capture = ",L_cap)
#                   print("luminosity from fusion = ",L_NR)
#                   print("T_eff  = ", T_eff)
#                   print("Radius = ", R_photo_final/10**11)
#                   # print('C_cap  = ', C_cap  )
#                   # print('C_self = ', C_self )
#                   # print('C_ann  = ', C_ann  )
#                   # print('Regime = ', C_self/np.sqrt(C_cap*C_ann))
#                   # print('tau_EQ = ', tau_EQ/10**15)
#                   # print('tau_TH = ', tau_TH)
#                   # print('G_ann  = ', Gamma_ann)
#                   #exit()


#         final_data = [R_final, M_tot, T_eff, L_final]
#         with open('STELLAR_DATA.txt', 'a') as f:
#              for item in final_data:
#                  f.write("%s, " % item)
#              f.write('\n')

#         # calculate Eddington luminosity and compare it to current luminosity
#         # if current luminosity is greater than Eddington luminosity, stop the mass building
#         L_Edd = 1.216*10**(31.0)*M/solarmass
#         print("L_cap    = ", L_cap)
#         print("L_fusion = ", L_NR)
#         print("L_grav   = ", L_GC)
#         print("L_final  = ", L_final)
#         print("L_Edd    = ", L_Edd)
#         if L_final < L_Edd:
#                 total_time = total_time + DltTm/yearinsec
#                 if printing :
#                   print("total time = ",total_time," years")
#                   print("current mass = ",M/solarmass," solar masses")
#                   print("corressponding radius = ",R_final)
#                   print("effective temperature = ",T_eff)
#                   print("luminosity = ",L_final/solarluminosity," solar luminosity")
#                   print("luminosity from DM = ",L_DM/solarluminosity," solar luminosity")
#                   print("luminosity from capture = ",L_cap/solarluminosity," solar luminosity")
#                   print("luminosity from fusion = ",L_NR/solarluminosity," solar luminosity")
#                   print("Eddington luminosity = ",L_Edd/solarluminosity," solar luminosity")
#                 M = M + DeltaM
#         else:
#                 print("life time of dark star phase =",total_time," years")
#                 print("maximum mass reached, M = ",M/solarmass," solar mass")
#                 print("corressponding radius = ",R_final)
#                 print("current luminosity = ",L_final/solarluminosity," solar luminosity")
#                 print("luminosity from DM = ",L_DM/solarluminosity," solar luminosity")
#                 print("luminosity from capture = ",L_cap/solarluminosity," solar luminosity")
#                 print("luminosity from fusion = ",L_NR/solarluminosity," solar luminosity")
#                 print("Eddington luminosity = ",L_Edd/solarluminosity," solar luminosity")
#                 big_enough = True

# exit()