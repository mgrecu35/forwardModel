from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np


def read_wrf(fname,it):
    f=Dataset(fname)
    qv=f['QVAPOR'][it,:,:,:]     # water vapor
    qr=f['QRAIN'][it,:,:,:]      # rain mixing ratio
    qs=f['QICE'][it,:,:,:]       # snow mixing ratio
    qc=f['QCLOUD'][it,:,:,:]     # cloud mixing ratio
    qg=f['QGRAUP'][it,:,:,:]     # graupel mixing ratio
    ncr=f['QNRAIN'][it,:,:,:]    # rain mixing ratio
    ncs=f['QNICE'][it,:,:,:]     # snow mixing ratio
    ncg=f['QNGRAUPEL'][it,:,:,:] # graupel mixing ratio
    
    th=f['T'][it,:,:,:]+300    # potential temperature (K)
    prs=f['P'][it,:,:,:]+f['PB'][it,:,:,:]  # pressure (Pa)
    T=th*(prs/100000)**0.286  # Temperature
    t2c=T-273.15
    #stop
    z=(f['PHB'][it,:,:,:]+f['PH'][it,:,:,:])/9.81/1000.
    xlat=f['XLAT'][0,:,:]
    xlong=f['XLONG'][0,:,:]
    R=287.058  #J*kg-1*K-1
    rho=prs/(R*T)
    return qr,qs,qg,ncr,ncs,ncg,rho,z,T,prs,f,qv
it=0
fname="../DA_MCS/wrfout_d03_2017-07-04_00:00:00"
it=0

import combAlg as sdsu
sdsu.mainfortpy()
sdsu.initp2()

def process_time_step(fname,it):
    qr,qs,qg,ncr,ncs,ncg,rho,z,T,prs,fh,qv=read_wrf(fname,it)
    nwL=[]
    dmL=[]
    rwc=qr*rho*1e3
    swc=qs*rho*1e3
    gwc=qg*rho*1e3
    wv=qv*rho
    hint=np.arange(120)*0.125+0.125/2
    
    sl=2.0
    dr=0.125
    pRateL=[]
    nz,ny,nx=rwc.shape
    zKuL=[]
    zKaL=[]
    piaKaL=[]
    piaKuL=[]
    for i in range(nx):
        for j in range(ny):
            dn1=np.zeros((120),float)
            if rwc[0,j,i]>0.05:
                hm=0.5*(z[1:,j,i]+z[:-1,j,i])
                rwc1=np.interp(hint,hm,rwc[:,j,i])
                swc1=np.interp(hint,hm,swc[:,j,i])
                gwc1=np.interp(hint,hm,gwc[:,j,i])
                temp=np.interp(hint,hm,T[:,j,i])
                press=np.interp(hint,hm,prs[:,j,i])
                wv1=np.interp(hint,hm,wv[:,j,i])
                swc1+=gwc1
                zku_m ,zku_t, attku, piaku, \
                    kext,salb,asym,kext_,salb_,asym_,pRate,\
                    dm_out=sdsu.reflectivity_ku(rwc1,swc1,wv1,dn1,temp,press,dr,sl)
                a=np.nonzero(rwc1[:36]>0.1)
                dn1[a]+=-3.19/2*np.log10(dm_out[a]/1.2)
                #print(dm_out[a])
                dn1+=np.random.randn()*0.75
                zku_m ,zku_t, attku, piaku, \
                    kext,salb,asym,kext_,salb_,asym_,pRate,\
                    dm_out=sdsu.reflectivity_ku(rwc1,swc1,wv1,dn1,temp,press,dr,sl)
                nwL.extend(dn1[a])
                dmL.extend(dm_out[a])
                #print(dn1[a])
                #print(dm_out[a])
                zka_m ,zka_t, attka, piaka, \
                    kext_ka,salb_ka,asym_ka,kext_ka_,salb_ka_,asym_ka_,pRateKa\
                    =sdsu.reflectivity_ka(rwc1,swc1,wv1,dn1,temp,press,dr,sl)
                
                noms=0
                alt=400.
                freq=13.8
                nonorm=0
                theta=0.35/2.0
                freqKa=35.5
                zms = sdsu.multiscatterf(kext[::-1],salb[::-1],asym[::-1],\
                                         zku_t[::-1],dr,noms,alt,\
                                         theta,freq,nonorm)
                zms_ka = sdsu.multiscatterf(kext_ka[::-1],salb_ka[::-1],\
                                            asym_ka[::-1],\
                                            zka_t[::-1],dr,noms,alt,\
                                            theta,freqKa,nonorm)
                zKuL.append(zku_m)
                zKaL.append(zms_ka)
                piaKaL.append(piaka)
                piaKuL.append(piaku)
                pRateL.append(pRate)
    return zKuL,zKaL,pRateL,piaKuL,piaKaL,dmL,nwL
                


it=0
import xarray as xr

for it in range(0,5):
    zKuL,zKaL,pRateL,piaKuL,piaKaL=[],[],[],[],[]
    print("it=%2.2i"%it)
    for iseed in range(3):
        zKuL1,zKaL1,pRateL1,piaKuL1,piaKaL1,dmL,nwL=process_time_step(fname,it)
        zKuL.extend(zKuL1)
        zKaL.extend(zKaL1)
        pRateL.extend(pRateL1)
        piaKuL.extend(piaKuL1)
        piaKaL.extend(piaKaL1)
        
    zKu=xr.DataArray(zKuL)
    zKa=xr.DataArray(zKaL)
    pRate=xr.DataArray(pRateL)
    piaKu=xr.DataArray(piaKuL)
    piaKa=xr.DataArray(piaKaL)
    #stop
    ds=xr.Dataset({"zKu":zKu,"zKa":zKa,"pRate":pRate,"piaKu":piaKu,"piaKa":piaKa})
    ds.to_netcdf("%s.%2.2i.nc"%(fname[10:],it))
    
