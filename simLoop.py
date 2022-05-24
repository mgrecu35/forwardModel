import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter

def simLoop(f,zKuL,zKu_msL,xL,beta,maxHL,zKu_cL,sdsu,zku_3d,
            zku_3d_ms, zku_3d_true,\
            pRateL,attKuL,rwcL,swcL,f1L,f2L,ijL,dnL,
            zKaL,zKa_msL,zKa_cL,attKaL,sl):
    R=287
    h=(250/2.+np.arange(76)*250)/1e3
    h1=(0.+np.arange(77)*250)/1e3
    fh=Dataset(f)
    z=fh['zh'][:]
    qr=fh["qr"][0,:,:,:]
    qg=fh["qg"][0,:,:,:]
    qs=fh["qs"][0,:,:,:]
    a=np.nonzero(qr[0,:,:]>-0.00015)
    print(a[0].shape[0],qr[0,:,:].max())
    T=fh["th"][0,:,:,:]
    prs=fh["prs"][0,:,:,:]
    T=T*(prs/1e5)**(0.286)
    rho=prs/(R*T)
    rwc=qr*rho*1e3
    gwc=(qg+qs)*rho*1e3
    qv=fh['qv'][0,:,:,:]
    wv=qv*rho
    zku_3d=zku_3d*0.0-99
    zku_3d_ms=zku_3d_ms*0.0-99
    zku_3d_true=zku_3d_true*0.0-99
    zka_3d=zku_3d*0.0-99
    zka_3d_ms=zku_3d_ms*0.0-99
    zka_3d_true=zku_3d_true*0.0-99
    pRate3D=zku_3d_true*0.0
    #plt.figure(figsize=(8,11))
    nz,nx,ny=rwc.shape
    pia2d=np.zeros((nx,ny),float)
    pia2d_ka=np.zeros((nx,ny),float)
    f1_=np.exp(np.random.randn(76,200,200)*1.5)
    f2_=np.exp(np.random.randn(76,200,200)*1.5)
    f1_=gaussian_filter(f1_, sigma=2)
    f2_=gaussian_filter(f2_, sigma=2)
    dn3d=np.zeros((76,200,200),float)+np.random.randn()*0.75
    dn3d=gaussian_filter(dn3d, sigma=2)
    for i1,i2 in zip(a[0],a[1]):
        zKu_1D=[]
        att_1D=[]
        pwc_1D=[]
        prate_1D=[]
        f1=f1_[:,i1,i2]
        f2=f2_[:,i1,i2]
        piaKu=0
        dr=0.250
        dn1=0*dn3d[:,i1,i2]-0.05
        rwc1=np.interp(h,z[:],rwc[:,i1,i2])*1
        swc1=np.interp(h,z[:],gwc[:,i1,i2])*1
        temp=np.interp(h,z[:],T[:,i1,i2])
        temp1=np.interp(h1,z[:],T[:,i1,i2])
        press=np.interp(h,z[:],prs[:,i1,i2])
        wv1=np.interp(h,z[:],wv[:,i1,i2])
        a=np.nonzero(swc1>0.01)
        if len(a[0])>5:
            ht1=h[a[0][-1]]
            ht2=ht1+np.random.random()*1
            ht2=min(h[-1],ht2)
            hb1=h[a[0][0]]
            dh=(ht2-ht1)/len(a[0]-1)*np.arange(len(a[0]))
            hint=h[a[0]]+dh
            ntop=min(int(ht2/dr)+1,75)
            swc11=np.interp(h[a[0][0]:ntop],h[a[0]]+dh,swc1[a[0]])
            #for k in range(a[0][0],ntop):
            #    swc1[k]=max(swc1[k],swc11[k-a[0][0]])

        rwc1=(0.4+0.6*f1)*rwc1
        swc1=(0.4+0.6*f2)*swc1
        twc=rwc1
        a1=np.nonzero(twc>0.1)
        dn1[a1]=dn1[a1]-0.1*np.log10(twc[a1]/0.1)
        #if 16 in a1[0]:
        #    dn1[17:20]+=-0.1*np.log10(twc[16]/0.1)
        #sl=2.0
        zku_m ,zku_t, attku, piaku, \
            kext,salb,asym,kext_,salb_,asym_,pRate\
            =sdsu.reflectivity_ku(rwc1,swc1,wv1,dn1,temp,press,dr,sl)
        zka_m ,zka_t, attka, piaka, \
            kext_ka,salb_ka,asym_ka,kext_ka_,salb_ka_,asym_ka_,pRateKa\
            =sdsu.reflectivity_ka(rwc1,swc1,wv1,dn1,temp,press,dr,sl)
        dr=0.25
        noms=0
        alt=400.
        freq=13.8
        nonorm=0
        theta=0.35/3.0
        zku_3d_true[:,i1,i2]=zku_t.copy()
        zka_3d_true[:,i1,i2]=zka_t.copy()
        zku_3d[:,i1,i2]=zku_m.copy()
        zku_3d_ms[:,i1,i2]=zku_m.copy()
        zka_3d[:,i1,i2]=zka_m.copy()
        zka_3d_ms[:,i1,i2]=zka_m.copy()
        pia2d[i1,i2]=piaku
        pia2d_ka[i1,i2]=piaka
        freqKa=35.5
        pRate3D[:,i1,i2]=pRate
        if zku_m.max()>25 and i1>2 and i1<nx-2 and i2>2 and i2<ny-2:
            ijL.append([i1,i2])
            zms = sdsu.multiscatterf(kext[::-1],salb[::-1],asym[::-1],\
                                 zku_t[::-1],dr,noms,alt,\
                                 theta,freq,nonorm)
            zms_ka = sdsu.multiscatterf(kext_ka[::-1],salb_ka[::-1],\
                                        asym_ka[::-1],\
                                        zka_t[::-1],dr,noms,alt,\
                                        theta,freqKa,nonorm)
            
            pRateL.append(pRate)
            zKuL.append(zku_m)
            zKaL.append(zka_m)
            f1L.append(f1)
            f2L.append(f2)
            dnL.append(dn1)
            zKu_msL.append(zms[::-1])
            zKa_msL.append(zms_ka[::-1])
            zku_3d_ms[:,i1,i2]=zms[::-1].copy()
            zka_3d_ms[:,i1,i2]=zms_ka[::-1].copy()
            attKuL.append(attku)
            attKaL.append(attka)
            rwcL.append(rwc1)
            swcL.append(swc1)
            zKu_cL.append(zku_t)
            zKa_cL.append(zka_t)
            ind=np.argmax(zku_m[12:24])
            maxHL.append([zku_m[12+ind],ind,zku_t[12+ind]-zku_m[12+ind]])
            

        #zms = sdsu.multiscatterf(kext[::-1],salb[::-1],asym[::-1],\
        #                         zku_t[::-1],dr,noms,alt,\
        #                         theta,freq,nonorm)
        #kext1=np.zeros((75),float)
        #salb1=np.zeros((75),float)
        #asym1=np.zeros((75),float)
        #for k in range(75):
        #    kext1[k]=kext[2*k:2*k+2].mean()
        #    salb1[k]=salb[2*k:2*k+2].mean()
        #    asym1[k]=asym[2*k:2*k+2].mean()

        #kext1_=np.zeros((75),float)
        #salb1_=np.zeros((75),float)
        #asym1_=np.zeros((75),float)
        #for k in range(75):
        #    kext1_[k]=kext_[2*k:2*k+2].mean()
        #    salb1_[k]=salb_[2*k:2*k+2].mean()
        #    asym1_[k]=asym_[2*k:2*k+2].mean()
        #emis=0.85+0.1*np.random.rand()
        #ebar=emis
        #lambert=1
        #salb1[salb1>0.99]=0.99
        #tb = sdsu.radtran(umu,temp1[0],temp1,h1/1000.,kext1,salb1,asym1,\
        #                  fisot,emis,ebar,lambert)
    return zku_3d,zku_3d_ms,pia2d,zka_3d,zka_3d_ms,pia2d_ka,zku_3d_true,zka_3d_true,pRate3D
