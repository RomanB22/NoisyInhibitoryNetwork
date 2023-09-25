import numpy as np

def compute_SPC_R2_meanFR(
    spbins,      #STH 1ms/bin
    ksigma= 3.,  #sigma for Gaussian to smooth FR
    kwidth=25.,  #kernel width
    peaks=False  #set True if peaks are neede
    ):
    """
    The function get network spike-time histogram as numpy array
    with 1ms bin, smooth it with Gaussian kernel, find peaks (peak detector)
    and computes participation (spikes-per-cycle, spc), R2, and mean network period,
    Returns tuple:
        spc,R2,mean_Tnet,[peakmarks]
    
    Ruben Tikidji-Hamburyan (LSUHSC, GWU, 2013-2021)
    """
    print("==================================")
    print("===         Peak Detector      ===")
    kernel = np.arange(-kwidth,kwidth,1.)
    kernel = np.exp(kernel**2/(-ksigma**2))
    module = np.convolve(spbins,kernel,mode='same')
    
    spbinmark = []
    for idx in (np.diff(np.sign(np.diff(module))) < 0).nonzero()[0] + 1:
        spbinmark.append([idx,1])
    for idx in (np.diff(np.sign(np.diff(module))) > 0).nonzero()[0] + 1:
        spbinmark.append([idx,-1])
    peakmark  = []
    spc,ccnt = 0.,0.
    if len(spbinmark) > 2:
        spbinmark.sort()
        spbinmark = np.array(spbinmark)
        for mx in np.where( spbinmark[:,1] > 0 )[0]:
            if mx <= 2 or mx >= (spbinmark.shape[0] -2):continue
            if spbinmark[mx-1][1] > 0 or spbinmark[mx+1][1] > 0 or spbinmark[mx][1] < 0:continue
            peakmark.append(spbinmark[mx].tolist()[0]-1)
            ccnt    += 1
            curespc = np.sum(spbins[spbinmark[mx-1][0]:spbinmark[mx+1][0]])
            spc += curespc
    else:
        spbinmark = None
        
    print("==================================")
    print("===             R2             ===")

    if ccnt > 0:
        spc = spc /ccnt
    else:
        spc = None
    X,Y,Rcnt,netpermean,netpercnt=0.,0.,0.,0.,0.
    if not spbinmark is None:
        for mx in np.where( spbinmark[:,1] > 0 )[0]:
            #if mx >= (spbinmark.shape[0]/2 - 3):continue
            if mx >= (spbinmark.shape[0] - 3):continue
            if spbinmark[mx+1][1] > 0 or spbinmark[mx+2][1] < 0 or spbinmark[mx][1] < 0:continue
            Pnet = float(spbinmark[mx+2][0] - spbinmark[mx][0])
            netpermean += Pnet
            netpercnt  += 1.
            for i,n in enumerate(spbins[spbinmark[mx][0]:spbinmark[mx+2][0]]):
                phyX = np.cos(np.pi*2.*float(i)/Pnet)
                phyY = np.sin(np.pi*2.*float(i)/Pnet)
                X += n*phyX
                Y += n*phyY
                Rcnt += n
    if Rcnt > 0.:
        R2 = (X/Rcnt)**2+(Y/Rcnt)**2
    else:
        R2 = None
    if netpercnt > 1.:
        netpermean /= ( netpercnt - 1)
    else:
        netpermean = None
    print("  > R2       =           ",R2)
    print("  > SPC      =           ",spc)
    print("  > netPmean =           ",netpermean)
    print("==================================\n")
    return (spc,R2,netpermean,peakmark) if peaks else (spc,R2,netpermean)

if __name__ == "__main__":
    from matplotlib.pyplot import *
    sth = np.zeros(2000)
    nt  = np.random.exponential(1000/30)
    na  = np.random.randn()*30 + 300
    while nt < 2000:
        p   = np.arange(2000)
        p   = np.exp(-(p-nt)**2/900)*na#+np.random.rand(2000)*30
        p[np.where(p<0)] = 0
        sth+= p
        nt += np.random.exponential(1000/30)
        na  = np.random.randn()*30 + 300
    p   = np.arange(2000)
    plot(p,sth,'k-')
    spc,R2,netpermean,peakmark = compute_SPC_R2_meanFR(sth,peaks=True)
    for p in peakmark:
        a = sth[int(p)]
        plot([p,p],[a-30,a+30],'r--')
    
    show()

