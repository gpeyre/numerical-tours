import numpy as np

def compute_wavelet_filter(type,par):
    """
        compute_wavelet_filter - Generate Orthonormal QMF Filter for Wavelet Transform
        
        
           [h,g] = compute_wavelet_filter(Type,Par)
        
         Inputs
           Type   string, 'Haar', 'Beylkin', 'Coiflet', 'Daubechies',
                  'Symmlet', 'Vaidyanathan','Battle'
           Par    integer, it is a parameter related to the support and vanishing
                  moments of the wavelets, explained below for each wavelet.
        
        Outputs
          h   low pass quadrature mirror filter
          g   high pass
        
         Description
           The Haar filter (which could be considered a Daubechies-2) was the
           first wavelet, though not called as such, and is discontinuous.
    
           The Beylkin filter places roots for the frequency response function
           close to the Nyquist frequency on the real axis.
         
           The Coiflet filters are designed to give both the mother and father
           wavelets 2*Par vanishing moments; here Par may be one of 1,2,3,4 or 5.
         
           The Daubechies filters are minimal phase filters that generate wavelets
           which have a minimal support for a given number of vanishing moments.
           They are indexed by their length, Par, which may be one of
           2,4,6,8,10,12,14,16,18 or 20. The number of vanishing moments is par/2.
         
           Symmlets are also wavelets within a minimum size support for a given
           number of vanishing moments, but they are as symmetrical as possible,
           as opposed to the Daubechies filters which are highly asymmetrical.
           They are indexed by Par, which specifies the number of vanishing
           moments and is equal to half the size of the support. It ranges
           from 4 to 10.
         
           The Vaidyanathan filter gives an exact reconstruction, but does not
           satisfy any moment condition.  The filter has been optimized for
           speech coding.
         
           The Battle-Lemarie filter generate spline orthogonal wavelet basis.
           The parameter Par gives the degree of the spline. The number of
           vanishing moments is Par+1.
         
        See Also
           FWT_PO, IWT_PO, FWT2_PO, IWT2_PO, WPAnalysis
    
        References
            The books by Daubechies and Wickerhauser.
            
        Warning : only Daubechies implemented for the moment !
    """
     
    if type == 'Daubechies':
        
        if par == 1:
            f = [1,1]/np.sqrt(2)

        if par == 4:
            f = [.482962913145,.836516303738,
                .224143868042,-.129409522551]
        
        if par == 6:
            f = [.332670552950,.806891509311,
            .459877502118,-.135011020010,
                -.085441273882,.035226291882]
        
        if par == 8:
            f = [ .230377813309,.714846570553,
                .630880767930,-.027983769417,
                -.187034811719,.030841381836,
                .032883011667,-.010597401785]
        
        if par == 10:
            f = [.160102397974,.603829269797,.724308528438,
                .138428145901,-.242294887066,-.032244869585,
                .077571493840,-.006241490213,-.012580751999,
                .003335725285]
        
        if par == 12:
            f = [.111540743350,.494623890398,.751133908021,
                .315250351709,-.226264693965,-.129766867567,
                .097501605587,.027522865530,-.031582039317,
                .000553842201,.004777257511,-.001077301085]
        
        if par == 14:
            f = [.077852054085,.396539319482,.729132090846,
                .469782287405,-.143906003929,-.224036184994,
                .071309219267,.080612609151,-.038029936935,
                -.016574541631,.012550998556,.000429577973,
                -.001801640704,.000353713800]
        
        if par == 16:
            f = [.054415842243,.312871590914,.675630736297,
                .585354683654,-.015829105256,-.284015542962,
                .000472484574,.128747426620,-.017369301002,
                -.044088253931,.013981027917,.008746094047,
                -.004870352993,-.000391740373,.000675449406,
                -.000117476784]
        
        if par == 18:
            f = [.038077947364,.243834674613,.604823123690,
                .657288078051,.133197385825,-.293273783279,
                -.096840783223,.148540749338,.030725681479,
                -.067632829061,.000250947115,.022361662124,
                -.004723204758,-.004281503682,.001847646883,
                .000230385764,-.000251963189,.000039347320]
        
        if par == 20:
            f = [.026670057901,.188176800078,.527201188932,
                .688459039454,.281172343661,-.249846424327,
                -.195946274377,.127369340336,.093057364604,
                -.071394147166,-.029457536822,.033212674059,
                .003606553567,-.010733175483,.001395351747,
                .001992405295,-.000685856695,-.000116466855,
                .000093588670,-.000013264203]
    
    else:
        raise ValueError("Wrong arguments, see comments for acceptable values")
        
    f = list(f/np.linalg.norm(f))
    
    if len(f)%2 == 0:
        f = [0] + f
    
    return f
