import numpy as np

#### Uncertainty
def sigma_err(sigma_data, oversample):

    '''Calculate uncertainy for bin considering only y data'''
    id_nonnan = ~np.isnan(sigma_data)
    sigma_data_ = sigma_data[id_nonnan]
    n = len(sigma_data_)

    sigma_err_var = np.sqrt(np.nansum(sigma_data_**2 / n))
    sigma_err_mean = sigma_err_var / np.sqrt(n / oversample)
    sigma_err_median = sigma_err_mean * 1.25

    return sigma_err_mean, sigma_err_median

#not used yet...
def sigma_err_2d(data1_err, data2_err):

    '''Calculate uncertainy for bin considering both x and y data'''

    n = len(data1_err)

    sigma1_err_ = np.nansum(data1_err**2)**0.5 / n
    sigma2_err_ = np.nansum(data2_err**2)**0.5 / n
    sigma_err = (sigma1_err_**2 + sigma2_err_**2)**0.5

    return sigma_err
####

### Get binned y data, with defined bins in x data
def get_bins_1d(data1, data2, data1_err, data2_err, oversample=1, bins='', nbins=10, logbins=False):

    '''Takes two numpy arrays of data and computes the binned data.
        INPUT:
            data1 = x-axis data - accepts multi-dimentional arrays
            data2 = y-axis data - accepts multi-dimentional arrays
            data1_err = x-axis data err (not used)
            data2_err = y-axis data err
            oversample = data oversampling ratio, accounting that data points oversample beam, and are therefore correlated and not independent measuremnets- usually N_oversample = 1.13 * beam_area / pixel_area; Default = 1
            bins = bin ranges, i.e. edges of bins; Default = determine bins from min and max x-data values
            nbins = if bins not defined, then number of bins to deteremine; Default = 10
            logbins = if bins not defined, log spaced bins for data >0
        OUTPUT:
            x = median binned data x
            y = median binned data y
            stats = statistics of binned data:                
                ['significance'] = signficance of each bin (i.e. S/N) - only y-data error taken into account
                ['ybin_err'] = error in each bin - only y-data error taken into account
                ['-1sigma'] = -1 sigma (15.9 pecentile) of y-data within bin
                ['+1sigma'] = +1 sigma (84.1 pecentile) of y-data within bins
                ['nneg'] = number of negative data points within with bins (useful for plotting)
                ['npos'] = number of posative data points within with bins (useful for plotting)
                ['ntot'] = number of all data points within with bins (useful for plotting)
                ['isnan'] = number of nan data points within with bins (useful for error checks)
                ['isnotnan'] = number of non-nan data points within with bins (useful for error checks)
            bins = bin ranges, i.e. edges of bins
        '''

    data1 = data1.flatten()
    data2 = data2.flatten()
    data1_err = data1_err.flatten()
    data2_err = data2_err.flatten()
    
    #Remove nan values in data2 -> may screw with stats!
    id_nonnan = np.where(~np.isnan(data2))
    data1 = data1[id_nonnan]
    data2 = data2[id_nonnan]
    data1_err = data1_err[id_nonnan]
    data2_err = data2_err[id_nonnan]
    
    id_nonnan = np.where(~np.isnan(data2_err))
    data1 = data1[id_nonnan]
    data2 = data2[id_nonnan]
    data1_err = data1_err[id_nonnan]
    data2_err = data2_err[id_nonnan]

    if bins=='':
        min=np.nanmin(data1)
        max=np.nanmax(data1)
        if logbins:
            min = np.nanmin(data1[data1>0])
            bins = np.logspace(np.log10(min), np.log10(max), nbins+1)
        else:
            bins = np.linspace(min, max, nbins+1)
    else:
        nbins = len(bins)-1

    x = np.empty([nbins]) *np.nan
    y = np.empty([nbins]) *np.nan
    p1 = np.empty([nbins]) *np.nan
    p2 = np.empty([nbins]) *np.nan

    sigma = np.empty([nbins]) *np.nan
    significant = np.empty([nbins]) *np.nan
    neg = np.empty([nbins]) *np.nan
    pos = np.empty([nbins]) *np.nan
    ntot = np.empty([nbins]) *np.nan
    isnan = np.empty([nbins]) *np.nan
    isnotnan = np.empty([nbins]) *np.nan
    
    for i in range(nbins):

        ids_ = (data1 > bins[i]) & (data1 < bins[i+1])
        ids = np.where(ids_)

        y[i] = np.nanmedian(data2[ids])
        x[i] = np.nanmedian([data1[ids], data1[ids]])

        p1[i] = np.nanpercentile(data2[ids], 50 - 34.1)
        p2[i] = np.nanpercentile(data2[ids], 50 + 34.1)

        sigma_err_mean, sigma_err_median = sigma_err(data2_err[ids], oversample)

        sigma[i] = sigma_err_median
        significant[i] = y[i] / sigma[i]

        neg[i] = len(np.where(data2[ids]<=0)[0])
        pos[i] = len(np.where(data2[ids]>0)[0])
        ntot[i] = len(data2[ids])
        isnan[i] = len(np.where(np.isnan(data2[ids]))[0])
        isnotnan[i] = len(np.where(~np.isnan(data2[ids]))[0])
        
        if isnan[i] > 0: 
            print('[Warning] Bin contains nan values that should have been taken care of!')
            
    id_ = ~np.isnan(significant)
    significant = significant[id_]
    sigma = sigma[id_]
    p1 = p1[id_]
    p2 = p2[id_]
    x = x[id_]
    y = y[id_]
    ntot = ntot[id_]
    isnan = isnan[id_]
    isnotnan = isnotnan[id_]
    neg = neg[id_]
    pos = pos[id_]
    
    keys = ['significance', 'ybin_err', '-1sigma', '+1sigma', 'nneg', 'npos', 'ntot', 'isnotnan', 'isnan']

    stats = dict.fromkeys(keys)

    stats['significance'] = significant
    stats['ybin_err'] = sigma
    stats['-1sigma'] = p1
    stats['+1sigma'] = p2
    stats['nneg'] = neg
    stats['npos'] = pos
    stats['ntot'] = ntot
    stats['isnan'] = isnan
    stats['isnotnan'] = isnotnan
    
#     stats = np.array([significant, p1, p2, neg, pos])

    return x, y, stats, bins
