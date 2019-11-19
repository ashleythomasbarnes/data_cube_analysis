import numpy as np

#### Uncertainty
def sigma_err(sigma_data, oversample):

    '''Calculate uncertainy for bin considering only y data'''

    n = len(sigma_data)

    sigma_err_var = np.sqrt(np.nansum(sigma_data**2 / n))
    sigma_err_mean = sigma_err_var / np.sqrt(n / oversample)
    sigma_err_median = sigma_err_mean * 1.25

    return sigma_err_mean, sigma_err_median

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
            data1_err = x-axis data err
            data2_err = x-axis data err
            oversample = data oversampling ratio, accounting that data points oversample beam, and are therefore correlated and not independent measuremnets- usually N_oversample = 1.13 * beam_area / pixel_area; Default = 1
            bins = bin ranges, i.e. edges of bins; Default = determine bins from min and max x-data values
            nbins = if bins not defined, then number of bins to deteremine; Default = 10
            logbins = if bins not defined, log spaced bins for data >0
        OUTPUT:
            x = median binned data x
            y = median binned data y
            stats = statistics of binned data:
                [0]: signficance of each bin (i.e. S/N) - only y-data error taken into account
                [1]: -1 sigma (15.9 pecentile) of y-data within bin
                [2]: +1 sigma (84.1 pecentile) of y-data within bins
                [3]: number of negative data points within with bins (useful for plotting)
                [4]: number of posative data points within with bins (useful for plotting)
            bins = bin ranges, i.e. edges of bins
        '''

    data1 = data1.flatten()
    data2 = data2.flatten()
    data1_err = data1_err.flatten()
    data2_err = data2_err.flatten()

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

    for i in range(nbins):

        ids_ = (data1 > bins[i]) & (data1 < bins[i+1])
        ids = np.where(ids_)

        y[i] = np.nanmedian(data2[ids])
        x[i] = np.nanmedian([data1[ids], data1[ids]])

        p1[i] = np.nanpercentile(data2[ids], 50 - 34.1)
        p2[i] = np.nanpercentile(data2[ids], 50 + 34.1)

        sigma_err_mean, sigma_err_median = sigma_err(data1_err[ids], oversample)

        sigma[i] = sigma_err_median
        significant[i] = y[i] / sigma[i]

        neg[i] = len(np.where(data2[ids]<=0)[0])
        pos[i] = len(np.where(data2[ids]>0)[0])

    id_ = ~np.isnan(significant)
    significant = significant[id_]
    p1 = p1[id_]
    p2 = p2[id_]
    x = x[id_]
    y = y[id_]

    stats = np.array([significant, p1, p2, neg, pos])

    return x, y, stats, bins
