import numpy as np

# TO DO: calculated value weighted histogram

def get_hist(data, bins, norm=True, cum=False):
    
    '''Get (area weighted) histograms from data, and return bin centre and histogram data.
    Input:
        data = input data 
        bins = input bin edges for histogram calculaiton
        norm = normalise such that max is equal to unity 
        cum = cumulative distobution; otherwise probability distobution
    Return: 
        bin_cent = bin centres, for easy plotting in matplotlib
        hist = histogram data for each bin centre
        '''
    
    data = data.flatten()
    bins_cent = np.empty([len(bins)-1])
    for i in range(len(bins)-1):
        bins_cent[i] = np.nanmean([bins[i], bins[i+1]])
    
    hist = np.histogram(data.flatten(), bins)[0]
    
    if cum: 
        hist = np.cumsum(hist)
    if norm: 
        hist = hist/np.nanmax(hist)
    
    return(bins_cent, hist)