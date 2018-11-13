import numpy as np
import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube

def get_momentmaps(cube, outputfile, velocity_range = '', 
            moments = ['mom0', 'max', 'mom1', 'sigma', 'mom9', 'mom9-mom1'], 
            singlefile = True,
            ):

    if velocity_range != '':
        cube_slab = cube.spectral_slab(velocity_range[0] * u.km / u.s, velocity_range[1] * u.km / u.s)
    else: 
        cube_slab = cube
    
    moment_0 = cube_slab.moment(order=0, axis=0)
    maximum = cube_slab.max(axis = 0)

    if outputfile != '':

        print 'Creating intensity moments'

        hdu = fits.HDUList()
        primary = fits.PrimaryHDU()

        if 'mom0' in moments:
            image1 = fits.ImageHDU(data = np.asarray(moment_0), header = moment_0.header, name = 'mom0')
            hdu.append(image1)
            
        if 'max' in moments:
            image2 = fits.ImageHDU(data = np.asarray(maximum), header = maximum.header, name = 'max')
            hdu.append(image2)

        if ['mom1', 'sigma', 'mom9', 'mom9-mom1'] in moments:
            
            print 'Creating velocity moments'
            
            moment_1 = cube_slab.moment(order=1, axis=0)
            sigma_map = cube_slab.linewidth_sigma()

            spectral_axis = cube_slab.spectral_axis
            moment_9 = np.empty([cube_slab.shape[1], cube_slab.shape[2]]) * np.nan
            moment_9_ID = cube_slab.argmax(axis = 0)

            for velo in ProgressBar(range(1, len(spectral_axis))):
                for x in range(cube_slab.shape[1]):
                    for y in range(cube_slab.shape[2]):

                        if moment_9_ID[x,y] == 0:
                            moment_9[x,y] = np.nan
                        else:
                            moment_9[x,y] = spectral_axis[moment_9_ID[x,y]].value

            velo_complexity = np.absolute( (np.asarray(moment_1) - moment_9) )
            velo_complexity = velo_complexity / np.nanmean(np.asarray( sigma_map ))

            if 'mom1' in moments:
                image3 = fits.ImageHDU(data = np.asarray(moment_1), header = moment_1.header, name = 'mom1')
                hdu.append(image3)
            if 'sigma' in moments:
                image4 = fits.ImageHDU(data = np.asarray(sigma_map), header = sigma_map.header, name = 'sigma')
                hdu.append(image4)
            if 'mom9'in moments:
                image5 = fits.ImageHDU(data = np.asarray(moment_9), header = maximum.header, name = 'mom9')
                hdu.append(image5)
            if 'mom9-mom1' in moments:
                image6 = fits.ImageHDU(data = np.asarray(velo_complexity), header = moment_1.header, name = 'mom9-mom1')
                hdu.append(image6)        
            
        if ['variance', 'fwhm'] in moments:
            
            print 'Creating extra moments.'
            
            variance_map = cube_slab.moment(order=2, axis=0)
            fwhm_map = cube_slab.linewidth_fwhm()
            
            if 'variance' in moments:
                image7 = fits.ImageHDU(data = np.asarray(variance_map), header = variance_map.header, name = 'variance')
                hdu.append(image7)
            
            if 'fwhm' in moments:
                image8 = fits.ImageHDU(data = np.asarray(fwhm_map), header = fwhm_map.header, name = 'fwhm')
                hdu.append(image8)

        if singlefile:
            hdu.writeto(outputfile, overwrite = True)
            
        else:            
            if 'mom0' in moments:
                moment_0.write(outputfile+'_mom0.fits', overwrite = True)
                
            if 'max' in moments:
                maximum.write(outputfile+'_mom9.fits', overwrite = True)
                
            if 'mom1' in moments:
                moment_1.write(outputfile+'_mom1.fits', overwrite = True)
                
            if 'sigma' in moments:
                sigma_map.write(outputfile+'_sigma.fits', overwrite = True)
                
            if 'mom9-mom1' in moments:
                velo_complexity.write(outputfile+'_mom9-mom1.fits', overwrite = True)
                
            if 'variance' in moments:
                variance_map.write(outputfile+'_variance.fits', overwrite = True)
                
            if 'fwhm' in moments:
                fwhm_map.write(outputfile+'_fwhm.fits', overwrite = True)
                
    return hdu