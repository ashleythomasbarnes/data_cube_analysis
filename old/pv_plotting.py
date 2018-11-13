def pv_plot(comps, output, xlim = [], ylim = [], size_factor = 5, arm1 = [], arm2 = []):

    fig = plt.figure(figsize = (5.5,4.5))
    ax = fig.add_subplot(111)

    size_max = np.nanmax(comps[:,:,3])
    for m in range(number_comps):
        ax.scatter(comps[m,:,2], comps[m,:,5], s = (comps[m,:,3] * size_factor) / size_max, c = 'grey', edgecolor='none', zorder = 10)

    if xlim != [] or ylim != []:
        ax.set_xlim(ylim); ax.set_ylim(zlim)

    ax.set_xlabel(r'$\Delta$ Dec (J2000, arcsec)')
    ax.set_ylabel(r'Velocity (km\,s$^{-1}$)')

    if arm1 != [] or arm2 != []:
        ax.fill_between(arm1[0], arm1[2] - 5., arm1[2] + 5., color = 'blue', alpha = 0.1)
        ax.fill_between(arm2[0], arm2[2] - 5., arm2[2] + 5., color = 'green', alpha = 0.1)

        ax.plot(arm1[0], arm1[2], zorder = 1, color = 'blue')
        ax.plot(arm2[0], arm2[2], zorder = 1, color = 'green')

    # ax.text(250, 31, 'Sgr near', color = 'green', fontweight = 'bold', ha = 'left')
    # ax.text(250, 52, 'Sgr far', color = 'blue', fontweight = 'bold', ha = 'left')
    # ax.text(0, 65, 'GRS velocity components \n (Barnes et al.) ', color = 'grey', ha = 'left')

    plt.tight_layout()
    rasterize = True
    fig.set_rasterized(rasterize)
    plt.savefig(output,  dpi = 300)
    print 'Saving: '+output
    plt.close('all')



def pvv_plot(comps, arr_file, output, aspect = 1.0, elev=0, azim=0, xlim = [], ylim = [], zlim = [], size_factor = 100, x_inv = True):

    fig = plt.figure(figsize = (5.5,4.5))
    ax = fig.add_subplot(111)

    ax.pbaspect = [aspect, 1.0, 1.0]
    if xlim != [] and ylim != [] and zlim != []:
        ax.set_xlim(ylim); ax.set_ylim(zlim)
        ax.set_xlim(xlim[0], xlim[1]); ax.set_ylim(ylim[0], ylim[1]); ax.set_zlim(zlim[0], zlim[1])

    if x_inv:
        ax.invert_xaxis()
    ax.set_xlabel(r'$\Delta$ RA (J2000, arcsec) '); ax.set_ylabel(r'$\Delta$ Dec (J2000, arcsec)'); ax.set_zlabel(r'Velocity (km\,s$^{-1}$)')

    arr = read_png(arr_file)
    stepX, stepY = (np.absolute(xlim[0])+np.absolute(xlim[1])) / arr.shape[0], (np.absolute(ylim[0]) + np.absolute(ylim[1])) / arr.shape[1]
    X1 = np.arange(xlim[0], xlim[1], stepX)
    Y1 = np.arange(ylim[0], ylim[1], stepY)
    X1, Y1 = np.meshgrid(X1, Y1)
    ax.plot_surface(np.transpose(X1), np.transpose(Y1), zlim[0], rstride=5, cstride=5, facecolors = arr , zorder = 0, antialiased=True, norm = 0.5, cmap = 'gist_yarg')

    #### IRAM
    number_comps = len(comps)
    colours = iter(cm.gist_rainbow(np.linspace(0,1,number_comps)))
    size_max = np.nanmax(comps[:,:,3])

    for m in range(number_comps):
        c = next(colour)
        ax.scatter(comps[m,:,1] + offset[0], comps[m,:,2] + offset[1], comps[m,:,5], s = (comps[m,:,3] * size_factor) / size_max, color=c, edgecolor='none', zorder = 10, depthshade = 1)
        ax.collections[-1].__class__ = FixZorderScatter

    ax.view_init(elev=elev, azim=azim)
    plt.tight_layout()
    rasterize = True
    fig.set_rasterized(rasterize)
    plt.savefig(output,  dpi = 300)
    print 'Saving: '+output
    plt.close('all')
