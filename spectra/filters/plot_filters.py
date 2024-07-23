import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

SPLUS = ['SPLUS_U', 'SPLUS_F378', 'SPLUS_F395', 'SPLUS_F410',
         'SPLUS_F430', 'SPLUS_G', 'SPLUS_F515', 'SPLUS_R',
         'SPLUS_F660', 'SPLUS_I', 'SPLUS_F861', 'SPLUS_Z']

SDSS = ['SDSS_U', 'SDSS_G', 'SDSS_R', 'SDSS_I', 'SDSS_Z']

PS = ['PS_G', 'PS_R', 'PS_I', 'PS_Z']

DES = ['DES_G', 'DES_R', 'DES_I', 'DES_Z']

SM = ['SM_U', 'SM_V', 'SM_G', 'SM_R', 'SM_I', 'SM_Z']

GAIA_DR2 = ['GAIA_BP_DR2', 'GAIA_G_DR2', 'GAIA_RP_DR2']

GAIA_DR2Revised = ['GAIA_G', 'GAIA_BP', 'GAIA_RP']

GALEX = ['GALEX_FUV', 'GALEX_NUV']

# Estimate central wavelength <- very approximate (not using interpolation)

def central_wl(data):
    l = data.iloc[:,0].values
    S = data.iloc[:,1].values
    
    S_peak = S.max()
    l_peak = l[S == S_peak][0]
    
    dS = S - S_peak/2
    
    # Find l1
    
    l1 = l[l < l_peak][abs(dS[l < l_peak]) == abs(dS[l < l_peak]).min()]
    l2 = l[l > l_peak][abs(dS[l > l_peak]) == abs(dS[l > l_peak]).min()]
    
    lcenter = (l2 + l1)/2
    
    return lcenter

def wavelength_to_hex(wavelength, gamma=0.8):
    '''This converts a given wavelength of light to an
    approximate hex color value. The wavelength must be given
    in angstrons in the range from 3000 A through 10000 A
    
    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    
    wavelength = float(wavelength) / 10.
    
    if wavelength >= 200 and wavelength < 300:
        attenuation = 0.3 + 0.7 * (wavelength - 200) / (440 - 200)
        R = ((-(wavelength - 440) / (440 - 200)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 300 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 300) / (440 - 300)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 1000:
        attenuation = 0.3 + 0.7 * (1000 - wavelength) / (1000 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    R *= 255
    G *= 255
    B *= 255
    
    color = '#%02x%02x%02x' % (R, G, B)
    return color

################################################################################

def plot_filters(filter_list, ax, zorder = 1, alpha = 1, legend = False, loc = 6):
    
    for filter in filter_list:
        
        
        # load data
        filter_data = pd.read_csv(filter+".dat", delim_whitespace=True)
        
        Norm = np.nanmax(filter_data.iloc[:,1])
        
        # calculate central wl
        wl_central = central_wl(filter_data)
        
        # get color for the plot
        color = wavelength_to_hex(wl_central[0])
        print "filter: {} | color: {}".format(filter, color)
        
        linewidth = 1
        
        # plot the filter
        x = filter_data.iloc[:,0]
        y = filter_data.iloc[:,1]/Norm
        x = np.concatenate(([x[0]], x, [list(x)[-1]]))
        y = np.concatenate(([0], y, [0]))
        
        if zorder != -1:
            ax.fill_between(x, 0, y,
                color = color, zorder = zorder, alpha = 0.3)
            
            linewidth = 1
        
        ax.plot(x, y,
                color = color, zorder = zorder, alpha = alpha,
                label = filter, linewidth = linewidth)
        
    
    if legend:
        ncol = 2 if len(filter_list) > 4 else 1
        ax.legend(loc = loc, ncol = ncol, fontsize = 9)
        
    
    ax.set_yticks([])
    
    ax.set_xlim((800, 11100))
    ax.set_ylim((0, 1.1))



fig = plt.figure(figsize=(12, 12)) 

gs = gridspec.GridSpec(6, 2, width_ratios=[0.99, 0.01]) 



ax1 = plt.subplot(gs[0,0])
plot_filters(SPLUS, ax = ax1, legend = True)
ax1.text(11000, 0.83, 'S-PLUS', fontsize = 14, fontweight = 'bold', color = "#666666", ha = 'right')


ax2 = plt.subplot(gs[1,0])
plot_filters(SM, ax = ax2, legend = True)
plot_filters(SPLUS, ax = ax2, zorder = -1, alpha = 0.4)
ax2.text(11000, 0.83, 'Skymapper', fontsize = 14, fontweight = 'bold', color = "#666666", ha = 'right')

ax3 = plt.subplot(gs[2,0])
plot_filters(SDSS, ax = ax3, legend = True)
plot_filters(SPLUS, ax = ax3, zorder = -1, alpha = 0.4)
ax3.text(11000, 0.83, 'SDSS', fontsize = 14, fontweight = 'bold', color = "#666666", ha = 'right')


ax4 = plt.subplot(gs[3,0])
plot_filters(PS, ax = ax4, legend = True)
plot_filters(SPLUS, ax = ax4, zorder = -1, alpha = 0.4)
ax4.text(11000, 0.83, 'PanSTARRS/\nATLAS RefCat', fontsize = 14, fontweight = 'bold', color = "#666666", ha = 'right')


#ax5 = plt.subplot(gs[4,0])
#plot_filters(DES, ax = ax5, legend = True)
#plot_filters(SPLUS, ax = ax5, zorder = -1, alpha = 0.2)
#ax5.text(11000, 0.83, 'DES', fontsize = 14, fontweight = 'bold', color = "#666666", ha = 'right')


ax6 = plt.subplot(gs[4,0])
plot_filters(GAIA_DR2Revised, ax = ax6, legend = True)
plot_filters(SPLUS, ax = ax6, zorder = -1, alpha = 0.4)
ax6.text(11000, 0.83, 'Gaia DR2 Revised', fontsize = 14, fontweight = 'bold', color = "#666666", ha = 'right')


ax7 = plt.subplot(gs[5,0])
plot_filters(GALEX, ax = ax7, loc = 7, legend = True)
plot_filters(SPLUS, ax = ax7, zorder = -1, alpha = 0.4)
ax7.text(11000, 0.83, 'GALEX', fontsize = 14, fontweight = 'bold', color = "#666666", ha = 'right')



#ax4.set_ylabel("Normalized transmission")


ax7.set_xlabel("Wavelength (A)")

ax1.axes.get_xaxis().set_visible(False)
ax2.axes.get_xaxis().set_visible(False)
ax3.axes.get_xaxis().set_visible(False)
ax4.axes.get_xaxis().set_visible(False)
#ax5.axes.get_xaxis().set_visible(False)
ax6.axes.get_xaxis().set_visible(False)

plt.subplots_adjust(top=0.98,bottom=0.05,left=0.02,right=0.99,hspace=0.0,wspace=0.0)
plt.savefig("filters_v2.png")
plt.clf()
plt.close()

