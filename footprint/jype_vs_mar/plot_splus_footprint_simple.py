################################################################################
# Extinction Correction
import pandas as pd
import numpy as np
import matplotlib.colors as colors

from cmcrameri import cm

from matplotlib import pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.patches import Polygon

#fields_splus = ["Main.csv", "STRIPE82.csv", "Galactic.csv", "MC.csv"]
#fields_color = ["#ec5649ff", "#bd49c5", "#2c2cf6", "#34ca5e"]
#fields_alpha = [0.5, 0.6, 0.7, 0.8]
#fields_label = ["Main Survey", "STRIPE-82", "Galactic Survey", "Magellanic Clouds"]
fields_splus = ["iDR5_JYPE_full.csv", "iDR5_MAR_full.csv"]
fields_color = ["#57de70", "#a657de"]
fields_alpha = [0.8, 0.8]
fields_label = ["JYPE", "MAR"]

all_splus = "../all_splus_current.csv"
all_color = "#000000"
all_alpha = 0.05

plot_cbar = False
galactic  = False

plot_legend = True

save_name = "footprint_iDR5_reduction.png"

def plot_splus_field(ra, dec, galactic = True, ax = plt, **kargs):
    
    if ra > 180:
        ra = ra - 360

    # Calculate corners
    dec0 = dec - 0.7
    dec1 = dec - 0.7
    dec2 = dec + 0.7
    dec3 = dec + 0.7

    ra0 = ra - 0.7/np.cos(dec0*np.pi/180)
    ra1 = ra + 0.7/np.cos(dec1*np.pi/180)
    ra2 = ra + 0.7/np.cos(dec2*np.pi/180)
    ra3 = ra - 0.7/np.cos(dec3*np.pi/180)
    
    ra = np.array([ra0, ra1, ra2, ra3, ra0])
    dec = np.array([dec0, dec1, dec2, dec3, dec0])
    
    c = SkyCoord(ra*u.degree, dec*u.degree)
    l = c.galactic.l.degree
    b = c.galactic.b.degree
    
    l[l > 180] = l[l > 180] - 360
    
    if not galactic:
        x = -ra
        y = dec
    else:
        x = -l
        y = b
    
    
    y[y > 90] = 90
    y[y < -90] = -90
    
    # Check if crosses border
    f1 = x > 1/np.cos(y*np.pi/180)
    f2 = x < -1/np.cos(y*np.pi/180)
    f3 = f1.sum() * f2.sum()
    if f3 == 0:
        cross_border = False
    else:
        cross_border = True
            
    if not cross_border:   
        ax.fill(x*np.pi/180, y*np.pi/180, **kargs)
    
    if cross_border:
        # Plot posivite side
        x1 = x + 0
        x1[x1 < 0] = 180
        ax.fill(x1*np.pi/180, y*np.pi/180, **kargs)
        
        # Plot negative side
        x2 = x + 0
        x2[x2 > 0] = -180
        ax.fill(x2*np.pi/180, y*np.pi/180, **kargs)

        

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

ebv_cmap = truncate_colormap(cm.grayC_r, 0, 0.5)

# Get all_splus
all_splus = pd.read_csv(all_splus)

# Get eb_v for ism extinction map

ra_list = np.linspace(-180, 180, 1000)
dec_list = np.linspace(-90, 90, 1000)

plt.figure(figsize=(12,6))
plt.subplot(111, projection="aitoff")
plt.grid(True, color = "#444444", linewidth = 0.5)

for i in range(len(all_splus)):
    plot_splus_field(all_splus.loc[i,"RA_d"],all_splus.loc[i,"DEC_d"],galactic=galactic,
                     color = all_color, alpha = all_alpha)

for i in range(len(fields_splus)):
    fields = pd.read_csv(fields_splus[i])
    for j in range(len(fields)):
        plot_splus_field(fields.loc[j,"RA_d"],fields.loc[j,"DEC_d"],galactic=galactic,
                         color = fields_color[i], alpha = fields_alpha[i])
    
    plt.fill_between([], [], [], color = fields_color[i], alpha = fields_alpha[i], label = fields_label[i])

if plot_legend:
    if fields_splus:
        plt.legend(loc = 1, fontsize = 16)
    
if plot_cbar:
    cax = plt.axes([0.92, 0.1, 0.02, 0.8])
    plt.colorbar(cmap, cax = cax)
    plt.ylabel(r"$\log_{10}(E_{B-V})$")
    plt.subplots_adjust(left = 0.05, top = 0.95, bottom = 0.05, right = 0.9)
else:
    plt.subplots_adjust(left = 0.05, top = 0.95, bottom = 0.05, right = 0.95)


labels = ['150°', '120°', '90°', '60°', '30°', '0°', '-30°', '-60°', '-90°', '-120°', '-150°']
plt.gca().set_xticklabels(labels)

plt.savefig(save_name, transparent=False)
plt.clf()
plt.close()

