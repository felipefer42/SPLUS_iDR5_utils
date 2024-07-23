"""
Created on 03/11/2017

@Author: Carlos Eduardo Barbosa

Calculates the area of observations from a selection of titles.

- Felipe Almeida Fernandes: Bug fixed in 2020
    Updated make_polygon adding a 1/cos(DEC) factor
    to take into account projection effects
    
    Updated the script to use sys.argv input
"""

from __future__ import division, print_function

import sys
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from matplotlib.collections import PatchCollection
from matplotlib import patches
from spherical_geometry.polygon import SphericalPolygon

class SurveyArea():
    def __init__(self, coords):
        self.coords = coords
        self.make_polygons()
        self.calc_area()
        return

    def make_polygons(self):
        """ Produces squares with the vertices of the polygons used to
        construct a footprint. """
        delta = np.sqrt(2) / 2  # * u.degree
        self.polygons = []
        for i, coord in enumerate(self.coords):
            ra0 = coord.ra.to(u.degree).value
            dec0 = coord.dec.to(u.degree).value

            dec_low = np.radians(dec0 - delta)
            dec_high = np.radians(dec0 + delta)

            ras = np.array([ra0 - delta / np.cos(dec_low),
                            ra0 - delta / np.cos(dec_high),
                            ra0 + delta / np.cos(dec_high),
                            ra0 + delta / np.cos(dec_low)])

            decs = np.array([dec0 - delta,
                             dec0 + delta,
                             dec0 + delta,
                             dec0 - delta])

            vertices = [(r, d) for r, d in zip(ras, decs)]
            self.polygons.append(vertices)
        return
        
    def calc_area(self):
        """ Determines the area of a survey footprint. """
        polys = [Polygon(_) for _ in self.polygons]
        self.footprint = cascaded_union(polys)
        area = []
        if isinstance(self.footprint, Polygon):
            self.footprint = [self.footprint]
        for poly in self.footprint:
            ra, dec = np.array(poly.exterior.xy)
            sptot = SphericalPolygon.from_radec(ra, dec)
            area.append(sptot.area())
        self.areas = np.array(area) * np.power(180 / np.pi, 2)
        self.total_area = self.areas.sum() * u.degree**2
        return

    def make_collection(self, c="k"):
        """ Determines the projection. """
        xys = []
        for poly in self.polygons:
            coords = np.array(poly)
            x, y = self.radec2projection(coords[:,0], coords[:,1])
            xys.append(np.column_stack((x,y)))
        xys = [Polygon(_) for _ in xys]
        union = cascaded_union(xys)
        if isinstance(union, Polygon):
            union = [union]
        list_of_patches = []
        for i, poly in enumerate(union):
            ext = np.array(poly.exterior.xy)
            list_of_patches.append(patches.Polygon(ext.T))
        p = PatchCollection(list_of_patches, facecolors=c, alpha=0.8)
        return p

    def radec2projection(ra, dec):
        """ Converts coordinates to Mollweide/Aitoff coordinate system in
        matplotlib. """
        xx = np.deg2rad(ra - 180.)
        yy = np.deg2rad(dec)
        return xx, yy

if __name__ == "__main__":
    ############################################################################
    # Point to tiles file here
    filename = sys.argv[1]
    ############################################################################
    table = Table.read(filename, format="csv")
    coords = SkyCoord(ra=table["RA"], dec=table["DEC"],
                           unit=(u.hourangle, u.degree))
    # Example: print area of each sub-survey with different PID
    subsurveys = np.unique(table["PID"])
    for survey in subsurveys:
        idx = np.where(table["PID"] == survey)
        footprint = SurveyArea(coords[idx])
        print(survey, footprint.total_area)
    footprint = SurveyArea(coords)
    print("Total area: ", footprint.total_area)
