#! /usr/bin/env python
import numpy as np

#------------------------------------------------------------------------
#---- function
#------------------------------------------------------------------------

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    lon1, lat1 := scalars
    lon2, lat2 := 1D arrays
    """

    Re = 6378.137

    # convert decimal degrees to radians
    deg2rad = np.pi / 180.
    lon1 = np.array(lon1) * deg2rad
    lat1 = np.array(lat1) * deg2rad
    lon2 = np.array(lon2) * deg2rad
    lat2 = np.array(lat2) * deg2rad

    if lon2.shape:
        N = lon2.shape[0]
        lon1 = np.repeat(lon1,N)
        lat1 = np.repeat(lat1,N)

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2.)**2. + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.)**2.
    c = 2. * np.arcsin(np.sqrt(a))
    km = Re * c
    return km

#------------------------------------------------------------------
#--- FUNCTION
#------------------------------------------------------------------

def inpolygon(polygon, xp, yp):
    from shapely.geometry import Point
    return np.array([Point(x, y).intersects(polygon) for x, y in zip(xp, yp)],
                    dtype=np.bool)

#------------------------------------------------------------------
#--- FUNCTION
#------------------------------------------------------------------
def points_in_range(clon,clat,plon,plat,range_km):
    if hasattr(clon, '__len__'):
        m = np.zeros(plon.shape,dtype=bool)
        for cx,cy in zip(clon,clat):
            mi = points_in_range(cx,cy,plon,plat,range_km)
            m = (m | mi)
    else:
        mask = np.array((haversine(clon,clat,plon.ravel(),plat.ravel()) <= range_km))
        if mask.ndim != plon.ndim:
            m = mask.reshape(plon.shape)
        else:
            m = mask
    return m


#------------------------------------------------------------------
#--- FUNCTION
#------------------------------------------------------------------

def ocean_points(x,y):
    from cartopy.feature import NaturalEarthFeature
    from shapely.ops import cascaded_union
    shp = NaturalEarthFeature(category='physical', name='ocean',
                              scale='50m')
    geoms = shp.geometries()
    polygon = cascaded_union(list(geoms))
    mask = inpolygon(polygon, x.ravel(), y.ravel())
    if mask.ndim != x.ndim:
        m = mask.reshape(x.shape)
    else:
        m = mask
    return m

#------------------------------------------------------------------
#--- test
#------------------------------------------------------------------

if __name__ == '__main__':
    import matplotlib as mpl
    #mpl.use('Agg')
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    x = np.arange(-180.,185.,5)
    y = np.arange(-90,95,5)
    X,Y = np.meshgrid(x,y)

    import xarray as xr
    ds = xr.open_dataset('/glade/p/work/mclong/orcas/gv_merged_data/ORCASall.merge10.20170526.nc')
    X = ds.GGLON.values
    Y = ds.GGLAT.values

    m = ocean_points(X,Y)
    #clon = []
    #clat = []
    #for k,v in named_points.items():
    #    clon.append(v[1])
    #    clat.append(v[0])

    #m = points_in_range(clon,clat,X,Y,100)

    plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()

    #ax.set_extent([-180, 180, -90, 90],crs=ccrs.PlateCarree())

    ax.plot(X[~m],Y[~m],'r.',alpha=1,transform=ccrs.PlateCarree(),markersize=1)
    ax.plot(X[m],Y[m],'b.',alpha=1,transform=ccrs.PlateCarree(),markersize=1)

    plt.savefig('earth_geometry_test_oceanland.png',bbox_inches='tight')
    plt.close()
