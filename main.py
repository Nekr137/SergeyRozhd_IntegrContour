﻿from pandas import read_table
from numpy import linspace, interp, meshgrid, asfarray, abs, mean
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp2d
from math import sqrt

def find_nearest(array, value):
    array = asfarray(array)
    return (abs(array - value)).argmin()

def load_columns(fname, delimiter='\t', decimal=','):
    data = read_table(fname, sep=delimiter,header=None, dtype=float, decimal=decimal)
    return [list(col) for col in data.values.T]

def interp(ax, cnt, xdata, ydata):
    xx = linspace(xdata[0], xdata[-1], cnt)
    yy = interp(xx, xdata, ydata)
    return list(xx),list(yy)

class Perenos:
    def __init__(self, xPos) -> None:
        self._xPos = xPos
        self._depths = []
        self._values = []

    def append(self, depth, value):
        self._depths.append(depth)
        self._values.append(value)
    
    def min_depth(self):
        return min(self._depths)
    
    def max_depth(self):
        return max(self._depths)

def load_perenos(fname):
    data = load_columns(fname, decimal=',')
    perenos = []
    i = 0
    while i < len(data[0]):
        x = data[0][i]
        perenos.append(Perenos(x))
        while i < len(data[0]) and (data[0][i] - x) < 1e-6:
            depth = data[1][i]
            value = data[2][i]
            perenos[-1].append(depth, value)
            i += 1
    return perenos

def show_perenos_pnts(ax, perenos):
    for p in perenos:
        for d in p._depths:
            ax.plot([p._xPos], [d], 'k.', markersize=1)


def interpoate_2d(X, Y, Z):
    Nx, Ny = 90,60
    xx,yy = linspace(X[0][0], X[0][-1], Nx),linspace(Y[0][1], Y[-1][0], Ny)
    xnew, ynew = meshgrid(xx,yy)
    f = interp2d(X,Y,Z, kind='cubic')
    znew = f(xx,yy)
    return xnew, ynew, znew

def build_perenos_data(aPerenos):
    xmin = aPerenos[ 0]._xPos
    xmax = aPerenos[-1]._xPos
    dmin = min([min(p._depths) for p in aPerenos])
    dmax = max([max(p._depths) for p in aPerenos])
    Nx = len(aPerenos)
    Ny = len(aPerenos[0]._depths)
    x = linspace(xmin, xmax, Nx)
    y = linspace(dmin, dmax, Ny)
    X,Y = meshgrid(x,y)
    m = mean([mean(p._values) for p in aPerenos]) # mean values
    D = [[m for i in range(len(x))] for j in range(len(y))]
    for p in aPerenos:
        i = find_nearest(x,p._xPos)
        for idx in range(len(p._depths)):
            value = p._values[idx]
            depth = p._depths[idx]
            j = find_nearest(y, depth)
            D[j][i] = value
    X,Y,D = interpoate_2d(X,Y,D)
    return X,Y,D

class P2D:
    def __init__(self, x=0.0, y=0.0) -> None:
        self.x, self.y = x, y
    
    def dist2(self, other):
        return (other.x - self.x)**2 + (other.y - self.y)**2
    
    def dist(self, other):
        return sqrt(self.dist2(other))

def find_depth_border_polyline(aPerenos):
    pnts = []
    for p in aPerenos:
        x = p._xPos
        y = p.min_depth()
        pnts.append(P2D(x, y))
    return pnts

def find_border_polyline(aPerenos):
    pnts = []
    ul = P2D(aPerenos[0]._xPos, aPerenos[0].max_depth())
    ur = P2D(aPerenos[-1]._xPos, aPerenos[-1].max_depth())
    pnts.append(ul)
    for p in find_depth_border_polyline(aPerenos):
        pnts.append(p)
    pnts.append(ur)
    pnts.append(ul)
    return pnts

def treat_day(ax, fname):
    """
    return: contour length
    """
    perenos = load_perenos(fname)
    show_perenos_pnts(ax, perenos)
    X,Y,D = build_perenos_data(perenos)
    cs = ax.contour(X, Y, D, colors='k', levels=15, linewidths=1, linestyles='solid')
    ax.clabel(cs, inline=True, fontsize=12)
    ext = (X[0][0], X[0][-1], Y[0][0], Y[-1][0])
    im = ax.pcolor(X, Y, D, cmap=cm.jet)

    # find and show borders
    pnts = find_border_polyline(perenos)
    perim = sum([pnts[i].dist(pnts[i+1]) for i in range(len(pnts) - 1)])
    # ax.plot([p.x for p in pnts], [p.y for p in pnts], 'r.-')

    # draw blue polygons
    pnts = find_depth_border_polyline(perenos)
    def make_lower_poly(p1, p2):
        li = ax.get_ylim()
        lo = li[0] #lowest point on the ax
        x = [p1.x, p1.x, p2.x, p2.x]
        y = [lo, p1.y, p2.y, lo]
        ax.fill(x, y, color='k')
    for i in range(len(pnts) - 1):
        make_lower_poly(pnts[i], pnts[i+1])

    return perim

def main():
    fig, axs = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(20,6)

    axs[0].set_xlabel('Distance (km)')
    axs[1].set_xlabel('Distance (km)')
    axs[0].set_title('J(m/sec)*(mg/l)')
    axs[1].set_title('J(m/sec)*(mg/l)')
    axs[0].set_ylabel('Depth (m)')

    perim10 = treat_day(axs[0], '10.08.copy.perenos.xlsx.dat')
    perim11 = treat_day(axs[1], '11.08.copy.perenos.xlsx.dat')

    # fig0.tight_layout()
    fig.savefig('output', dpi=300)
    pass

if __name__ == '__main__':
    main()
