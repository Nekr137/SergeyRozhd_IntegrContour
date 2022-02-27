from pandas import read_table
from numpy import linspace, interp, meshgrid, asfarray, abs, mean
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp2d
from math import sqrt
from copy import deepcopy

from Perenos import Perenos

DARK_BLUE = (0.0, 0.0, 0.3)

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

def diff(aValues):
    return [aValues(i+1) - aValues(i) for i in range(len(aValues) - 1)]

def make_perenos_squared(aPerenos):
    total_depth_bottom = min([p.get_min_depth()[0] for p in aPerenos])
    aPerenosCopy = deepcopy(aPerenos)
    for p in aPerenosCopy:
        depth_step = p.get_depth_step()
        sorted_by_depths = p.get_sorted_by_depths()
        curr_bottom = sorted_by_depths[0][0]
        curr_value = sorted_by_depths[0][1]
        new_bottom = curr_bottom - depth_step
        while new_bottom > total_depth_bottom:
            p.append(new_bottom, curr_value)
            new_bottom -= depth_step
    return aPerenosCopy

def show_perenos_pnts(ax, perenos, style='k.'):
    for p in perenos:
        for v in p.get_sorted_by_depths():
            ax.plot([p._xPos], [v[0]], style, markersize=3)

def interpoate_2d(X, Y, Z):
    Nx, Ny = 90, 60
    xx,yy = linspace(X[0][0], X[0][-1], Nx),linspace(Y[0][1], Y[-1][0], Ny)
    xnew, ynew = meshgrid(xx,yy)
    f = interp2d(X,Y,Z, kind='cubic')
    znew = f(xx,yy)
    return xnew, ynew, znew

def build_perenos_data(aPerenos):
    xmin = aPerenos[ 0]._xPos
    xmax = aPerenos[-1]._xPos
    dmin = min([min(p.get_depths()) for p in aPerenos])
    dmax = max([max(p.get_depths()) for p in aPerenos])
    Nx = len(aPerenos)
    Ny = len(aPerenos[0].get_depths())
    x = linspace(xmin, xmax, Nx)
    y = linspace(dmin, dmax, Ny)
    X,Y = meshgrid(x,y)
    D = [[1e6 for i in range(Nx)] for j in range(Ny)]
    for i in range(Nx):
        for j in range(Ny):
            p1 = P2D(x[i], y[j])
            min_dist2 = 1e12
            for p in aPerenos:
                for v in p.get_sorted_by_depths():
                    p2 = P2D(p._xPos, v[0])
                    dist2 = p1.dist2(p2)
                    if dist2 < min_dist2:
                        min_dist2 = dist2
                        D[j][i] = v[1]
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
        y = p.get_min_depth()[0]
        pnts.append(P2D(x, y))
    return pnts

def find_border_polyline(aPerenos):
    pnts = []
    ul = P2D(aPerenos[0]._xPos, aPerenos[0].get_max_depth()[0])
    ur = P2D(aPerenos[-1]._xPos, aPerenos[-1].get_max_depth()[0])
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
    perenos_orig = load_perenos(fname)
    perenos = make_perenos_squared(perenos_orig)
    # show_perenos_pnts(ax, perenos,'g.')
    # show_perenos_pnts(ax, perenos_orig,'r.')
    X,Y,D = build_perenos_data(perenos)
    cs = ax.contour(X, Y, D, colors='k', levels=15, linewidths=1, linestyles='solid', use_clabeltext=True)
    cl = ax.clabel(cs, inline=True, fontsize=12)
    pc = ax.pcolor(X, Y, D, cmap=cm.jet)
    ax.figure.colorbar(pc, ax=ax)

    # find and show borders
    pnts = find_border_polyline(perenos_orig)
    perim = sum([pnts[i].dist(pnts[i+1]) for i in range(len(pnts) - 1)])
    # ax.plot([p.x for p in pnts], [p.y for p in pnts], 'r.-')

    # draw blue polygons
    pnts = find_depth_border_polyline(perenos_orig)
    ax.set_xlim(pnts[0].x, pnts[-1].x) # fix limits a little bit
    def make_lower_poly(p1, p2):
        li = ax.get_ylim()
        lo = li[0] #lowest point on the ax
        x = [p1.x, p1.x, p2.x, p2.x]
        y = [lo, p1.y, p2.y, lo]
        ax.fill(x, y, color=DARK_BLUE, zorder=10)
    for i in range(len(pnts) - 1):
        make_lower_poly(pnts[i], pnts[i+1])

    return perim

def main():
    fig, axs = plt.subplots(nrows=1, ncols=2)
    w = 18.0 
    h = w / 10.0 * 3.0
    fig.set_size_inches(w,h)

    fontsize = 14
    axs[0].set_xlabel('Distance (km)', fontsize=fontsize)
    axs[1].set_xlabel('Distance (km)', fontsize=fontsize)
    axs[0].set_title('J(m/sec)*(mg/l)   10.08.2021', fontsize=fontsize)
    axs[1].set_title('J(m/sec)*(mg/l)   11.08.2021', fontsize=fontsize)
    axs[0].set_ylabel('Depth (m)', fontsize=fontsize)
    axs[1].set_ylabel('Depth (m)', fontsize=fontsize)

    perim10 = treat_day(axs[0], '10.08.copy.perenos.xlsx.dat')
    perim11 = treat_day(axs[1], '11.08.copy.perenos.xlsx.dat')

    # fig0.tight_layout()
    fig.savefig('output', dpi=300)
    pass

if __name__ == '__main__':
    main()
