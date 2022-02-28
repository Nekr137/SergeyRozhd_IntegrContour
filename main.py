from pandas import read_table
from numpy import linspace, interp, meshgrid, asfarray, abs, mean
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp2d
from math import sqrt
from copy import deepcopy
from adjustText import adjust_text

from sr_perenos import Perenos
from sr_contour import sr_contour

M2KM = 1e-3

FONTSIZE = 14
DPI = 300
COLORBAR_LIMITS = (-0.4, 0.4)
POLYGON_COLOR = (0.0, 0.0, 0.3)

INTERP_NX, INTERP_NY = 90, 60 
CONTOUR_INTERP = 500

def find_nearest(array, value):
    array = asfarray(array)
    return (abs(array - value)).argmin()

def load_columns(fname, delimiter='\t', decimal=','):
    data = read_table(fname, sep=delimiter,header=None, dtype=float, decimal=decimal)
    return [list(col) for col in data.values.T]

def interp_line(ax, cnt, xdata, ydata):
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
    xx,yy = linspace(X[0][0], X[0][-1], INTERP_NX),linspace(Y[0][1], Y[-1][0], INTERP_NY)
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

def interp_line(cnt, xdata, ydata):
    xx = linspace(xdata[0], xdata[-1], cnt)
    yy = interp(xx, xdata, ydata)
    return list(xx),list(yy)

def interp_polyline(cnt, polyline):
    x = [p.x for p in polyline]
    y = [p.y for p in polyline]
    tmpx, xx = interp_line(cnt, list(range(len(x))), x)
    tmpy, yy = interp_line(cnt, list(range(len(y))), y)
    return [P2D(xx[i], yy[i]) for i in range(len(xx))]

def find_contour_integr(polyline, X, Y, D, ax):
    pnts = interp_polyline(CONTOUR_INTERP, polyline)
    xx = X[0]
    yy = [v[0] for v in Y]
    sums = [] 
    
    for idx, p in enumerate(pnts):
        i = find_nearest(xx, p.x)
        j = find_nearest(yy, p.y)
        pkm = P2D(p.x, M2KM*p.y)
        value = D[j][i]
        if not idx == 0:
            dist = pkm.dist(prev_pnt_km)
            if dist < 1e-10:
                print('c', end='')
                continue
            sums.append(0.5 * (value + prev_val) * dist)
        prev_val = value
        prev_pnt_km = pkm
    return sum(sums)


def treat_day(ax, fname):
    """
    return: contour length
    """
    perenos_orig = load_perenos(fname)
    perenos = make_perenos_squared(perenos_orig)
    # show_perenos_pnts(ax, perenos,'g.')
    # show_perenos_pnts(ax, perenos_orig,'r.')
    X,Y,D = build_perenos_data(perenos)

    sr_contour(ax, X, Y, D)
    
    pc = ax.pcolor(X, Y, D, cmap=cm.jet, vmin=COLORBAR_LIMITS[0], vmax=COLORBAR_LIMITS[1])
    ax.figure.colorbar(pc, ax=ax)

    polyline = find_border_polyline(perenos_orig)
    # ax.plot([p.x for p in polyline], [p.y for p in polyline], 'r.-')
    integr = find_contour_integr(polyline, X, Y, D, ax)

    # draw blue polygons
    pnts = find_depth_border_polyline(perenos_orig)
    ax.set_xlim(pnts[0].x, pnts[-1].x) # fix limits a little bit
    def make_lower_poly(p1, p2):
        li = ax.get_ylim()
        lo = li[0] #lowest point on the ax
        x = [p1.x, p1.x, p2.x, p2.x]
        y = [lo, p1.y, p2.y, lo]
        ax.fill(x, y, color=POLYGON_COLOR, zorder=100)
    for i in range(len(pnts) - 1):
        make_lower_poly(pnts[i], pnts[i+1])

    return integr

def main():
    fig, axs = plt.subplots(nrows=1, ncols=2)
    w = 18.0 
    h = w / 10.0 * 3.0
    fig.set_size_inches(w,h)

    axs[0].set_xlabel('Distance (km)', fontsize=FONTSIZE)
    axs[1].set_xlabel('Distance (km)', fontsize=FONTSIZE)
    axs[0].set_title('J(m/sec)*(mg/l)   10.08.2021', fontsize=FONTSIZE)
    axs[1].set_title('J(m/sec)*(mg/l)   11.08.2021', fontsize=FONTSIZE)
    axs[0].set_ylabel('Depth (m)', fontsize=FONTSIZE)
    axs[1].set_ylabel('Depth (m)', fontsize=FONTSIZE)

    contour_integr10 = treat_day(axs[0], '10.08.copy.perenos.xlsx.dat')
    contour_integr11 = treat_day(axs[1], '11.08.copy.perenos.xlsx.dat')

    s = 'Contour integr (10.08)\tContour integr (11.08)\n'
    s += '{}\t{}\n'.format(contour_integr10, contour_integr11)
    print(s)
    with open('output.txt', 'w', encoding='utf-8') as f:
        f.write(s)

    # fig0.tight_layout()
    fig.savefig('output', dpi=DPI)
    pass

if __name__ == '__main__':
    main()
