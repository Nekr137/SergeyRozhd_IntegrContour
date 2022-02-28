
class Cntr:
    def to_float(s):
        assert(len(s) == 4 or len(s) == 5)
        if len(s) == 5:
            return -1.0 * float(s[1:])
        if len(s) == 4:
            return float(s)
        return 0.0

    def get_bbox(t, fig, ax):
        transf = ax.transData.inverted()
        bb = t.get_window_extent(renderer = fig.canvas.renderer)
        return bb.transformed(transf)

    def intersected(b1, b2):
        return not (b2.x0>b1.x1 or b2.x1<b1.x0 or b2.y1<b1.y0 or b2.y0>b1.y1)
    
    def adjust_labels(texts, ax):
        fig = ax.figure
        n = len(texts)
        texts = sorted(texts, key=lambda it: -Cntr.to_float(it.get_text()))
        for i in range(n-1):
            ti = texts[i]
            if not ti.get_visible():
                continue
            for j in range(i+1,n):
                tj = texts[j]
                if not tj.get_visible():
                    continue
                bi = Cntr.get_bbox(ti, fig, ax)
                bj = Cntr.get_bbox(tj, fig, ax)
                if Cntr.intersected(bi, bj):
                    vi = Cntr.to_float(ti.get_text())
                    vj = Cntr.to_float(tj.get_text())
                    if vi > vj:
                        tj.set_visible(False)
                    else:
                        ti.set_visible(False)

def sr_contour(ax, X, Y, D):
    levels=20
    linewidth=1
    linestyle='solid'
    use_clabeltext=True
    fontsize=9
    inline=True
    colors='k'

    # building fake contour to find labels
    cs = ax.contour(X, Y, D, colors=colors, levels=levels, linewidths=linewidth, linestyles=linestyle, use_clabeltext=use_clabeltext)
    cl1 = ax.clabel(cs, inline=inline, fontsize=fontsize)

    # removing overlapped labels
    Cntr.adjust_labels(cl1, ax)
    manual_pos = [c.get_position() for c in cl1 if c.get_visible() == True]

    # build real contour
    cs2 = ax.contour(X, Y, D, colors=colors, levels=levels, linewidths=linewidth, linestyles=linestyle, use_clabeltext=use_clabeltext)
    cl2 = ax.clabel(cs2, inline=inline, fontsize=fontsize, manual=manual_pos)

