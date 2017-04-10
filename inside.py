from numpy import *
from numpy.random import *
import multiprocessing
from scipy.spatial import ConvexHull
from time import *

def inside(x,y,xpts,ypts):
	xpts=list(xpts)
	ypts=list(ypts)
	#xpts.append(xpts[0])
	#ypts.append(ypts[0])
	n = len(xpts)
	ins = False
	p1x,p1y = xpts[0],ypts[0]
	for i in range(n):
		p2x,p2y = xpts[i],ypts[i]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xints:
						ins = not ins
		p1x,p1y = p2x,p2y

	return ins

def inside_worker(r):
	x=r[0]
	y=r[1]
	xpts=r[2]
	ypts=r[3]
	xpts=list(xpts)
	ypts=list(ypts)
	#xpts.append(xpts[0])
	#ypts.append(ypts[0])
	n = len(xpts)
	ins = False
	p1x,p1y = xpts[0],ypts[0]
	for i in range(n):
		p2x,p2y = xpts[i],ypts[i]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xints:
						ins = not ins
		p1x,p1y = p2x,p2y

	return ins

def point_in_poly(x,y,poly):
    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside
	
	
def sort_hull_results(sarr,xarr,yarr):
	#Inputs: 
	#sarr: Xx2 array of simplices
	#xarr: Xx2 array of x-coordinates
	#yarr: Xx2 array of y-coordinates
	#Outputs:
	#xpts,ypts: 1D arrays of x,y-coordinates
	spts=[]
	xpts=[]
	ypts=[]
	spts.append(sarr[0][0])
	xpts.append(xarr[0][0])
	ypts.append(yarr[0][0])
	spts.append(sarr[0][1])
	xpts.append(xarr[0][1])
	ypts.append(yarr[0][1])


	while len(spts) < shape(sarr)[0]:
		for i in arange(shape(sarr)[0]-1)+1:
			if sarr[i][0] not in spts or sarr[i][1] not in spts:
				if sarr[i][0]==spts[-1]:
					spts.append(sarr[i][1])
					xpts.append(xarr[i][1])
					ypts.append(yarr[i][1])
				elif sarr[i][1]==spts[-1]:
					spts.append(sarr[i][0])
					xpts.append(xarr[i][0])
					ypts.append(yarr[i][0])
	spts.append(spts[0])
	xpts.append(xpts[0])
	ypts.append(ypts[0])
	return array(xpts),array(ypts)

def check_inside(x,y):
	points = rand(30, 2)*10.-5.
	hull = ConvexHull(points)

	slist=[]
	xlist=[]
	ylist=[]
	for simplex in hull.simplices:
		slist.append(simplex)
		xlist.append(points[simplex,0])
		ylist.append(points[simplex,1])

	sarr=array(slist)
	xarr=array(xlist)
	yarr=array(ylist)

	perim_x,perim_y=sort_hull_results(sarr,xarr,yarr)

	tf1=inside(x,y,perim_x,perim_y)
	#tf2=point_in_poly(x,y,npts)

	print tf1

#if tf1==False:
#	##Scatter Plot
	#app = QtGui.QApplication([])
	#mw = QtGui.QMainWindow()
	#mw.resize(800,800)
	#view = pg.GraphicsLayoutWidget()  # GraphicsView with GraphicsLayout inserted by default
	#mw.setCentralWidget(view)
	#mw.show()
	#mw.setWindowTitle('pyqtgraph example: ScatterPlot')
	#w1=view.addPlot()
	#s1=pg.ScatterPlotItem(pen='b')
	#s1.addPoints(x=xpts,y=ypts)
	#s2=pg.ScatterPlotItem(pen='r')
	#s2.addPoints(x=[x,x],y=[y,y])
	#s3=pg.ScatterPlotItem(pen='g')
	#s3.addPoints(x=pts[:,0],y=pts[:,1])
	#w1.addItem(s1)
	#w1.addItem(s3)
	#w1.addItem(s2)
	
#	##Line Plot
	win=pg.GraphicsWindow()
	pl=win.addPlot()
	pl.plot(perim_x,perim_y,pen='b')
	xx=[x-1,x+1]
	yy=[y,y]
	pl.plot(xx,yy)
	xx=[x,x]
	yy=[y-1,y+1]
	pl.plot(xx,yy)
	
#	##Actually Plots the thing
	if __name__=='__main__':
		import sys
		if(sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
			QtGui.QApplication.instance().exec_()
