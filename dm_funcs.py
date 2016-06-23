import numpy as np
from scipy.misc import imread


def load_images():
	""" Load the HST images and make them presentable."""
	im1 = imread('lens1.jpg')
	im1 = [ im1[:,:,2], im1 ]
	im1[0][300:450,400:550] = im1[0].mean()
	im1[0]=im1[0][100:-100,100:-100]
	im1[1] = im1[1][100:-100,100:-100,:]

	im2 = imread('lens2.jpg')
	im2 = [ im2[500:-500,500:-500,0], im2[500:-500,500:-500,:] ]

	im3 = imread('lens3.jpg')
	im3 = [ im3[200:500,200:600,2], im3[200:500,200:600,:] ]
	im3[0][100:225,125:250] = im3[0][:50,:50].mean()

	im4 = imread('lens4.jpg')
	im4 = [ im4[:,:,2], im4 ]

	return {'Exercise 1':im1,'Exercise 2':im2,'Exercise 3':im3,'Exercise 4':im4}

def create_galaxy(x,y,**kwargs):
	""" Create a mock galaxy using a 2D Gaussian. """
	posangle = kwargs['posangle']
	xcen = kwargs['xcen']
	ycen = kwargs['ycen']
	peak = kwargs['peak']
	sigma = kwargs['sigma']
	axratio = kwargs['axratio']

	phi = np.deg2rad(posangle)
	x1 = (x - xcen)* np.cos(phi) + ( y- ycen)*np.sin(phi)
	y1 = (y-ycen)*np.cos(phi) - (x-xcen)*np.sin(phi)

	return peak * np.exp( -.5*( x1*x1*axratio + y1*y1/axratio)/(sigma*sigma))


def load_grid(Npoints):
	""" Create a 2D grid of points. """
	nx = Npoints
	ny = Npoints
	xlims = (-2.5,2.5)
	ylims = (-2.5,2.5)

	x = (xlims[1] - xlims[0]) * np.outer(np.ones(ny),np.arange(nx))/float(nx-1) + xlims[0]
	y = (ylims[1] - ylims[0]) * np.outer(np.arange(ny),np.ones(nx))/float(ny-1) + ylims[0]

	return x,y


def compute_deflection(x,y,**kwargs):
	""" Compute the deflection and shear from a singular isothermal ellipsoid potential """
	mass = kwargs['Mass']
	xcen = kwargs['xcen']
	ycen = kwargs['ycen']
	axratio = np.abs(kwargs['axratio'])
	posangle = kwargs['posangle']

	tol = .001

	if axratio>1:
		axratio = 1.0/axratio
		posangle += 90
	phi = np.deg2rad(posangle)

	x_new = (x-xcen)*np.cos(phi) + (y-ycen)*np.sin(phi)
	y_new = (y-ycen)*np.cos(phi) - (x-xcen)*np.sin(phi)

	r = np.sqrt(axratio*x_new*x_new + y_new * y_new /axratio)

	q = np.sqrt(1./axratio - axratio)
    softening = (r==0).astype(float)
    fac_x = x_new/(r + softening)
	fac_y = y_new/(r + softening)



	if q>=tol:
		x_def = (mass/q)*np.arctan(q*fac_x)
		y_def = (mass/q)*np.arctanh(q*fac_y)
	else:
		x_def = mass*fac_x
		y_def = mass*fac_y

	x_out = x_def *np.cos(phi) - y_def*np.sin(phi)
	y_out = y_def * np.cos(phi) + x_def * np.sin(phi)

	return x_out,y_out

def draw_ellipse(scale=1,xc=0,yc=0,r=1,phi=0):
	"""" Draw an ellipse given some orientation and axis ratio. """
	phi = np.deg2rad(phi)
	theta = np.linspace(0,2*np.pi,100)
	xe =  np.cos(theta)*np.cos(phi) - r*np.sin(theta)*np.sin(phi)
	ye =  np.cos(theta)*np.sin(phi) + r*np.sin(theta)*np.cos(phi)
	xe = xc + xe/scale
	ye = yc + ye/scale
	return xe,ye
