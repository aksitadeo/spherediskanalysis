import matplotlib.pyplot as plt
import numpy as np
import fnmatch
import os
from astropy.io import fits
from skimage.measure import EllipseModel
from matplotlib.patches import Ellipse
from scipy import interpolate
import cv2
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from skimage.measure import EllipseModel
from matplotlib.patches import Ellipse

def rotate_image(image, angle,xc,yc):
    #angle in deg
    #positive angle rotates image conter-clockwise but we are plotting image upside down(I don`t know why, it is caused by extent parameter in the imshow) and thats why it seems like positive angle corresponds to the clockwise direction. In this case image looks similar to plotted with ds9 after irdap.  
    image_center = (xc,yc)
    (h, w) = image.shape[:2]

    rot_mat = cv2.getRotationMatrix2D(image_center, angle, 1)
    result = cv2.warpAffine(image, rot_mat, (w, h))#,flags=cv2.INTER_LINEAR)
    return result


def rotate_image0(img, ang,xc,yc):

    rows,cols= img.shape

    # Create the transformation matrix
    angle = np.radians(-ang)
    x0, y0 = xc,yc
    M = np.array([[np.cos(angle), -np.sin(angle), x0*(1-np.cos(angle))+ y0*np.sin(angle)],
                  [np.sin(angle), np.cos(angle), y0*(1-np.cos(angle))- x0*np.sin(angle)]])
    # get the coordinates in the form of (0,0),(0,1)...
    # the shape is (2, rows*cols)
    orig_coord = np.indices((cols, rows)).reshape(2,-1)
    # stack the rows of 1 to form [x,y,1]
    orig_coord_f = np.vstack((orig_coord, np.ones(rows*cols)))
    transform_coord = np.dot(M, orig_coord_f)
    # Change into int type
    transform_coord = transform_coord.astype(np.int)
    # Keep only the coordinates that fall within the image boundary.
    indices = np.all((transform_coord[1]<rows, transform_coord[0]<cols, transform_coord[1]>=0, transform_coord[0]>=0), axis=0)
    # Create a zeros image and project the points
    img1 = np.zeros_like(img)
    img1[transform_coord[1][indices], transform_coord[0][indices]] = img[orig_coord[1][indices], orig_coord[0][indices]]
    # Display the image
    #out = cv2.hconcat([img,img1])
    #cv2.imshow('a2',out)
    #cv2.waitKey(0)
    return img1




def rotate_points_el(x,y,angle,xc,yc):
    rotx=xc+(x-xc)*np.cos(angle)-(y-yc)*np.sin(angle)
    roty=yc+(y-yc)*np.cos(angle)+(x-xc)*np.sin(angle)
    return rotx,roty


def Loadimages(star,fittype,annulus,dirdat):
#This function download fits files.
#star - star name that used as a folder name for finding the right way
#fittype - type of file (Q_phi, PI, etc)

    dir =dirdat + star +'/deconvolution/deconvolved_'+fittype+'/'
    ifile = '*_decon.fits'
    files = os.listdir(dir)
    for file in files:
        if fnmatch.fnmatch(file, ifile):
            hduli = fits.open(dir + file)
            image = hduli['Primary'].data 
            
    return image,hduli

    
def cart2polar(x, y,xc,yc,startangle):
    r = np.sqrt((x-xc)**2 + (y-yc)**2)
    theta = np.arctan2((y-yc), (x-xc))-startangle
    for i in range(0,len(theta)):
        if theta[i]<0:
            theta[i]=2*np.pi+theta[i]
        

    return r, theta

def cart2polar_for_mask_defining(x, y):
    r = np.sqrt((x)**2 + (y)**2)
    theta = np.arctan2(y, x)


    return r, theta


def northeast_pix(lim,xc,yc):
    x_ar=(xc+lim-2)
    y_ar=(yc-lim+2)
    plt.arrow(x_ar, y_ar, 0, 3,color='white',width=0.1, length_includes_head=True)
    plt.text(x_ar-0.5,y_ar+3.5,'N',fontsize='small',color="white")        
    plt.arrow(x_ar, y_ar, -3, 0, color='white',width=0.1, length_includes_head=True)
    plt.text(x_ar-4,y_ar-0.3,'E',fontsize='small',color="white") 

def northeast_rotated_pix(angle1,xc,yc,lim): 
    angle=np.deg2rad(-angle1)              
    x_ar=xc+(lim-7)
    y_ar=(yc+-lim)+7
    
    x_len=-(3)*np.sin(angle)
    y_len=(3)*np.cos(angle)
    x_len_t=-(4.5)*np.sin(angle)
    y_len_t=(4.5)*np.cos(angle)
    
    plt.arrow(x_ar, y_ar, x_len, y_len, color='white',width=0.1, length_includes_head=True)
    plt.text(x_ar+x_len_t,y_ar+y_len_t ,'N',fontsize='small',color="white") 
    x_len=(-3)*np.cos(angle)
    y_len=(-3)*np.sin(angle)
    x_len_t=(-4.5)*np.cos(angle)
    y_len_t=(-4.5)*np.sin(angle)   
    plt.arrow(x_ar, y_ar, x_len, y_len, color='white',width=0.1, length_includes_head=True)

    plt.text(x_ar+x_len_t,y_ar+y_len_t ,'E',fontsize='small',color="white")   

def northeast_rotated_depr_pix(angle1,xc,yc,lim,cosi): 
    angle=np.deg2rad(-angle1)              
    x_ar=xc+(lim/cosi-7)
    y_ar=(yc+-lim)/cosi+7
    
    x_len=-(3)*np.sin(angle)
    y_len=(3)*np.cos(angle)/cosi
    x_len_t=-(4.5)*np.sin(angle)
    y_len_t=(4.5)*np.cos(angle)/cosi
    
    plt.arrow(x_ar, y_ar, x_len, y_len, color='white',width=0.1, length_includes_head=True)
    plt.text(x_ar+x_len_t,y_ar+y_len_t ,'N',fontsize='small',color="white") 
    x_len=(-3)*np.cos(angle)
    y_len=(-3)*np.sin(angle)/cosi
    x_len_t=(-4.5)*np.cos(angle)
    y_len_t=(-4.5)*np.sin(angle)/cosi      
    plt.arrow(x_ar, y_ar, x_len, y_len, color='white',width=0.1, length_includes_head=True)

    plt.text(x_ar+x_len_t,y_ar+y_len_t ,'E',fontsize='small',color="white")   



def northeast2(lim,ps):
    x_ar=(lim-2) * ps
    y_ar=(-lim+2) * ps
    plt.arrow(x_ar, y_ar, 0, 3*ps, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)
    plt.text(x_ar-0.5 * ps,y_ar+3.5 * ps,'N',fontsize='small',color="white")        
    plt.arrow(x_ar, y_ar, -3*ps, 0, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)
    plt.text(x_ar-4 * ps,y_ar-0.5 * ps,'E',fontsize='small',color="white") 


def northeast(lim,ps,ax):
    x_ar=(lim-2) * ps
    y_ar=(-lim+2) * ps
    ax.arrow(x_ar, y_ar, 0, 3*ps, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)
    ax.text(x_ar-0.5 * ps,y_ar+3.5 * ps,'N',fontsize='small',color="white")        
    ax.arrow(x_ar, y_ar, -3*ps, 0, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)
    ax.text(x_ar-4 * ps,y_ar-0.5 * ps,'E',fontsize='small',color="white") 
   

def northeast_rotated(angle1,xc,yc,lim,ps): 
    angle=np.deg2rad(-angle1)              
    x_ar=(lim-5) * ps
    y_ar=(-lim+5) * ps
    
    x_len=-(3)*np.sin(angle)*ps
    y_len=(3)*np.cos(angle)*ps
    x_len_t=-(4.5)*np.sin(angle)*ps
    y_len_t=(4.5)*np.cos(angle)*ps
    
    plt.arrow(x_ar, y_ar, x_len, y_len, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)
    plt.text(x_ar+x_len_t,y_ar+y_len_t ,'N',fontsize='small',color="white") 
    x_len=(-3)*np.cos(angle)*ps
    y_len=(-3)*np.sin(angle)*ps
    x_len_t=(-4.5)*np.cos(angle)*ps
    y_len_t=(-4.5)*np.sin(angle)*ps   
    plt.arrow(x_ar, y_ar, x_len, y_len, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)
    plt.text(x_ar+x_len_t,y_ar+y_len_t ,'E',fontsize='small',color="white")   

def northeast_rotated_depr(angle1,xc,yc,lim,cosi,ps): 
    angle=np.deg2rad(-angle1)              
    x_ar=(lim-5) * ps
    y_ar=(-lim+5) * ps
    
    x_len=-(3)*np.sin(angle)*ps
    y_len=(3)*np.cos(angle)/cosi*ps
    x_len_t=-(4.5)*np.sin(angle)*ps
    y_len_t=(4.5)*np.cos(angle)/cosi*ps
    
    plt.arrow(x_ar, y_ar, x_len, y_len, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)
    plt.text(x_ar+x_len_t,y_ar+y_len_t ,'N',fontsize='small',color="white") 
    x_len=(-3)*np.cos(angle)*ps
    y_len=(-3)*np.sin(angle)/cosi*ps
    x_len_t=(-4.5)*np.cos(angle)*ps
    y_len_t=(-4.5)*np.sin(angle)/cosi*ps      
    plt.arrow(x_ar, y_ar, x_len, y_len, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)

    plt.text(x_ar+x_len_t,y_ar+y_len_t ,'E',fontsize='small',color="white")   

def northeast_rotated_depr_ax(angle1,xc,yc,lim,cosi,ps,ax): 
    angle=np.deg2rad(-angle1)              
    x_ar=(lim-5) * ps
    y_ar=(-lim+5) * ps
    
    x_len=-(3)*np.sin(angle)*ps
    y_len=(3)*np.cos(angle)/cosi*ps
    x_len_t=-(4.5)*np.sin(angle)*ps
    y_len_t=(4.5)*np.cos(angle)/cosi*ps
    
    ax.arrow(x_ar, y_ar, x_len, y_len, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)
    ax.text(x_ar+x_len_t,y_ar+y_len_t ,'N',fontsize='small',color="white") 
    x_len=(-3)*np.cos(angle)*ps
    y_len=(-3)*np.sin(angle)/cosi*ps
    x_len_t=(-4.5)*np.cos(angle)*ps
    y_len_t=(-4.5)*np.sin(angle)/cosi*ps      
    ax.arrow(x_ar, y_ar, x_len, y_len, color='white',width=1,head_width=10, length_includes_head=True,head_length=15)

    ax.text(x_ar+x_len_t,y_ar+y_len_t ,'E',fontsize='small',color="white")   


def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)


def plottingroutinemas(image,lim,ps,n,star,ax):
    max = np.max(image)
    min=np.min(image)
    d = (n-1) * ps / 2

    plt.imshow(image, vmin=min, vmax=max, extent=(-d, d, d, -d))
    #plt.plot(0, 0, "+", color="white")
    plt.xlim(-lim * ps, lim * ps)
    plt.ylim(-lim * ps, lim * ps)
    plt.xlabel('mas')
    plt.ylabel("mas")
    plt.colorbar()
    plt.tight_layout      
    northeast(lim,ps,ax)
    scale_mas(star,ax)


def plotImage(image, lim):
    n = image.shape[0]
    
    fig, ax = plt.subplots()
    max = np.max(image)
    min=np.min(image)
    ps = 12.27 #mas per pixel for IRDIS
    d = n * ps / 2
    plt.imshow(image, vmin=min, vmax=max, extent=(-d, d, d, -d))
    #plt.plot(0, 0, "+", color="red")
    plt.xlim(-lim * ps, lim * ps)
    plt.ylim(-lim * ps, lim * ps)
    plt.colorbar()
    plt.tight_layout
   
def cutimage(image,cutlimit):
    n=image.shape[0]
    adc=int(cutlimit)
    bdc=int(n-cutlimit)
    image= image[adc:bdc,adc:bdc]
    return image


def plotImages(image1, image2, lim):
    n = image1.shape[0]

    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,6))
    image1 = np.arcsinh(image1)
    image2 = np.arcsinh(image2)
    max1 = np.max(image1)
    max2 = np.max(image2)
    min1 = np.min(image1)
    min2 = np.min(image2)
    ps = 12.27
    d = n*ps/2
    im1=ax1.imshow(image1, vmin=min1, vmax = max1, extent=(-d, d, d, -d) )
    im2=ax2.imshow(image2, vmin=min2, vmax = max2, extent=(-d, d, d, -d) )
    ax1.set_xlim(-lim*ps, lim*ps)
    ax1.set_ylim(-lim*ps, lim*ps)
    ax2.set_xlim(-lim*ps, lim*ps)
    ax2.set_ylim(-lim*ps, lim*ps)
    
    ax1.set_xlabel('mas')
    ax1.set_ylabel("mas")
    ax2.set_xlabel('mas')
    ax2.set_ylabel("mas")
    plt.tight_layout(pad=3.0)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax, orientation='vertical')    
    #plt.show()
    #plt.close()
def plotImages3(image1, image2,image3, lim):
    n = image1.shape[0]
    fig=plt.figure(figsize=(12, 12))
    G = gridspec.GridSpec(2, 2)
    ax1 = plt.subplot(G[0, :])
    ax2=plt.subplot(G[1, 0])
    ax3=plt.subplot(G[1, 1])
    image1 = np.arcsinh(image1)
    image2 = np.arcsinh(image2)
    image3 = np.arcsinh(image3)
    max1 = np.max(image1)
    max2 = np.max(image2)
    min1 = np.min(image1)
    min2 = np.min(image2)
    max3 = np.max(image3)
    min3 = np.min(image3)
    ps = 12.27
    d = n*ps/2
    im1=ax1.imshow(image1, vmin=min1, vmax = max1, extent=(-d, d, d, -d) )
    im2=ax2.imshow(image2, vmin=min2, vmax = max2, extent=(-d, d, d, -d) )
    im3=ax3.imshow(image3, vmin=min3, vmax = max3, extent=(-d, d, d, -d) )
    ax3.set_xlim(-lim*ps, lim*ps)
    ax3.set_ylim(-lim*ps, lim*ps)
    
    ax1.set_xlim(-lim*ps, lim*ps)
    ax1.set_ylim(-lim*ps, lim*ps)
    ax2.set_xlim(-lim*ps, lim*ps)
    ax2.set_ylim(-lim*ps, lim*ps)
    ax1.set_title('original_image')
    ax2.set_title('image*squared separation in pix')
    ax3.set_title('image*squared separation in arcsec') 

    plt.tight_layout(pad=3.0)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    fig.colorbar(im2, cax=cax, orientation='vertical')    
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im3, cax=cax, orientation='vertical')    
       
    #plt.show()
    #plt.close()

def plotImagesdepr(image1, image2, lim,cosi,xc,yc):
    n = image1.shape[0]

    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,6))
    #image1 = np.arcsinh(image1)
    #image2 = np.arcsinh(image2)
    max1 = np.max(image1)
    max2 = np.max(image2)
    min1 = np.min(image1)
    min2 = np.min(image2)
    ps = 12.27
    d = n*ps/2
    im1=ax1.imshow(image1, vmin=min1, vmax = max1, extent=(-d, d, d/cosi, -d/cosi) )
    im2=ax2.imshow(image2, vmin=min2, vmax = max2, extent=(-d, d, d/cosi, -d/cosi) )
    ax1.set_xlim(-lim*ps, lim*ps)
    ax1.set_ylim(-lim*ps, lim*ps)
    ax2.set_xlim(-lim*ps, lim*ps)
    ax2.set_ylim(-lim*ps, lim*ps)
    plt.tight_layout(pad=3.0)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax, orientation='vertical') 


     
    
def scale_mas(star,ax):


    starsdict = {'iras08544-4431_calib':1470,'hr4049':1574,'iras15469-5311':3179,'iras17038-4815':4330,'iw_car_calib':1811,'ru_cen_calib':1822,'u_mon_2019-01-03_calib':1067,'u_mon_2019-01-14_calib':1067,'u_mon_combined':1067,'ac_her':1231}
    """
    Draw a horizontal bar with length of in data coordinates,
    with a fixed label underneath.
    """
    #lenghtau=50
    #lenghtmas=lenghtau*3600000/np.pi*180/starsdict[star]/206264.806
    
    #from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    #asb = AnchoredSizeBar(ax.transData,
     #                     lenghtmas,
     #                     r"50 au",
     #                     loc='lower center',
     #                     pad=0.1, borderpad=0.5, sep=5,
     #                     frameon=False, color='white')
    #ax.add_artist(asb)
                        


def boarder_points_on_a_line(xst,yst,xend,yend,length):
    k=(yend-yst)/(xend-xst)
    topup=(length-math.dist((xst,yst),(xend,yend)))/2
    dx=(topup/np.sqrt(1+k*k))
    dy=k*dx
    if xst>xend:
        dx=-dx
        dy=-dy
    x1=xst-dx
    x2=xend+dx
    y1=yst-dy
    y2=yend+dy

    return x1,y1,x2,y2
  

def ellipsefit(image,pointvalue,shift,ps):
    
    data_p= np.where(image >= pointvalue)
    data_p1=np.vstack((data_p[1],data_p[0])).transpose()    
    a_points = data_p1
    x = a_points[:, 0]
    y = a_points[:, 1]

    #fiting an ellipse
    ell = EllipseModel()
    ell.estimate(a_points)
    xc, yc, a, b, theta = ell.params
    #calculating points of ellipse based on the previously defined parameters
    el_points=ell.predict_xy(np.linspace(0, 2 * np.pi, 50),params=(xc,yc,a,b,theta)) #without scaling into mas, only pixels` numbers
    el_points_mas=ell.predict_xy(np.linspace(0, 2 * np.pi, 50),params=((xc-shift)*ps,(yc-shift)*ps,a*ps,b*ps,theta)) #in mas
    sigma_2=np.sum(ell.residuals(a_points)**2)
    return xc, yc, a, b, theta,el_points,el_points_mas,sigma_2,x,y,ell

def significance_array_aolp(aolp,R):
    n = aolp.shape[0]
    lim=n/2
    ps = 12.27 #mas per pixel for IRDIS
    d = n * ps / 2
        
    
    phi = (aolp+90)    
    critarray=np.zeros_like(phi)
    for ix in range (0,n):
        for iy in range(0,n):
            if phi[ix,iy]>180:
                phi[ix,iy]=phi[ix,iy]-180
            
            
            
    for ix in range (490+2,n-2-490):
        for iy in range(490+2,n-2-490):
            if R[ix,iy]>=1:            
                datapix=[]
                for (iix,iiy) in [(ix,iy),(ix-1,iy),(ix+1,iy),(ix,iy-1),(ix,iy+1)]: 
                #for iix in range (ix-1,ix+1):
                    #for iiy in range(iy-1,iy+1):
                    if R[iix,iiy]>=1:
                        datapix.append(abs(phi[iix,iiy]-90))
                        
                
                crit=np.std(datapix)               
                critarray[ix,iy]=crit
                print('working'+str(ix)+str(iy))
    for ix in range (0,n):
        for iy in range(0,n):
            if critarray[ix,iy]==0:
                critarray[ix,iy]=np.max(critarray)
    medianstd=np.nanmedian(critarray[490+2:n-2-490,490+2:n-2-490])
    return medianstd, critarray


def arc_fitting(star,imagefull,ps):
   

    n = imagefull.shape[0]

    #Creating grid for the multiplying the image by the separation from star.
    xr = np.linspace(-n/2, n/2, num=n)
    yr = np.linspace(-n/2, n/2, num=n)
    x0 = 0.5
    y0 = 0.5
    xr = xr-x0
    yr = yr-y0
    Xr, Yr = np.meshgrid(xr, yr)
    R = np.sqrt(Xr**2 + Yr**2)
    if star=='iras08544-4431_calib' or star=='iras15469-5311' or star=='iw_car_calib' or star=='ac_her':
        r, pos_angle = cart2polar_for_mask_defining(Xr, Yr)
        pos_angle=np.rad2deg(pos_angle)+180
    if star=='iras08544-4431_calib':
        image_arc = imagefull*(R<16)*(R>6)*(Xr>-10)*(Xr<10)*(Yr<9)*(pos_angle<270)*(pos_angle>35)
        q_lim=1.1
        image_arc2 = imagefull*(R<15)*(R>6)*(pos_angle<280)*(pos_angle>220)

    if star=='iras15469-5311':
        Xr15, Yr15 = np.meshgrid(xr+1, yr+1)
        R15 = np.sqrt(Xr15**2 + Yr15**2)
        image_arc = imagefull*(R15<12)*(R15>5.5)*(imagefull<10)
        q_lim=4
    if star=='iw_car_calib':
        image_arc = imagefull*(R<14)*(R>5)#(Xr<-4)
        q_lim=2#1.5
    if star=='ac_her':
        image_arc = imagefull*(R<7)*(R>4)+imagefull*(R>3)*(R<6)*(pos_angle<60)+imagefull*(R>2)*(R<5)*(pos_angle>200)+imagefull*(R>5)*(R<15)*(pos_angle>90)*(pos_angle<180)#(Xr<-4)#*(R>5)
        q_lim=4


    #selecting points that will be used for the ellipse fitting 
    data_p= np.where(image_arc >= q_lim)
    data_p1=np.vstack((data_p[1],data_p[0])).transpose()    
    a_points_arc = data_p1
    x_arc = a_points_arc[:, 0]
    y_arc = a_points_arc[:, 1]

    #fiting an ellipse
    ell_arc = EllipseModel()
    ell_arc.estimate(a_points_arc)
    xc_arc, yc_arc, a_arc, b_arc, theta_arc = ell_arc.params

    shift=n/2-0.5
    #calculating points of ellipse based on the previously defined parameters
    points_arc=ell_arc.predict_xy(np.linspace(0, 2 * np.pi, 50),params=(xc_arc,yc_arc,a_arc,b_arc,theta_arc)) #without scaling into mas, only pixels` numbers
    xc_arc2, yc_arc2, a_arc2, b_arc2, theta_arc2,points_arc2, points_mas_arc2 =(0,0,0,0,0,0,0)
    if star=='iras08544-4431_calib':
        #points_mas_arc=ell_arc.predict_xy(np.linspace(1.2* np.pi, 2.1 * np.pi, 50),params=((xc_arc-shift)*ps,(yc_arc-shift)*ps,a_arc*ps,b_arc*ps,theta_arc)) #in mas
        points_arc=ell_arc.predict_xy(np.linspace(0.7* np.pi, 2.1 * np.pi, 50),params=(xc_arc,yc_arc,a_arc,b_arc,theta_arc)) #without scaling into mas, only pixels` numbers
        points_mas_arc=ell_arc.predict_xy(np.linspace(0.7* np.pi, 2.1 * np.pi, 50),params=((xc_arc-shift)*ps,(yc_arc-shift)*ps,a_arc*ps,b_arc*ps,theta_arc)) #in mas
        
        #selecting points that will be used for the second arc fitting
        data_p_arc2= np.where(image_arc2 >= q_lim)
        data_p2=np.vstack((data_p_arc2[1],data_p_arc2[0])).transpose()    
        a_points_arc2 = data_p2
        x_arc2 = a_points_arc2[:, 0]
        y_arc2 = a_points_arc2[:, 1]

        #fiting 
        ell_arc2 = EllipseModel()
        ell_arc2.estimate(a_points_arc2)
        xc_arc2, yc_arc2, a_arc2, b_arc2, theta_arc2 = ell_arc2.params
        #calculating points of arc2 based on the previously defined parameters
        points_arc2=ell_arc2.predict_xy(np.linspace(0.8* np.pi, 2.1* np.pi, 50),params=(xc_arc2,yc_arc2,a_arc2,b_arc2,theta_arc2)) #without scaling into mas, only pixels` numbers
        points_mas_arc2=ell_arc2.predict_xy(np.linspace(0.8* np.pi, 2.1* np.pi, 50),params=((xc_arc2-shift)*ps,(yc_arc2-shift)*ps,a_arc2*ps,b_arc2*ps,theta_arc2)) #in mas

    if star=='iras15469-5311':
        points_mas_arc=ell_arc.predict_xy(np.linspace(0* np.pi, 2 * np.pi, 50),params=((xc_arc-shift)*ps,(yc_arc-shift)*ps,a_arc*ps,b_arc*ps,theta_arc)) #in mas
    if star=='iw_car_calib':
        points_mas_arc=ell_arc.predict_xy(np.linspace(0* np.pi, 2 * np.pi, 50),params=((xc_arc-shift)*ps,(yc_arc-shift)*ps,a_arc*ps,b_arc*ps,theta_arc)) #in mas
    if star=='ac_her':
        #points_mas_arc=ell_arc.predict_xy(np.linspace(0* np.pi, 2 * np.pi, 50),params=((xc_arc-shift)*ps,(yc_arc-shift)*ps,a_arc*ps,b_arc*ps,theta_arc)) #in mas
        a_arc=70
        b_arc=40
        points_mas_arc=ell_arc.predict_xy(np.linspace(0* np.pi, 2 * np.pi, 50),params=(2,5,a_arc,b_arc,1.33*np.pi)) #in mas
        xc_arc=2/ps+shift
        yc_arc=5/ps+shift
        points_arc=ell_arc.predict_xy(np.linspace(0* np.pi, 2 * np.pi, 50),params=(xc_arc,yc_arc,a_arc,b_arc,1.33*np.pi)) 
        
        cosi_arc=b_arc/a_arc
        incl_arc= np.rad2deg(np.arccos(cosi_arc))
        ecc_arc=np.sqrt(1-b_arc**2/a_arc**2)
        print('incl %f, eccentr %f' %(incl_arc,ecc_arc))

        sigma_2_arc=np.sum(ell_arc.residuals(a_points_arc)**2)
    
    return xc_arc, yc_arc, a_arc, b_arc, theta_arc, points_arc, points_mas_arc, xc_arc2, yc_arc2, a_arc2, b_arc2, theta_arc2,points_arc2, points_mas_arc2
        


