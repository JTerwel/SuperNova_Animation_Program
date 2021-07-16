
#This program retrieves a ref im and difference images of requsted epochs from ZTF
#It then cuts out the requested region and puts it in an animation cycling through
#the images chronologically.
#Version 2 has a better way of handling objects close to the edge of the image.
#If centra pixel (where [ra,dec] is focussed on) is outside image, it is skipped in the animation
#If central pixel is at the edge of the image, it is padded.
#Version 3 has improved the image skipping condition to allow the skipping condition on a larger  area than just the middle 1x1
#Generator function is used to be more memory friendly, also should make selecting usable images easier
#wcs is used to make sure all images have the correct orientation
#Improved readability of plots
#Added small compass plot

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
import pandas as pd
from ztfquery import query
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy import time
from astropy.coordinates import SkyCoord
import astropy.units as u

def main(): #main function
	#General parameters one might want to change are before the line calling to download the data
	#Object parameters
	obj_name = 'ZTF18aazhwnh' #Name of object shown
	obj_ra = 233.011667 #RA in degrees
	obj_dec = 23.5200056 #dec in degrees
	obj_lc = '/Users/terwelj/Documents/ZTFData/marshal/fp_noise_test/ZTF18aazhwnh_offset_SNT_1e-08.csv' #Location of the light curve for this object
	lc_dates = [58250, 58800] #The region of mjd that will be shown
	#query parameters
	size = 0.01 #radius in which to search for obsrvations in degrees
	#Make sure all dates are in jd before making the sql query
	t = time.Time([58288, 28292, 58360, 58770], format='mjd').jd
	sql_query = 'fid=2 and (obsjd BETWEEN {} and {} or obsjd BETWEEN {} and {})'.format(t[0], t[1], t[2], t[3])
	download_data = False #If false, all required data is assumed to be stored locally already -> no downloading required or attempted
	#Image showing and animation parameters
	t_frame = 500 #amount of time each frame is visible in the animation in ms
	cutout_dims = 41 #tuple or int (single value generates square) best to have odd value(s)-> equal nr of pixels on each side
	min_visible_region = 9 #Like cutout_dims. controls size of region that has to be inside the image for it to be included in the animation. if 0 no check will be done.
	central_region = (5,5) #tuple (height, width), sets rectangle size around image centre and how many rows / colums to average over in central row / column plots
	save_name = 'ZTF18aazhwnh_offset.mp4' #Name under which the animation is saved
	#Download the data and list image locations
	refquery, difquery = download_images([obj_ra, obj_dec], size, sql_query, download_data)
	im_locs = refquery.get_local_data() + difquery.get_local_data('scimrefdiffimg.fits.fz')
	#Load lc data
	lc_dat = get_lc(obj_lc, lc_dates)
	#Make animation
	anim = animate_observations([(obj_ra, obj_dec)], obj_name, im_locs, lc_dat, t_frame, min_visible_region, cutout_dims, central_region)
	#save animation
	anim.save('/Users/terwelj/Movies/{}'.format(save_name))
	return

def animate_observations(obj_loc, obj_name, im_locs, lc_dat, t_frame, min_visible_region, cutout_dims, central_region):
	#Open 1st image so all plots can be initiated properly
	im0_full = fits.open(im_locs[0])[-1]
	wcs = WCS(im0_full.header)
	pix_xy = wcs.all_world2pix(obj_loc, 1)[0] #Pixel coords of obj
	im0 = Cutout2D(im0_full.data, pix_xy, cutout_dims, mode='partial', wcs=wcs)
	x_mesh, y_mesh = np.meshgrid(np.arange(0,np.shape(im0.data)[1]), np.arange(0,np.shape(im0.data)[0]))
	#Set figure and subplot locations
	gs = gridspec.GridSpec(60,140)
	fig = plt.figure(figsize=(14,6))
	ax_im = fig.add_subplot(gs[11:, :50]) #Image cutout
	ax_line_hor = fig.add_subplot(gs[:10, :50]) #line through central row of im
	ax_line_ver = fig.add_subplot(gs[11:, 50:60]) #line through central col of im
	ax_wire = fig.add_subplot(gs[:28, 65:100], projection='3d') #wireframe
	ax_text = fig.add_subplot(gs[:28, 101:]) #text plot
	ax_lc = fig.add_subplot(gs[31:, 67:]) #lightcurve
	ax_cb = fig.add_axes([0.1, 0.11, 0.015, 0.625])
	ax_ne = fig.add_subplot(gs[:10, 50:60]) #Arrows pointing N, E
	#Set static parts of the plots (what is seen in all frames)
	ax_im.set_xticklabels([])
	ax_im.set_yticklabels([])
	ax_im.set_xticks([int(np.shape(im0)[0]/2)])
	ax_im.set_yticks([int(np.shape(im0)[1]/2)])
	ax_line_hor.set_xticks([])
	ax_line_ver.set_yticks([])
	ax_wire.set_xticks([])
	ax_wire.set_yticks([])
	ax_wire.view_init(45, 225) #elev, azimut (45,225) points bottom left corner of image down
	ax_text.set_xticks([])
	ax_text.set_yticks([])
	ax_text.set_xlim(0, 1)
	ax_text.set_ylim(0, 1)
	ax_text.spines['top'].set_visible(False)
	ax_text.spines['bottom'].set_visible(False)
	ax_text.spines['left'].set_visible(False)
	ax_text.spines['right'].set_visible(False)
	for filt in ['g', 'r', 'i']:
		selection = lc_dat[lc_dat.obs_filter.str.contains(filt)]
		if filt=='i': #Now filt can be used to set the color of the points
			filt='y'
		ax_lc.scatter(selection[selection.mag==99].obsmjd, selection[selection.mag==99].upper_limit, color=filt, alpha=0.3, marker='v', s=4, zorder=1) #upper limits
		ax_lc.errorbar(selection[selection.mag!=99].obsmjd, selection[selection.mag!=99].mag, yerr=selection[selection.mag!=99].mag_err, fmt='o', color=filt, ms=4, zorder=2) #detections
	ax_lc.set_xlim(min(lc_dat.obsmjd)-1, max(lc_dat.obsmjd)+1)
	ax_lc.set_ylim(23, 15)
	ax_lc.set_xlabel('mjd')
	ax_lc.set_ylabel('mag')
	ax_ne.set_xlim(0, 1)
	ax_ne.set_ylim(0, 1)
	ax_ne.spines['top'].set_visible(False)
	ax_ne.spines['bottom'].set_visible(False)
	ax_ne.spines['left'].set_visible(False)
	ax_ne.spines['right'].set_visible(False)
	#Initiate objects holding the changing parts of the animation animation
	im = ax_im.imshow(im0.data, cmap='gray', origin='lower')
	cb = fig.colorbar(im, cax = ax_cb, ticklocation='left') #colorbar
	line_hor, = ax_line_hor.plot([], [], 'k')
	line_ver, = ax_line_ver.plot([], [], 'k')
	wire = ax_wire.plot_wireframe(x_mesh, y_mesh, im0.data, color='k', linewidth=0.5)
	text = ax_text.text(0, 1, ' ', va='top') #Position is locked here (lower left of text block)
	lc = ax_lc.scatter([], [], edgecolors='k', facecolors='none', linewidths=3, s=64, zorder=3) #Hollow circle to place over associated point if found
	lc_line = ax_lc.axvline(0, ls='--', color='k')
	movers = [im, line_hor, line_ver, wire, text, lc, lc_line]
	#animate and return animation
	anim = FuncAnimation(fig, animate_step, frames=get_next_im(im_locs, obj_loc, min_visible_region, cutout_dims), interval=t_frame, fargs=(movers, ax_im, ax_wire, ax_ne, lc_dat, central_region, obj_name), save_count=len(im_locs), blit=True, repeat=False)
	return anim

def animate_step(im_dat, movers, ax_im, ax_wire, ax_ne, lc_dat, central_region, obj_name):
	#im_dat is the following collection: [cutout, imtype, date, obsfilter]
	#Get values needed to plot this frame
	data, theta = find_orientation(im_dat[0]) #theta to be used for arrows to N, E
	vmin = np.nanmin(data) #min value of image ignoring nan
	vmax = np.nanmax(data) #max value of image ignoring nan
	central_bot = int(np.shape(data)[0]/2) - int(central_region[0]/2)
	central_top = central_bot + central_region[0]
	central_left = int(np.shape(data)[1]/2) - int(central_region[1]/2)
	central_right = central_left + central_region[1]
	points_hor = np.average(data[central_bot:central_top,:], axis=0) #Get averaged middle rows (middle rounded down)
	points_ver = np.average(data[:, central_left:central_right], axis=1) #Get averaged middle columns (rounded down)
	horver_minmax = [np.nanmin(points_hor), np.nanmin(points_ver), np.nanmax(points_hor), np.nanmax(points_ver)]
	for i in range(len(horver_minmax)): #Change unexpected output to 0
		if np.isnan(horver_minmax[i]) or np.isinf(horver_minmax[i]):
			horver_minmax[i]=0
	rect = Rectangle((central_left-0.5, central_bot-0.5), central_region[1], central_region[0], edgecolor='g', facecolor='none', lw=1.5) #rectangle around the middle of the object
	x_mesh, y_mesh = np.meshgrid(np.arange(0,np.shape(data)[1]), np.arange(0,np.shape(data)[0])) #make the meshgrid for the wireframe
	text = gen_im_text(obj_name, im_dat[1], im_dat[2], im_dat[3], np.nanmax(data), np.nanmean(data), np.nanstd(data)) #generate text for this image
	if im_dat[1]=='ref':
		lc_points = [[],[]]
	else:
		lc_points = find_lc_point(im_dat[2], im_dat[3], lc_dat)
	#Set the axes that have changed
	movers[0].set_clim(vmin, vmax)
	movers[1].axes.set_xlim(-0.5, len(points_hor)-0.5)
	movers[1].axes.set_ylim(horver_minmax[0]-1, horver_minmax[2]+1)
	movers[2].axes.set_xlim(horver_minmax[1]-1, horver_minmax[3]+1)
	movers[2].axes.set_ylim(-0.5, len(points_ver)-0.5)
	#Place the data in the plots
	movers[0].set_array(data)
	ax_im.patches = [] #remove old patch
	ax_im.add_patch(rect)
	movers[1].set_data(np.arange(0,len(points_hor)), points_hor)
	movers[2].set_data(points_ver, np.arange(0,len(points_ver)))
	ax_wire.collections.remove(movers[3])
	movers[3] = ax_wire.plot_wireframe(x_mesh, y_mesh, data, color='k', linewidth=0.5)
	movers[4].set_text(text)
	movers[5].set_offsets(lc_points)
	if im_dat[1]=='ref':
		movers[6].set_xdata([0, 0])
	else:
		movers[6].set_xdata([im_dat[2].mjd, im_dat[2].mjd])
	ax_ne.clear()
	ax_ne.arrow(0.5, 0.5, 0.3*np.cos(theta), 0.3*np.sin(theta), color='k', length_includes_head=True, head_width=0.1, head_length=0.1)
	ax_ne.arrow(0.5, 0.5, -0.3*np.sin(theta), 0.3*np.cos(theta), color='k', length_includes_head=True, head_width=0.1, head_length=0.1)
	ax_ne.text(0.5+0.4*np.cos(theta), 0.5+0.4*np.sin(theta), 'E', color='k', va='center', ha='center')
	ax_ne.text(0.5-0.4*np.sin(theta), 0.5+0.4*np.cos(theta), 'N', color='k', va='center', ha='center')
	ax_ne.set_xticks([])
	ax_ne.set_yticks([])
	return movers

def find_orientation(im): #Checks if im is mirrored, and finds direction of non-mirrored image east
	wcs = im.wcs
	pix00 = wcs.pixel_to_world(0,0)
	to_east = wcs.world_to_pixel(SkyCoord(ra=pix00.ra+1*u.arcmin, dec=pix00.dec, frame='icrs'))
	to_north = wcs.world_to_pixel(SkyCoord(ra=pix00.ra, dec=pix00.dec+1*u.arcmin, frame='icrs'))
	if (((to_east[0]>0) & (to_north[1]<0)) | ((to_east[0]<0) & (to_north[1]>0))):
		data = np.flipud(im.data)
		theta = np.arctan2(-to_east[1], to_east[0])
	else: #image is not mirrored
		data = im.data
		theta = np.arctan2(to_east[1], to_east[0])
	if theta<0: #convert to 0-360 degrees
		theta = (theta+2*np.pi)*180/np.pi
	else:
		theta = theta*180/np.pi
	#return image rotated by the required amount of 90 deg + residual theta in rad
	if theta%90>45: #Make sure east is not to far up
		theta += 90
		resid = (theta%90-90)*np.pi/180
	else:
		resid = (theta%90)*np.pi/180
	rots = int(theta/90)
	return np.rot90(data, rots, (1,0)), resid

def find_lc_point(date, im_filter, lc_dat): #Finds points in light curve that belong to the image
	good_dates = lc_dat[np.abs(lc_dat.obsmjd-date.mjd)<1e-4] #dates are within 8.64 s of each other
	if len(good_dates[good_dates.obs_filter==im_filter])==0: #no matches
		return [[],[]]
	vals = []
	for i in range(len(good_dates[good_dates.obs_filter==im_filter])):
		xval = good_dates[good_dates.obs_filter==im_filter].obsmjd.iloc[i]
		yval = np.nanmin(good_dates[good_dates.obs_filter==im_filter][['mag','upper_limit']].iloc[i])
		vals.append([xval,yval])
	return vals

def gen_im_text(name, im_type, date, obs_filter, peak_val, mean_val, standev): #Generates text lines that are specific for the currently shown image
	text = 'Name {}\nImage type: {}\nFilter: {}\n\nDate: {}\nmjd: {}\n\nPeak value: {}\nMean value: {}\nStandard deviation: {}'
	if im_type=='ref':
		return text.format(name, im_type, obs_filter, ' ', ' ', peak_val, mean_val, standev)
	else:
		greg = date.isot
		date_greg = '{}-{}-{} {}:{}:{}'.format(greg[8:10], greg[5:7], greg[:4], greg[11:13], greg[14:16], greg[17:])
		return text.format(name, im_type, obs_filter, date_greg, date.mjd, peak_val, mean_val, standev)

def get_next_im(im_locs, obj_loc, min_visible_region, cutout_dims):
	for i in im_locs:
		im = fits.open(i)[-1]
		wcs = WCS(im.header)
		pix_xy = wcs.all_world2pix(obj_loc, 1)[0] #Pixel coords of obj
		if ((min_visible_region != 0)&(min_visible_region !=(0,0))): #minimal visible region wanted
			try: #Should only succeed if all pixels in this cutout have values
				dummy = Cutout2D(im.data, pix_xy, min_visible_region, mode='strict')
			except:
				#print('Image at obsjd {} rejected due to NaN close to im center'.format(im.header['OBSJD']))
				continue #Start next i in loop
		cutout = Cutout2D(im.data, pix_xy, cutout_dims, mode='partial', wcs=wcs)
		#Get relevand image data not available in the cutout
		if 'ZTFDataref' in i:
			imtype = 'ref'
			date = ' '
			if im.header['DBFID']==1:
				obsfilter = 'ZTF_g'
			elif im.header['DBFID']==2:
				obsfilter = 'ZTF_r'
			else:
				obsfilter = 'ZTF_i'
		else:
			imtype = 'dif'
			date = time.Time(im.header['OBSJD'], format='jd')
			obsfilter = im.header['FILTER']
		#yield 1 item which is a list of required items
		yield [cutout, imtype, date, obsfilter]

def get_lc(lc_loc, lc_dates):
	lc = pd.read_csv(lc_loc, header=0, usecols=['obsmjd', 'filter', 'upper_limit', 'mag', 'mag_err'])
	lc.rename(columns={'filter':'obs_filter'}, inplace=True) #To avoid confusion
	return lc[(lc.obsmjd>=min(lc_dates)) & (lc.obsmjd<=max(lc_dates))] #Return only the wanted dates

def download_images(radec, size, sql_query, download_data): #Downloads the (meta-)data, and returns the meta-data
	#Get the meta-data and sort them to be in chronologically order
	refquery = query.ZTFQuery()
	refquery.load_metadata(kind='ref', radec=radec)
	difquery = query.ZTFQuery()
	difquery.load_metadata(radec=radec, size=size, sql_query=sql_query)
	difquery.metatable.sort_values(by=['obsjd'], inplace = True)
	#get the images if wanted
	if download_data:
		refquery.download_data()
		difquery.download_data('scimrefdiffimg.fits.fz')
	return refquery, difquery

if (__name__ == '__main__'):
	main()
