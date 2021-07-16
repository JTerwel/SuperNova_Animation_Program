#This program retrieves a ref im and difference images of requsted epochs from ZTF
#It then cuts out the requested region and puts it in an animation cycling through
#the images chronologically.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
import pandas as pd
from ztfquery import query
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy import time

def main(): #main function
	#General parameters one might want to change are before the line calling to download the data
	#Object parameters
	#obj_name = 'ZTF18abmxfrc' #Name of object shown
	#obj_ra = 252.7534888 #RA in degrees
	#obj_dec = 25.8759437 #dec in degrees
	#obj_lc = '/Users/terwelj/Documents/ZTFData/marshal/SN_Ia/ZTF18abmxfrc_SNT_1e-08.csv' #Location of the light curve for this object
	#lc_dates = [58300, 58800] #The region of mjd that will be shown
	obj_name = 'ZTF18aasdted'
	obj_ra = 261.4131984
	obj_dec = 59.4467859
	obj_lc = '/Users/terwelj/Documents/ZTFData/marshal/SN_Ia/ZTF18aasdted_SNT_1e-08.csv'
	lc_dates = [58200, 58850]
	#query parameters
	size = 0.01 #radius in which to search for obsrvations in degrees
	#RELOCATE save_name!
	#save_name = 'ZTF18abmxfrc_late.mp4' #Name under which the animation is saved
	save_name = 'ZTF18aasdted_late.mp4'
	#Make sure all dates are in jd before making the sql query
	#t = time.Time([58503, 58504, 58543, 58544, 58577, 58578], format='mjd').jd
	#sql_query = 'fid=2 and (obsjd BETWEEN {} and {} or obsjd BETWEEN {} and {} or obsjd BETWEEN {} and {})'.format(t[0], t[1], t[2], t[3], t[4], t[5])
	t = time.Time([58250, 58280, 58700, 58730], format='mjd').jd
	sql_query = 'fid=2 and (obsjd BETWEEN {} and {} or obsjd BETWEEN {} and {})'.format(t[0], t[1], t[2], t[3])
	download_data = True #If false, all required data is assumed to be stored locally already -> no downloading required or attempted
	#Image showing and animation parameters
	t_frame = 500 #amount of time each frame is visible in the animation in ms
	cutout_dims = 41 #tuple or int (single value generates square) best to have odd value(s)-> equal nr of pixels on each side
	#Download the data
	refquery, difquery = download_images([obj_ra, obj_dec], size, sql_query, download_data)
	#Get cutout of object in each image + stats of each image
	cutouts, im_dat = gen_cutouts(cutout_dims, [(obj_ra, obj_dec)], refquery, difquery, obj_name)
	lc_dat = get_lc(obj_lc, lc_dates)
	#Make animation
	anim = animate_observations(cutouts, im_dat, lc_dat, t_frame)
	#save animation
	anim.save('/Users/terwelj/Movies/{}'.format(save_name))
	return

def animate_observations(cutouts, im_dat, lc_dat, t_frame):
	x_mesh, y_mesh = np.meshgrid(np.arange(0, np.shape(cutouts[0].data)[1]), np.arange(0, np.shape(cutouts[0].data)[0])) #Meshgrid for wireframe, for some reason the order of x,y needs to be reversed to get x_mesh and y_mesh dimensions correct (not sure why, but it works for now)
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
	#Set static parts of the plots (what is seen in all frames)
	ax_im.set_xticklabels([])
	ax_im.set_yticklabels([])
	ax_im.set_xticks([int(np.shape(cutouts[0])[0]/2)])
	ax_im.set_yticks([int(np.shape(cutouts[0])[1]/2)])
	ax_line_hor.set_xticks([])
	ax_line_ver.set_yticks([])
	ax_wire.set_xticks([])
	ax_wire.set_yticks([])
	ax_wire.view_init(45, 45) #elev, azimut (45,45) points bottom left corner of image down
	ax_wire.invert_xaxis()
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
	ax_lc.set_ylim(23, 15)
	ax_lc.set_xlabel('mjd')
	ax_lc.set_ylabel('mag')
	#Initiate objects holding the changing parts of the animation animation
	im = ax_im.imshow(cutouts[0].data, cmap='gray')
	cb = fig.colorbar(im, cax = ax_cb, ticklocation='left') #colorbar
	line_hor, = ax_line_hor.plot([], [], 'k')
	line_ver, = ax_line_ver.plot([], [], 'k')
	wire = ax_wire.plot_wireframe(x_mesh, y_mesh, cutouts[0].data, color='k', linewidth=0.5)
	text = ax_text.text(0, 1, ' ', va='top') #Position is locked here (lower left of text block)
	lc = ax_lc.scatter([], [], edgecolors='k', facecolors='none', linewidths=3, s=64, zorder=3) #Hollow circle to place over associated point if found
	movers = [im, line_hor, line_ver, wire, text, lc]
	#animate and return animation
	nr_frames = np.shape(cutouts)[0]
	anim = FuncAnimation(fig, animate_step, frames=nr_frames, interval=t_frame, fargs=(cutouts, im_dat, lc_dat, movers, ax_wire), blit=True)
	return anim

def animate_step(i, cutouts, im_dat, lc_dat, movers, ax_wire):
	current_cutout = cutouts[i].data #Get the current image to show
	vmin = np.nanmin(current_cutout) #min value of image ignoring nan
	vmax = np.nanmax(current_cutout) #max value of image ignoring nan
	points_hor = current_cutout[int(np.shape(current_cutout)[0]/2),:] #Get the middle row (middle rounded down)
	points_ver = current_cutout[:,int(np.shape(current_cutout)[1]/2)] #Get the middle column (rounded down)
	x_mesh, y_mesh = np.meshgrid(np.arange(0,np.shape(current_cutout)[1]), np.arange(0,np.shape(current_cutout)[0])) #make the meshgrid for the wireframe
	text = gen_im_text(im_dat.iloc[i]) #generate text for this image
	lc_points = find_lc_point(im_dat.iloc[i], lc_dat)
	#Set the axes that have changed
	movers[0].set_clim(vmin, vmax)
	movers[1].axes.set_xlim(-0.5, len(points_hor)-0.5)
	movers[1].axes.set_ylim(min(points_hor)-1, max(points_hor)+1)
	movers[2].axes.set_xlim(min(points_ver)-1, max(points_ver)+1)
	movers[2].axes.set_ylim(-0.5, len(points_ver)-0.5)
	#Place the data in the plots
	movers[0].set_array(current_cutout)
	movers[1].set_data(np.arange(0,len(points_hor)), points_hor)
	movers[2].set_data(points_ver, np.flip(np.arange(0,len(points_ver))))
	ax_wire.collections.remove(movers[3])
	movers[3] = ax_wire.plot_wireframe(x_mesh, y_mesh, current_cutout, color='k', linewidth=0.5)
	movers[4].set_text(text)
	movers[5].set_offsets(lc_points)
	return movers

def find_lc_point(im_dat, lc_dat):
	if im_dat.im_type=='ref':
		return [[],[]]
	date = time.Time(im_dat.obsjd, format='jd')
	good_dates = lc_dat[np.abs(lc_dat.obsmjd-date.mjd)<1e-4] #dates are within 8.64 s of each other
	if len(good_dates[good_dates.obs_filter==im_dat.obs_filter])==0: #no matches
		return [[],[]]
	vals = []
	for i in range(len(good_dates[good_dates.obs_filter==im_dat.obs_filter])):
		xval = good_dates[good_dates.obs_filter==im_dat.obs_filter].obsmjd.iloc[i]
		yval = np.nanmin(good_dates[good_dates.obs_filter==im_dat.obs_filter][['mag','upper_limit']].iloc[i])
		vals.append([xval,yval])
	return vals

def gen_im_text(dat): #Generates text lines that are specific for the currently shown image
	text = 'Name {}\nImage type: {}\nFilter: {}\n\nDate: {}\nmjd: {}\n\nPeak value: {}\nMean value: {}\nStandard deviation: {}'
	if dat.im_type=='ref':
		return text.format(dat.obj_name, dat.im_type, dat.obs_filter, ' ', ' ', dat.peak_val, dat.mean_val, dat.standev)
	else:
		date = time.Time(dat.obsjd, format='jd')
		greg = date.isot
		date_greg = '{}-{}-{} {}:{}:{}'.format(greg[8:10], greg[5:7], greg[:4], greg[11:13], greg[14:16], greg[17:])
		return text.format(dat.obj_name, dat.im_type, dat.obs_filter, date_greg, date.mjd, dat.peak_val, dat.mean_val, dat.standev)

def get_lc(lc_loc, lc_dates):
	lc = pd.read_csv(lc_loc, header=0, usecols=['obsmjd', 'filter', 'upper_limit', 'mag', 'mag_err'])
	lc.rename(columns={'filter':'obs_filter'}, inplace=True) #To avoid confusion
	return lc[(lc.obsmjd>=min(lc_dates)) & (lc.obsmjd<=max(lc_dates))] #Return only the wanted dates

def gen_cutouts(cutout_dims, radec, refquery, difquery, obj_name):
	nr_refims = len(refquery.metatable)
	nr_frames = nr_refims+len(difquery.metatable) #amount of images to go through
	#Make an empty list containing cutouts
	cutouts = []
	#Also define an array storing relevant data for each image
	im_dat = pd.DataFrame(columns=['obj_name', 'im_type', 'obs_filter', 'obsjd', 'peak_val', 'mean_val', 'standev'], index=np.arange(0, nr_frames)) #Probably still needs to be expanded
	#Loop through the images to get required data
	data_locs = refquery.get_local_data()+difquery.get_local_data('scimrefdiffimg.fits.fz')
	for i in range(nr_frames):
		im = fits.open(data_locs[i])
		wcs = WCS(im[-1].header)
		pix_xy = wcs.all_world2pix(radec, 1)[0] #Pixel location of object
		cutouts.append(Cutout2D(im[-1].data, pix_xy, cutout_dims)) #Do the cutout
		#record data about image
		im_dat.obj_name[i] = obj_name
		if i<nr_refims:
			im_dat.im_type[i] = 'ref'
			im_dat.obsjd[i] = ' '
			if refquery.metatable.fid[i] == 1:
				im_dat.obs_filter[i] = 'ZTF_g'
			elif refquery.metatable.fid[i] == 2:
				im_dat.obs_filter[i] = 'ZTF_r'
			else:
				im_dat.obs_filter[i] = 'ZTF_i'
		else:
			im_dat.im_type[i] = 'diff'
			im_dat.obsjd[i] = im[-1].header['OBSJD']
			im_dat.obs_filter[i] = im[-1].header['FILTER']
		im_dat.peak_val[i] = np.nanmax(cutouts[i].data)
		im_dat.mean_val[i] = np.nanmean(cutouts[i].data)
		im_dat.standev[i] = np.nanstd(cutouts[i].data)
		#Close fits file before moving on to the next
		im.close()
	return cutouts, im_dat

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
