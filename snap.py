'''
Name: SNAP - SuperNova Animation Program
Author: Jacco Terwel
Date: 16 - 07 - 2021
SNAP retreives images from ZTF & puts them in chronological order in an
animation. It is made to be easily customizable to show exactly what is wanted.
'''

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


def main():
	'''
	Main function of the program.
	
	This program is made to run on its own, instead to be imported into another
	script. As such, This function is used to call the others in the correct
	order.
	'''
	#Required parameters
	name = 'ZTF18aaqeygf'
	ra = 200.3500294
	dec = 68.82638685
	t = time.Time([58260, 58400], format='mjd').jd
	sql = 'fid=2 and obsjd BETWEEN {} and {}\
		'.format(t[0], t[1])
	
	#Optional paramteres, comment out when the default value is used.
	download_data = False
	lc_loc = '/Users/terwelj/Documents/ZTFData/marshal/'\
		'SN_Ia/ZTF18aaqeygf_SNT_1e-08.csv'
	lc_dates = [58200, 58400]
	#search size = 
	#t_frame = 
	#cutout_dims = 
	#min_visible = 
	#central_reg = 
	subplot_selection = [
		init_im,
		init_midlines,
		init_compass,
		init_wireframe,
		init_text,
		init_lc]
	save_name = 'ZTF18aaqeygf_late_snap.mp4'
	
	#Make the target object, adjust depending to the optional parameters used.
	target = anim_target(name, ra, dec, sql, download_data=download_data,
		lc_loc=lc_loc, lc_dates=lc_dates,
		save_name=save_name)

	#Make and save the animation.
	anim = animate_observations(target)
	anim.save('/Users/terwelj/Movies/{}'.format(target.save_name))
	return


#*---------*
#| Classes |
#*---------*

class anim_target:
	'''
	Holds all relevant data required to make an animation.
	
	Attributes:
		name (string): Name of the target.
		ra (float): right ascencion of the target in degrees.
		dec (float): declination of the target in degrees.
		sql (string): SQL query stating which images are wanted.
		download_data (bool): Do the images have to be downloaded?
		lc_loc (string): Location of accompanying lightcurve.
		lc_dates (list): Visible time period of the lightcurve.
		size (float): Region around target for which all images are used.
		t_frame (int): Time each image is shown in the animation in miliseconds.
		cutout_dims (tuple): Size of shown image.
		min_visible (tuple): Minimum image size visible in all images.
		central_reg (tuple): nr of cols & rows considered for midline plots.
		init_list (list): List of initators of subplots used in the animation.
		save_name (string): Name under which the animation is saved.
		lc (DataFrame): Lightcurve made using lc_loc & lc_dates.
		im_locs (list): Location of all specified images (after downloading).
	'''

	def __init__(
			self, name, ra, dec, sql, download_data=False, lc_loc=None,
			lc_dates=[58119, 59215], search_size=0.01, t_frame=500,
			cutout_dims=41, min_visible=9, central_reg=5,
			subplot_selection=None, save_name=None):
		'''
		The constructor for this class.

		Parameters:
			name (string): Name of the target.
			ra (float): right ascencion of the target in degrees.
			dec (float): declination of the target in degrees.
			sql (string): SQL query stating which images are wanted.
			download_data (bool): Do the images have to be downloaded?
			lc_loc (string): Location of accompanying lightcurve.
			lc_dates (list): Visible time period of the lightcurve.
			size (float): Region around target for which all images are used.
			t_frame (int): Time each image is shown in the animation in ms
			cutout_dims (int / tuple): Size of shown image, internally tuple
				(height, width).
			min_visible (int / tuple): Minimum image size visible in all images,
				internally tuple (height, width).
			central_reg (int / tuple): nr of cols & rows considered for midline
				plots, internally tuple (height, width).
			subplot_selection (list): Custom subplot initiator list to be used.
			save_name (string): Name under which the animation is saved.
		'''
		self.name = name
		self.ra = ra
		self.dec = dec
		self.sql = sql
		self.download_data = download_data
		self.lc_loc = lc_loc				#If None, no lc wanted in animation.
		self.lc_dates = lc_dates
		self.size = search_size
		self.t_frame = t_frame
		if type(cutout_dims) is tuple:
			self.cutout_dims = cutout_dims
		else:
			self.cutout_dims = (cutout_dims, cutout_dims)
		if type(min_visible) is tuple:
			self.min_visible = min_visible
		else:
			self.min_visible = (min_visible, min_visible)
		if type(central_reg) is tuple:
			self.central_reg = central_reg
		else:
			self.central_reg = (central_reg, central_reg)
		#Name of animation defaults to name.mp4 with name the target name
		self.save_name = save_name if save_name is not None else name+'.mp4'
		if subplot_selection is not None:
			self.init_list = subplot_selection
		elif self.lc_loc is None:		#Default list of init functions
			self.init_list = [init_im, init_midlines, init_compass,
				init_wireframe, init_text]
		else:		#Default list of init functions including init_lc
			self.init_list = [init_im, init_midlines, init_compass,
				init_wireframe, init_text, init_lc]
		self.lc = self.get_lc() if self.lc_loc is not None else None
		self.im_locs = self.download_images()

	def get_lc(self):
		'''
		Get the lightcurve according to lc_loc & lc_dates

		Returns:
			lc_region: Section of the data within the required time period.
		'''
		if self.lc_loc is None:
			print('Warning: No lightcurve location given, returning None')
			return None
		lc = pd.read_csv(self.lc_loc, header=0,
			usecols=['obsmjd', 'filename', 'filter', 'upper_limit', 'mag', 'mag_err',
					'target_x', 'target_y'])
		lc.rename(columns={'filter':'obs_filter'}, inplace=True)#Avoid confusion
		lc_region = lc[(lc.obsmjd>=min(self.lc_dates))
					& (lc.obsmjd<=max(self.lc_dates))]
		return lc_region

	def download_images(self):
		'''
		Get all relevant images, download if needed & wanted,
		and list them in chronological order.

		Refquery is for the refence images, difquery for the difference images.

		Returns:
			finallist: list of locations of each image
		'''
		refquery = query.ZTFQuery()
		refquery.load_metadata(kind='ref', radec=[self.ra, self.dec])
		difquery = query.ZTFQuery()
		difquery.load_metadata(radec=[self.ra, self.dec], size=self.size,
			sql_query=self.sql)
		difquery.metatable.sort_values(by=['obsjd'], inplace = True)
		if self.download_data:	#Download images if it is requested.
			refquery.download_data()
			difquery.download_data('scimrefdiffimg.fits.fz')
		finallist = refquery.get_local_data() + difquery.get_local_data(
			'scimrefdiffimg.fits.fz')
		return finallist


#*----------------------------------*
#| Functions to make the animation. |
#*----------------------------------*

def animate_observations(target):
	'''
	Put the images of the target object into a slideshow

	Parameters:
		target (anim_target): Object to animate.

	Returns:
		anim: The generated animation of the target.
	'''
	#Initiate the lists that keep track of everything about each subplot.
	#subplots[i] is / are updated in update_funcs[i].
	#Movers store the specific parts of the subplots to be updated.
	#update_funcs is a list of update functions to call each step.
	#extra_params[i] lists the extra parameters that update_funcs[i] needs.
	subplots = []
	movers = []
	update_funcs = []
	extra_params = []

	#Initiate figure and subplots
	gs = gridspec.GridSpec(60,140)
	fig = plt.figure(figsize=(14,6))
	for i in target.init_list:
		subplots, movers, update_funcs, extra_params = i(target, fig, gs,
			subplots, movers, update_funcs, extra_params)

	#Make and return the animation.
	anim = FuncAnimation(fig, animate_step, frames=get_next_im(target),
		interval=target.t_frame,
		fargs=(target, subplots, movers, update_funcs, extra_params),
		save_count=len(target.im_locs), blit=True, repeat=False)
	return anim

def animate_step(im_dat, target, subplots, movers, update_funcs, extra_params):
	'''
	Make the next image in the animation using the update functions provided.

	Parameters:
		im_dat (list): The next yield from the generator function.
		target (anim_target): Object to animate.
		subplots (list): List of subplots in the animation.
		movers (list): List of moving parts in each subplot.
		update_funcs (list): List of update functions for each subplot.
		extra_params (list): List of extra parameters needed in the update_func.

	Returns:
	flat_movers (list): list of artists that are updated.
	'''
	for i in range(len(update_funcs)):
		subplots[i], movers[i], extra_params[i] = update_funcs[i](im_dat,
			target, subplots[i], movers[i], extra_params[i])

	#Make a flat version of movers to return
	flat_movers = [item for sublist in movers for item in sublist]
	return flat_movers

def get_next_im(target):
	'''
	Generator function for making the animation

	The next image in the list is opened, the target location is found and cut
	out, and the cutout is tested. If the tests are passed successfully, the
	cutout will be used in the next frame of the animation. Else a warning is
	issued before continueing to the next image.

	1st test: Is the minimum visible region fully observed (e.g. no NaN)?

	2nd test: Are there positive values? (prevents weird negative frames found
	when using previous versions of the code)

	Parameters:
		target (anim_target): Object to animate.

	Yields:
		list: A list containing the rotated cutout, rotation details, and
		some metadata.
	'''
	for i in target.im_locs:
		#Open next image and check if it meets the animation requirements.
		im = fits.open(i)[-1]
		wcs = WCS(im.header)
		pix_xy = wcs.all_world2pix([(target.ra, target.dec)], 0)[0]

		#1st test
		if ((target.min_visible!=0) & (target.min_visible!=(0,0))):
			try: #Should only succeed if all pixels in this cutout have values
				dummy = Cutout2D(im.data, pix_xy, target.min_visible,
					mode='strict')
			except:		#Can't fully fill min visible region
				print('Image at obsjd {} rejected due to NaN close to im \
						center'.format(im.header['OBSJD']))
				continue #Start next i in loop

		cutout = Cutout2D(im.data, pix_xy, target.cutout_dims, mode='partial',
			wcs=wcs)

		#2nd test
		if np.max(cutout.data) < 0:
			print('Image at obsjd {} rejected due to all values being negative\
				'.format(im.header['OBSJD']))
			continue
		
		#Get needed info from the header (Give entire header in the future?)
		if 'ZTFData/ref' in i:
			imtype = 'ref'
		else:
			imtype = 'dif'
		
		#Correct image orientation / being mirrored if needed.
		immat, rots, flip, resid = find_orientation(cutout)
		im_center = [round(pix_xy[0]), round(pix_xy[1])]

		#yield a list of items required by the subplots
		yield [immat, imtype, im.header, rots, flip, resid, im_center]


#*-----------------------------------------------------------*
#| Initiators and update functions for all default subplots. |
#*-----------------------------------------------------------*

def init_im(target, fig, gs, subplots, movers, update_funcs, extra_params):
	'''
	Iniate the image cutout & colorbar.

	Parameters:
		target (anim_target): Object about to be animated.
		fig (figure): Figure on which the subplot is added.
		gs (GridSpec): Gridspec object associated with fig.
		subplots (list): List of subplots already in fig.
		movers (list): List of moving objects in subplots already in fig.
		update_funcs (list): List of update functions associated with subplots.
		extra_params (list): List of extra parameters in the update functions.

	Returns:
		subplots (list): Updated list of subplots.
		movers (list): Updated list of movers.
		update_funcs (list): Updated list of update_funcs.
		extra_params (list): Updated list of extra_params.
	'''

	#Initiate the subplot and colorbar.
	ax_im = fig.add_subplot(gs[11:, :50])
	ax_cb = fig.add_axes([0.1, 0.11, 0.015, 0.625])

	#Set ticks and labels
	ax_im.set_xticklabels([])
	ax_im.set_yticklabels([])
	ax_im.set_xticks([int(target.cutout_dims[1]/2)])
	ax_im.set_yticks([int(target.cutout_dims[0]/2)])
	central_bot = int(target.cutout_dims[0]/2) - int(target.central_reg[0]/2)
	central_left = int(target.cutout_dims[1]/2) - int(target.central_reg[1]/2)
	empty_im = np.zeros(target.cutout_dims)
	
	#Set the central region rectangle
	rect = Rectangle((central_left-0.5, central_bot-0.5), target.central_reg[1],
		target.central_reg[0], edgecolor='g', facecolor='none', lw=1.5, zorder=1)
	ax_im.add_patch(rect)

	#Set changing parts: the image (im), forced photometry location (lc_loc),
	#and colorbar (cb).
	im = ax_im.imshow(empty_im, cmap='gray', origin='lower', zorder=0)
	lc_loc = ax_im.scatter([], [], color='r', marker='+', zorder=2)
	cb = fig.colorbar(im, cax = ax_cb, ticklocation='left')

	#Add to the lists and return.
	subplots.append([ax_im, ax_cb])
	movers.append([im, lc_loc])
	update_funcs.append(update_im)
	extra_params.append([])
	return subplots, movers, update_funcs, extra_params

def update_im(im_dat, target, subplot, mover, extra_args):
	'''
	Update the image and colorbar.

	Parameters:
		im_dat ([immat, imtype, header, rots, flip, resid, im_center]):
			Required image data.
		target (anim_target): Object to animate.
		subplot ([ax_im, ax_cb]): Subplots to be updated.
		mover ([im, lc_loc]): Parts of the subplot that need to be updated.
		extra_args ([]): Extra arguments that are needed.

	Returns:
		subplot (list): Updated version of the subplot given.
		mover (list): Updated version of the mover ginve.
		extra_args (list): Updated version of the extra args given.
	'''
	#Place new image.
	mover[0].set_array(im_dat[0])

	#Update colorbar to min and max value in the image.
	mover[0].set_clim(np.nanmin(im_dat[0]), np.nanmax(im_dat[0]))

	#Mark forced photometry spot.
	lc_locs = find_fp_loc(target.lc, im_dat, target.cutout_dims)
	mover[1].set_offsets(lc_locs)

	#Return updated subplot.
	return subplot, mover, extra_args

def init_midlines(target, fig, gs, subplots, movers, update_funcs,
		extra_params):
	'''
	Initiate the lines through central rows and columns of the image.

	Parameters:
		target (anim_target): Object about to be animated.
		fig (figure): Figure on which the subplot is added.
		gs (GridSpec): Gridspec object associated with fig.
		subplots (list): List of subplots already in fig.
		movers (list): List of moving objects in subplots already in fig.
		update_funcs (list): List of update functions associated with subplots.
		extra_params (list): List of extra parameters in the update functions.

	Returns:
		subplots (list): Updated list of subplots.
		movers (list): Updated list of movers.
		update_funcs (list): Updated list of update_funcs.
		extra_params (list): Updated list of extra_params.
	'''

	#Initiate
	ax_line_hor = fig.add_subplot(gs[:10, :50])
	ax_line_ver = fig.add_subplot(gs[11:, 50:60])

	#Set static parts including dashed lines at 0
	ax_line_hor.set_xticks([])
	ax_line_ver.set_yticks([])
	ax_line_hor.axhline(0, ls='--', color='k', alpha=0.7, lw=0.5)
	ax_line_ver.axvline(0, ls='--', color='k', alpha=0.7, lw=0.5)
	hor_offset = (max(target.cutout_dims)-target.cutout_dims[1])/2
	ver_offset = (max(target.cutout_dims)-target.cutout_dims[0])/2
	ax_line_hor.axes.set_xlim(-0.5-hor_offset,
		max(target.cutout_dims)-0.5-hor_offset)
	ax_line_ver.axes.set_ylim(-0.5-ver_offset,
		max(target.cutout_dims)-0.5-ver_offset)
	#Set changing parts
	line_hor, = ax_line_hor.plot([], [], 'k')
	line_ver, = ax_line_ver.plot([], [], 'k')

	#Calculate the edges of the central regions & set some scale limits.
	bot = int(target.cutout_dims[0]/2) - int(target.central_reg[0]/2)
	left = int(target.cutout_dims[1]/2) - int(target.central_reg[1]/2)
	top = bot + target.central_reg[0]
	right = left + target.central_reg[1]

	#Add to the lists and return.
	subplots.append([ax_line_hor, ax_line_ver])
	movers.append([line_hor, line_ver])
	update_funcs.append(update_midlines)
	extra_params.append([top, bot, left, right])
	return subplots, movers, update_funcs, extra_params

def update_midlines(im_dat, target, subplot, mover, extra_args):
	'''
	Update the midline plots.

	Parameters:
		im_dat ([immat, imtype, header, rots, flip, resid, im_center]):
			required image data.
		target (anim_target): Object to animate.
		subplot ([ax_line_hor, ax_line_ver]): Subplots to be updated.
		mover ([line_hor, line_ver]): Parts of the subplot to be updated.
		extra_args ([top, bot, left, right]): Extra arguments that are needed.

	Returns:
		subplot (list): Updated version of the subplot given.
		mover (list): Updated version of the mover ginve.
		extra_args (list): Updated version of the extra args given.
	'''
	#Get the averaged middle rows and columns.
	points_hor = np.average(im_dat[0][extra_args[1]:extra_args[0],:], axis=0)
	points_ver = np.average(im_dat[0][:, extra_args[2]:extra_args[3]], axis=1)

	#get the min & max values for the scaling.
	horver_minmax = [np.nanmin(points_hor), np.nanmin(points_ver),
		np.nanmax(points_hor), np.nanmax(points_ver)]
	for i in range(len(horver_minmax)): 		#Change unexpected output to 0
		if np.isnan(horver_minmax[i]) or np.isinf(horver_minmax[i]):
			horver_minmax[i]=0
	
	#Set new scaling
	mover[0].axes.set_ylim(horver_minmax[0]-1, horver_minmax[2]+1)
	mover[1].axes.set_xlim(horver_minmax[1]-1, horver_minmax[3]+1)

	#Place the new lines
	mover[0].set_data(np.arange(0,len(points_hor)), points_hor)
	mover[1].set_data(points_ver, np.arange(0,len(points_ver)))

	#Return the updated subplot.
	return subplot, mover, extra_args


def init_compass(target, fig, gs, subplots, movers, update_funcs, extra_params):
	'''
	Initiate the plot with arrows pointing towards North and East.

	Parameters:
		target (anim_target): Object about to be animated.
		fig (figure): Figure on which the subplot is added.
		gs (GridSpec): Gridspec object associated with fig.
		subplots (list): List of subplots already in fig.
		movers (list): List of moving objects in subplots already in fig.
		update_funcs (list): List of update functions associated with subplots.
		extra_params (list): List of extra parameters in the update functions.

	Returns:
		subplots (list): Updated list of subplots.
		movers (list): Updated list of movers.
		update_funcs (list): Updated list of update_funcs.
		extra_params (list): Updated list of extra_params.
	'''

	#Initiate
	ax_ne = fig.add_subplot(gs[:10, 50:60])

	#Set static parts
	ax_ne.set_xlim(0, 1)
	ax_ne.set_ylim(0, 1)
	ax_ne.spines['top'].set_visible(False)
	ax_ne.spines['bottom'].set_visible(False)
	ax_ne.spines['left'].set_visible(False)
	ax_ne.spines['right'].set_visible(False)

	#No moving parts, everything is erased and replaced completely.
	#Add to the lists and return.
	subplots.append([ax_ne])
	movers.append([])
	update_funcs.append(update_compass)
	extra_params.append([])
	return subplots, movers, update_funcs, extra_params

def update_compass(im_dat, target, subplot, mover, extra_args):
	'''
	Reset the compass.

	Parameters:
		im_dat ([immat, imtype, header, rots, flip, resid, im_center]):
			Required image data.
		target (anim_target): Object to animate.
		subplot ([ax_ne]): Subplots to be updated.
		mover ([]): Parts of the subplot that need to be updated.
		extra_args ([]): Extra arguments that are needed.

	Returns:
		subplot (list): Updated version of the subplot given.
		mover (list): Updated version of the mover ginve.
		extra_args (list): Updated version of the extra args given.
	'''
	#Remove the old compass.
	subplot[0].clear()

	#Add the new arrows.
	subplot[0].arrow(0.5, 0.5, 0.3*np.cos(im_dat[5]), 0.3*np.sin(im_dat[5]),
		 color='k', length_includes_head=True, head_width=0.1, head_length=0.1)
	subplot[0].arrow(0.5, 0.5, -0.3*np.sin(im_dat[5]), 0.3*np.cos(im_dat[5]),
		color='k', length_includes_head=True, head_width=0.1,
		head_length=0.1)

	#Add the new markers.
	subplot[0].text(0.5+0.4*np.cos(im_dat[5]), 0.5+0.4*np.sin(im_dat[5]), 'E',
		color='k', va='center', ha='center')
	subplot[0].text(0.5-0.4*np.sin(im_dat[5]), 0.5+0.4*np.cos(im_dat[5]), 'N',
		color='k', va='center', ha='center')

	#Remove the ticks that have sprung up
	subplot[0].set_xticks([])
	subplot[0].set_yticks([])

	#Return updated subplot
	return subplot, mover, extra_args


def init_wireframe(target, fig, gs, subplots, movers, update_funcs,
		extra_params):
	'''
	Initiate the wireframe.

	Parameters:
		target (anim_target): Object about to be animated.
		fig (figure): Figure on which the subplot is added.
		gs (GridSpec): Gridspec object associated with fig.
		subplots (list): List of subplots already in fig.
		movers (list): List of moving objects in subplots already in fig.
		update_funcs (list): List of update functions associated with subplots.
		extra_params (list): List of extra parameters in the update functions.

	Returns:
		subplots (list): Updated list of subplots.
		movers (list): Updated list of movers.
		update_funcs (list): Updated list of update_funcs.
		extra_params (list): Updated list of extra_params.
	'''

	#Initiate
	ax_wire = fig.add_subplot(gs[:28, 65:100], projection='3d')

	#Set static parts
	ax_wire.set_xticks([])
	ax_wire.set_yticks([])
	#Elevation, azimut (45,225) points bottom left corner of image down
	#Might want to change this to a free parameter
	ax_wire.view_init(45, 225) 
	empty_im = np.zeros(target.cutout_dims)
	x_mesh, y_mesh = np.meshgrid(np.arange(0, target.cutout_dims[1]),
		np.arange(0, target.cutout_dims[0]))

	#Set changing parts
	wire = ax_wire.plot_wireframe(x_mesh, y_mesh, empty_im, color='k', linewidth=0.5)

	#Add to the lists and return.
	subplots.append([ax_wire])
	movers.append([wire])
	update_funcs.append(update_wireframe)
	extra_params.append([])
	return subplots, movers, update_funcs, extra_params

def update_wireframe(im_dat, target, subplot, mover, extra_args):
	'''
	Update the wireframe.

	Parameters:
		im_dat ([immat, imtype, header, rots, flip, resid, im_center]):
			Required image data.
		target (anim_target): Object to animate.
		subplot ([ax_wire]): Subplots to be updated.
		mover ([wire]): Parts of the subplot that need to be updated.
		extra_args ([]): Extra arguments that are needed.

	Returns:
		subplot (list): Updated version of the subplot given.
		mover (list): Updated version of the mover ginve.
		extra_args (list): Updated version of the extra args given.
	'''
	#Make the wireframe meshgrid.
	x_mesh, y_mesh = np.meshgrid(np.arange(0,np.shape(im_dat[0])[1]),
		np.arange(0,np.shape(im_dat[0])[0]))

	#Replace the wireframe
	subplot[0].collections.remove(mover[0])
	mover[0] = subplot[0].plot_wireframe(x_mesh, y_mesh, im_dat[0], color='k',
		linewidth=0.5)

	#Return updated subplot
	return subplot, mover, extra_args


def init_text(target, fig, gs, subplots, movers, update_funcs, extra_params):
	'''
	Initiate the text plot.

	Parameters:
		target (anim_target): Object about to be animated.
		fig (figure): Figure on which the subplot is added.
		gs (GridSpec): Gridspec object associated with fig.
		subplots (list): List of subplots already in fig.
		movers (list): List of moving objects in subplots already in fig.
		update_funcs (list): List of update functions associated with subplots.
		extra_params (list): List of extra parameters in the update functions.

	Returns:
		subplots (list): Updated list of subplots.
		movers (list): Updated list of movers.
		update_funcs (list): Updated list of update_funcs.
		extra_params (list): Updated list of extra_params.
	'''

	#Initiate
	ax_text = fig.add_subplot(gs[:28, 101:])

	#Set static parts
	ax_text.set_xticks([])
	ax_text.set_yticks([])
	ax_text.set_xlim(0, 1)
	ax_text.set_ylim(0, 1)
	ax_text.spines['top'].set_visible(False)
	ax_text.spines['bottom'].set_visible(False)
	ax_text.spines['left'].set_visible(False)
	ax_text.spines['right'].set_visible(False)

	#Set changing parts
	text = ax_text.text(0, 1, ' ', va='top')

	#Add to the lists and return.
	subplots.append([ax_text])
	movers.append([text])
	update_funcs.append(update_text)
	extra_params.append([target.name])
	return subplots, movers, update_funcs, extra_params

def update_text(im_dat, target, subplot, mover, extra_args):
	'''
	Update the text plot.

	Parameters:
		im_dat ([immat, imtype, header, rots, flip, resid, im_center]):
			Required image data.
		target (anim_target): Object to animate.
		subplot ([ax_text]): Subplots to be updated.
		mover ([text]): Parts of the subplot that need to be updated.
		extra_args ([name]): Extra arguments that are needed.

	Returns:
		subplot (list): Updated version of the subplot given.
		mover (list): Updated version of the mover ginve.
		extra_args (list): Updated version of the extra args given.
	'''
	#Generate the text
	text = gen_text(extra_args[0], im_dat)

	#Place text
	mover[0].set_text(text)

	#return updated subplot
	return subplot, mover, extra_args


def init_lc(target, fig, gs, subplots, movers, update_funcs, extra_params):
	'''
	Initiate the lightcurve plot.

	Parameters:
		target (anim_target): Object about to be animated.
		fig (figure): Figure on which the subplot is added.
		gs (GridSpec): Gridspec object associated with fig.
		subplots (list): List of subplots already in fig.
		movers (list): List of moving objects in subplots already in fig.
		update_funcs (list): List of update functions associated with subplots.
		extra_params (list): List of extra parameters in the update functions.

	Returns:
		subplots (list): Updated list of subplots.
		movers (list): Updated list of movers.
		update_funcs (list): Updated list of update_funcs.
		extra_params (list): Updated list of extra_params.
	'''

	#Initiate
	ax_lc = fig.add_subplot(gs[31:, 67:])

	#Set static parts
	for filt in ['g', 'r', 'i']:
		selection = target.lc[target.lc.obs_filter.str.contains(filt)]
		if filt=='i': #Now filt can be used to set the color of the points
			filt='y'
		ax_lc.scatter(selection[selection.mag==99].obsmjd,
			selection[selection.mag==99].upper_limit, color=filt, alpha=0.3,
			marker='v', s=4, zorder=1) #upper limits
		ax_lc.errorbar(selection[selection.mag!=99].obsmjd,
			selection[selection.mag!=99].mag,
			yerr=selection[selection.mag!=99].mag_err, fmt='o', color=filt,
			ms=4, zorder=2) #detections
	ax_lc.set_xlim(target.lc_dates[0]-1, target.lc_dates[1]+1)
	ax_lc.set_ylim(23, 15)
	ax_lc.set_xlabel('mjd')
	ax_lc.set_ylabel('mag')

	#Set changing parts
	lc = ax_lc.scatter([], [], edgecolors='k', facecolors='none', linewidths=3,
		s=64, zorder=3) #Hollow circle to place over associated point if found
	lc_line = ax_lc.axvline(0, ls='--', color='k')

	#Add to the lists and return.
	subplots.append([ax_lc])
	movers.append([lc, lc_line])
	update_funcs.append(update_lc)
	extra_params.append([])
	return subplots, movers, update_funcs, extra_params

def update_lc(im_dat, target, subplot, mover, extra_args):
	'''
	Update the highlighted data point(s) in the light curve.

	Parameters:
		im_dat ([immat, imtype, header, rots, flip, resid, im_center]):
			Required image data.
		target (anim_target): Object to animate.
		subplot ([ax_lc]): Subplots to be updated.
		mover ([lc, lc_line]): Parts of the subplot that need to be updated.
		extra_args ([]): Extra arguments that are needed.

	Returns:
		subplot (list): Updated version of the subplot given.
		mover (list): Updated version of the mover ginve.
		extra_args (list): Updated version of the extra args given.
	'''
	#Find points to highlight and update the line positioning.
	if im_dat[1] =='ref':
		lc_points = [[],[]]
		mover[1].set_xdata([0, 0])
	else:
		lc_points, mjd = find_fp_points(target.lc, im_dat[2])
		mover[1].set_xdata([mjd, mjd])

	#Update highlighted points.
	mover[0].set_offsets(lc_points)

	#Return updated subplot
	return subplot, mover, extra_args


#*--------------------------------------------*
#| Helper functions called in functions above |
#*--------------------------------------------*

def find_orientation(im):
	'''
	Find the orientation of im, rotate and/or flip it to make sure North
	is up and East is to the right.

	Parameters:
		im (astropy.fits): image that needs to be rotated

	Returns:
		data (2d array): rotated image matrix.
		rots (int): nr of clockwise rotations performed.
		flip (bool): flip needed True/False.
		resid (float): residual angle.
	'''

	#Find current direction of North and East
	wcs = im.wcs
	pix00 = wcs.pixel_to_world(0,0)
	to_east = wcs.world_to_pixel(SkyCoord(ra=pix00.ra+1*u.arcmin,
		dec=pix00.dec, frame='icrs'))
	to_north = wcs.world_to_pixel(SkyCoord(ra=pix00.ra,
		dec=pix00.dec+1*u.arcmin, frame='icrs'))

	#Flip image if needed and find resulting angle of East.
	if (((to_east[0]>0) & (to_north[1]<0))
			| ((to_east[0]<0) & (to_north[1]>0))):
		data = np.flipud(im.data)
		theta = np.arctan2(-to_east[1], to_east[0])
		flip = True

	else: #image is not mirrored
		data = im.data
		theta = np.arctan2(to_east[1], to_east[0])
		flip = False

	#convert to 0-2pi rad
	if theta<0:
		theta = theta+2*np.pi

	#Rotate image as required amount and find residual theta in rad
	if theta%(np.pi/2)>np.pi/4: #Make sure east is not to far up
		theta += np.pi/2
		resid = theta%(np.pi/2)-np.pi/2
	else:
		resid = theta%(np.pi/2)
	rots = int(theta/(np.pi/2)) #nr of clockwise rotations
	data = np.rot90(data, rots, (1,0)) #rotated data
	return data, rots, flip, resid

def find_fp_loc(lc, im_dat, cutout_dims):
	'''
	Find the location of the associated forced photometry point(s).

	Parameters:
		lc (DataFrame): lightcurve with forced photometry points.
		im_dat ([immat, imtype, header, rots, flip, resid, im_center]):
			Required image data.
		cutout_dims (tuple or int): Size of the cutout image.

	Returns:
		lc_locs (list): list of rotation and flip corrected lc coordinates
			found to match the image.
	'''
	#Check if a lc is given and there is a date from the reference image.
	if lc is None or im_dat[1] == 'ref':
		return [[None, None]]
	else:
		date = time.Time(im_dat[2]['OBSJD'], format='jd')

	#Find all lc points derived from this image
	good_points = lc[lc.filename == im_dat[2]['ORIGNAME']]
	if len(good_points)==0:	#no matches
		return [[None, None]]
	xvals = good_points.target_x.to_numpy().transpose()	#Put in simple arrays
	yvals = good_points.target_y.to_numpy().transpose()

	#Calculate x & y values with respect to cutout center while taking the flip
	#into account if needed.
	diff_x = xvals-im_dat[6][0]
	if im_dat[4]:
		diff_y = im_dat[6][1] - yvals
	else:
		diff_y = yvals - im_dat[6][1]

	#Add rotations and move origin to the actual origin.
	if im_dat[3]%4==1: #1 time clockwise
		locs = [diff_y+cutout_dims[1]/2, -diff_x+cutout_dims[0]/2]
	elif im_dat[3]%4==2: #im was upside down
		locs = [-diff_x+cutout_dims[1]/2, -diff_y+cutout_dims[0]/2]
	elif im_dat[3]%4==3: #1 time anticlockwise
		locs = [-diff_y+cutout_dims[1]/2, diff_x+cutout_dims[0]/2]
	else: #rots%4==0 means do nothing
		locs = [diff_x+cutout_dims[1]/2, diff_y+cutout_dims[0]/2]
	
	#Put in appropiate format to return
	lc_locs = []
	for i in range(len(diff_x)):
		lc_locs.append([locs[0][i], locs[1][i]])
	return lc_locs

def find_fp_points(lc, header):
	'''
	Find the forced photometry point(s) belonging to the image.

	Parameters:
		lc (DataFrame): Lightcurve with forced photometry points.
		header: image header

	Returns:
		lc_points (list): list of lightcurve points belonging to the image.
		mjd (Time): mjd of the image
	'''
	#Check if a lc is given and there is a date from the reference image.
	if lc is None:
		return [[],[]]
	date = time.Time(header['OBSJD'], format='jd')
	im_filter = header['FILTER']

	#Find all lc points originating from this image
	good_dates = lc[lc.filename==header['ORIGNAME']]
	if len(good_dates)==0:	#no matches
		return [[],[]]

	#List all good data points, take mag when possible, else take upper_limit.
	vals = []
	for i in range(len(good_dates[good_dates.obs_filter==im_filter])):
		xval = good_dates[good_dates.obs_filter==im_filter].obsmjd.iloc[i]
		yval = np.nanmin(good_dates[good_dates.obs_filter==im_filter][['mag',
			'upper_limit']].iloc[i])
		vals.append([xval,yval])
	
	#Return the list of lightcurve points belonging to the image.
	return vals, date.mjd

def gen_text(name, im_dat):
	'''
	Generate the text for the next frame.

	Parameters:
		name (string): Object name.
		im_dat ([immat, imtype, header, rots, flip, resid, im_center]):
			Required image data.
	Returns:
		text (string): The text to be placed in the sublplot
	'''
	#Set text preset
	text = 'Name {}\nImage type: {}\nFilter: {}\nExposure time (s): {}\n\nDate:'\
		' {}\nmjd: {}\n\nPeak value: {}\nMean value: {}\nStandard deviation: {}'

	#Calculate some image statistics
	peak = np.nanmax(im_dat[0])
	mean = np.nanmean(im_dat[0])
	stddev = np.nanstd(im_dat[0])

	#Fill in & return text preset based on image type (ref / diff).
	if im_dat[1] == 'ref':
		date = ' '
		exptime = '{} ({} coadded frames)'.format(im_dat[2]['TOTEXPT'],
			im_dat[2]['NFRAMES'])
		if im_dat[2]['DBFID'] == 1:
			obsfilter = 'ZTF_g'
		elif im_dat[2]['DBFID'] == 2:
			obsfilter = 'ZTF_r'
		elif im_dat[2]['DBFID'] == 3:
			obsfilter = 'ZTF_i'
		else:
			obsfilter = 'UNKNOWN'
		return text.format(name, im_dat[1], obsfilter, exptime, ' ', ' ', peak,
			mean, stddev)
	else:
		date = time.Time(im_dat[2]['OBSJD'], format='jd')
		obsfilter = im_dat[2]['FILTER']
		exptime = im_dat[2]['EXPOSURE']
		#Calculate the date in Gregorian format (easier to read than mjd).
		greg = date.isot
		date_greg = '{}-{}-{} {}:{}:{}'.format(greg[8:10], greg[5:7], greg[:4],
			greg[11:13], greg[14:16], greg[17:])
		return text.format(name, im_dat[1], obsfilter, exptime, date_greg,
			date.mjd, peak, mean, stddev)


#*---------------------------------------------*
#| Start the program through its main function |
#*---------------------------------------------*

if (__name__ == '__main__'):
	main()
