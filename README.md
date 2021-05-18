# animate_ZTFims
### Show ZTF objects evolving over time by showing the images one after another in an animation

This program (current version: animate_diffims_v5.py) makes an animation of the difference images of the chosen object, and shows additional data as requested by the user. The images are obtained using the ztfquery package. The animation meanwhile is made using matplotlib.animation.FuncAnimation.

## Basic usage
The settings are located in the main() function of the program and consist of 4 mandatory parameters (name, ra, dec, sql), and 10 optional parameters. Each parameter is described briefly:

* name: string, name of the object shown in the animation.
* ra: float, right ascension of the object in degrees.
* dec: float, declination of the object in degrees.
* sql: string, query stating which images should be selected by ztfquery.
* download_data: bool, default False. If False, all requested images are assumed to be stored locally and the download step will be skipped.
* lc_loc: string, default None. Location of the lightcurve csv file for the object.
* lc_dates: list of 2 floats, default [58119, 59215]. mjd of the lightcurve region to include in the animation.
* search_size: float, default 0.01. region around the target location in degrees in which to look for images.
* t_frame: int, default 500. Time in ms each frame is visible.
* cutout_dims: int / tuple, default 41. Pixel size of the image cutout shown in the animation. If float, width and height are the same. Use tuple to specify height and width separately. Note that this is done before any rotations and flips to put North up and East to the right, therefore the order in the tuple might change between images.
* min_visible: int / tuple, default 9. Size of the central region of the cutout where all pixels must have a proper value (e.g. no NaN or inf). If set to 0, this check is omitted. If float, width and height are the same. Use tuple to specify height and width separately.
* central_reg: int / tuple, default 5. Size of the central region of the cutout. It is used in finding the middle rows and columns to average for the midline plots (see below). If float horizontal and vertical sides are assumed to be the same, if tuple, the order is (height, width).
* subplot_selection: list, default None. List of init functions for the subplots to include in the animation. If None, all default subplots will be used, except for the lightcurve subplot is no lightcurve is given. Custom subplots can be included (see below).
* save_name: string, default None. Name under which the animation is to be saved. If None, the default save name will be used: 'name'.mp4 with 'name' being the name parameter.

Please do not forget to change call to the anim_target() object in line 61 according to the optional parameters used. To set the location where the animations are saved, please change the relevant part in line 67.

## Default subplots
There are 6 types of subplots that are used by default. They are listed below along with a short description:

* Image (init_im, update_im): The cutout image in grayscale, a red '+' marking the location of the associated lightcurve point (if found), and a green rectangle outlining the central region.
* Midlines (init_midlines, update_midlines): Line plots of the average value of the central columns (right of the image) and rows (above the image). The vertical lines of the green rectangle in the image mark the used rows, and the horizontal lines mark the used columns. (The image needs to be shown as well to see the green rectangle)
* Compass (init_compass, update_compass): As the axes of the pixels of the image are not necessarily aligned with North and East, this plot shows the actual cardinal directions.
* Wireframe (init_wireframe, update_wireframe): A 3D version of the image, where height is used instead of a grayscale to show the pixel values. The bottom corner in the wireframe corresponds to the bottom left corner of the image.
* Text (init_text, update_text): Some general information about the currently shown image: object name, image type (ref / diff), image filter, date (mjd and Gregorian), and image statistics (Peak value, mean value, and standar deviation).
* Lightcurve (init_lc, update_lc): The associated lightcurve plot of the object. A vertical dashed line shows the date of the current image, and if a lightcurve point is associated to the currently shown image (mjd differs by < 8.64 s & same filter), it is highlighted.

## Adding a custom subplot
Of course there is always the possibility that the subplot you want to be included in the animation is not (exactly) one of those described above. In that case, it is possible to make your own subplot by following the steps detailed below, and adding it to the animation using the subplot_selection parameter.

Each subplot consists of 2 functions: an init function to initiate the subplot and set its static parts, and an update function to change those parts that (might) need changing each frame. At the start of the animation the program goes through the list of init functions, after that the list of associated update function is cycled through for each frame. Because of this, the input and output for for these functions needs to be the same for all subplots.

First off, the init function, which looks like this:
```python
def init_subplot(target, fig, gs, subplots, movers, update_funcs, extra_params):
  #Initiate subplot
  subplot = fig.add_subplot(gs[11:, :50]) #Subplot is in the same spot as the image.
  
  #set static parts of the subplot
  subplot.set_xticks([]) #Remove the ticks on the x axis for the entire animation
  
  #Initiate the changing parts of the subplot without values.
  changing_line, = subplot.plot([], [], 'k') #A black line that will change every frame
  changing_line2, = subplot.plot([], [], 'r') #A red line that will change every frame
  
  #Add subplot and associate products to the lists and return them.
  subplots.append([subplot])
  movers.append([changing_line, changing_line2])
  update_funcs.append(update_subplot)
  extra_params.append([])
  return subplots, movers, update_funcs, extra_params
```
The arguments needed are the following:
* target: a anim_target object containing all settings and data related to the object.
* fig: The main figure
* gs: a GridSpec object used to easily specify subplot positions. It has dimensions (60, 140).
* subplots: a list of lists of subplots. (Usually 1 per init / update function, but midlines for instance controls 2 separate plots.)
* movers: a list of lists of moving parts in the subplots, in this case it is the line object.
* update_funcs: a list of update functions that has to be cycled through for each frame.
* extra_params: a list of lists of extra parameters that might be useful to have in the update function, if none, the added list is empty.

The function creates the subplot, sets the static and changing parts, and adds them to the list of subplots and movers. It also adds the associated update function to the list of update functions needed to call each step. Finally, if something extra is needed in the update function, it can be added in the extra_params list.

The location, type, and contents of the subplots is entirely customizable, as long as this general path is followed everything should work.

Each frame, the subplots are updated by calling the update functions from the list:
```python
for i in range(len(update_funcs)):
  subplots[i], movers[i], extra_params[i] = update_funcs[i](im_dat,
    target, subplots[i], movers[i], extra_params[i])
```

The associated update function has the following structure:
```python
def update_subplot(im_dat, target, subplot, mover, extra_params):
  #Update what needs to be updated each frame
  x = np.linspace(0,30,500)
  pow = np.nanmax(im_dat[0])
  mover[0].set_data(x, x**pow) #plot x to the power of the max value in the cutout image in black (ignoring NaNs)
  mover[1].set_data(x, x*pow) #plot x times the max value in the cutout image in red (ignoring NaNs)
  
  #Return the list items
  return subplot, mover, extra_params
```
The arguments that are given are the following:
* im_dat: A list of data on the image for the current frame, see below for more information.
* target: The anim_target object containing all settings and data related to the object.
* subplot: A list of the subplots in which the changes are being performed.
* mover: A list of the moving parts in the subplots.
* extra_params: A list of extra parameters needed in the update function, this list can be updated inside the update function before returning it.

The last three objects are returned at the end.

The reason subplot is passed into the update function as well is because while things images, scattered points, lines and text can be changed each step, this is not the case for things like Patches (e.g. a rectangle / circle) or arrows. They need to be removed and re-added to the subplot.

As said above, im_dat is a list of different things spit out by the generator function at the beginnig of making the next frame. The contents are as follows:
* im_dat[0]: A 2D numpy array containing the cutout image after the needed flips and rotations are performed for North to be up and East to be to the right.
* im_dat[1]: A string stating the type of image ('ref' or 'diff')
* im_dat[2]: The observation date as an astropy time object (For 'diff' images only, for 'ref' images imdat[2] = ' ')
* im_dat[3]: The filter in which the observation was taken ('ZTF_g', 'ZTF_r', 'ZTF_i').
* im_dat[4]: The number of 90 degrees clockwise rotations performed on the cutout image.
* im_dat[5]: A bool stating if the image was flipped vertically.
* im_dat[6]: The residual angle between the positive x axis and the East direction in radians.
* im_dat[7]: The coordinates of the cutout central pixel in the original image.

If the custom subplot is constructed according to the recipe above, and the init function is included in the subplot_selection list, it should be added to the animation.
