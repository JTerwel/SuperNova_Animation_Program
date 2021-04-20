#This can be used as a template showing how to get different images

from ztfquery import query
from astropy import time #For time conversions to jd

def main():
	obj_loc=[252.7534888, 25.8759437] #Holds ra, dec in degrees
	get_ref=False #Get reference image?
	get_dif=True #Get difference images?
	dist=0.01 #get all observations within dist from obj_loc (in degrees)
	#For the query note the following:
	#fid=1->zg, 2->zr, 3->zi filter
	#dates are in jd
	#obsjd BETWEEN start_jd and stop_jd, or obsjd>, obsjd=, etc.
	#seeing<2 for seeing below 2 arcsec
	t1=time.Time(58355, format='mjd').jd
	t2=time.Time(58365, format='mjd').jd
	t3=time.Time(58590, format='mjd').jd
	t4=time.Time(58750, format='mjd').jd
	q1='fid=2 and obsjd BETWEEN {} and {}'.format(t1,t2) #Get r im at peak for this obj
	#q2='fid=2 and obsjd BETWEEN {} and {}'.format(t3,t4) #Get r ims just before and during possible late-time detects
	#Try to get the 3 interesting observations in 1 go
	t=time.Time([58503, 58504, 58543, 58544, 58577, 58578], format='mjd').jd
	q1='fid=2 and (obsjd BETWEEN {} and {} or obsjd BETWEEN {} and {} or obsjd BETWEEN {} and {})'.format(t[0],t[1],t[2],t[3],t[4],t[5])
	print(q1)
	if get_ref:
		zquery=query.ZTFQuery()
		zquery.load_metadata(kind='ref', radec=obj_loc)
		zquery.download_data()
	if get_dif:
		zquery=query.ZTFQuery()
		#zquery2=query.ZTFQuery()
		zquery.load_metadata(radec=obj_loc, size=dist, sql_query=q1)
		#zquery2.load_metadata(radec=obj_loc, size=dist, sql_query=q2)
		zquery.download_data('scimrefdiffimg.fits.fz')
		#zquery2.download_data('scimrefdiffimg.fits.fz')
	return

if (__name__ == "__main__"):
	main()
