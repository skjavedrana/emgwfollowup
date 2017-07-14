"""
Apply the optimal algorithm on a GW-localization over setting and observation duration
constraints. It also check whether the GW patch has rised or not. The algorithm is 
optimal for multi-epoch observations.

print summary:
    # observation start time
    # ra & dec to start the pointing
    # the amount of area can be covered 
    # the covered probability in covered area

The output is saved as output.csv file.
"""

import argparse
from astropy.time import Time

parser = argparse.ArgumentParser()
parser.add_argument(
    'input', metavar='input.csv', type=argparse.FileType('rb'),
    default='-', help='Name of input file [default: stdin]')
parser.add_argument(
    '-o', '--output', metavar='output.csv', type=argparse.FileType('wb'),
    default='-', help='Name of output file [default: stdout]')
parser.add_argument(
    '--method', type=str, choices=['greedy', 'optimal', 'sear'],
    default='optimal', help='Scheduling method')
parser.add_argument(
    '--observetory-lat', type=float,
    default=33.355833, help='Latitude of the observetory [default: %(default)s]')
parser.add_argument(
    '--observetory-lon', type=float,
    default=-116.863889, help='Longitude of the Observetory [default: %(default)s]')
parser.add_argument(
    '--fov-width', type=float,
    default=3.5, help='Width of the FOV [default: %(default)s]')
parser.add_argument(
    '--fov-hight', type=float,
    default=2.31, help='Hight of the FOV [default: %(default)s]')
parser.add_argument(
    '--altitude', type=float,
    default=1712, help='Altitude of the Telescope [default: %(default)s]')
parser.add_argument(
    '--horizon', type=float,
    default=30., help='Visibility horizon for the Telescope [default: %(default)s]')
parser.add_argument(
    '--twilight', type=float,
    default=-12., help='Observation starts, when Sun goes belowtwilight [default: %(default)s]')
parser.add_argument(
    '--time', type=str,
    default=Time.now(), help='Observation start time [default: %(default)s]')
parser.add_argument(
    '--texpo', type=float,
    default=106, help='Time of exposure [default: %(default)s]')
parser.add_argument(
    '--numb-telescopes', type=int,
    default=1, help='Number of telescopes in the Observetory [default: %(default)s]')
parser.add_argument(
    '--epoch', type=int,
    default=1, help='Number of EPOCH must be taken [default: %(default)s]')
parser.add_argument(
    '--timegap', type=float,
    default=0., help='The time-gap between two epoch [default: %(default)s]')
parser.add_argument(
    '--cumplot', type=int,
    default=False, help='Make the cumulative probability plots [default: %(default)s]')
parser.add_argument(
    '--skyplot', type=int,
    default=False, help='Make the sky-plots of the tiles [default: %(default)s]')
opts = parser.parse_args()



import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from ephem import *


class Observetory:
    """The detail of the Ovservetory"""

    fov_w = opts.fov_width
    fov_h = opts.fov_hight
    lat = opts.observetory_lat
    lon = opts.observetory_lon
    altitude = opts.altitude
    horizon = opts.horizon
    twilight = opts.twilight

observetory = Observetory()
tileinfo = ascii.read(opts.input)


start_time = Date(str(opts.time))
texpo = opts.texpo * second
if opts.epoch: texpo = opts.texpo * opts.epoch * second

"""Get the observer
We are using ephem for defining the observer."""
observer = Observer()
observer.lat = str(observetory.lat)
observer.lon = str(observetory.lon)
observer.date = start_time
observer.horizon = str(observetory.twilight)
observer.elevation = observetory.altitude

# Sunset and Sunrise at observetory
observetory.type = 'equatorial'
sun = Sun()
observetory.sunset = observer.next_setting(sun)
observetory.sunrise = observer.next_rising(sun)
if observetory.sunset > observetory.sunrise:
    observetory.sunset = observer.previous_setting(sun)
observer.horizon = str(observetory.horizon)
sunrise_time = observetory.sunrise


class time:
    """Start ans End time"""

    start = start_time
    end = observetory.sunrise
    #if opts.epoch: end = observetory.sunrise - (opts.epoch-1) * opts.timegap * second
time = time()
    

def radec2altaz(tileinfo, observetory, time):
    """Convert the locations of the tiles in local coordinate."""

    time = Time(julian_date(time), format='jd')

    observetory.location = EarthLocation(lat=observetory.lat * u.deg, lon=observetory.lon * u.deg)
    cc=SkyCoord(tileinfo['ra']*u.deg, tileinfo['dec']*u.deg)
    aa=cc.transform_to(AltAz(location=observetory.location, obstime=time))
    return aa.alt.degree, aa.az.degree


def abovehorizon(tileinfo, time):
    """Select the tiles above and below Horizon"""

    alt, az = radec2altaz(tileinfo, observetory, time)
    onsky = alt > observetory.horizon # above Horizon
    if observetory.type == 'equatorial':
        tileinfo_up = tileinfo[onsky]
	tileinfo_down = tileinfo[True^onsky] 		
    return tileinfo_up, tileinfo_down


def setrise_tileinfo(tileinfo, time):
    """Get the rising and setting times of all of the tiles"""

    zeros = np.zeros(len(tileinfo))
    clm0, clm1, clm2, clm3, clm4 = Column(data=zeros, name='set_time'), \
            Column(data=zeros, name='set_exp'), Column(data=zeros, name='rise_time'), \
            Column(data=zeros, name='rise_exp'), Column(data=zeros, name='obsrv_time')
    tileinfo.add_columns([clm0, clm1, clm2, clm3, clm4])

    """Get rise and set time of tiles"""
    fxdbdy = FixedBody()
    for tile in tileinfo:
        fxdbdy._ra = degrees(str(tile['ra']))
        fxdbdy._dec = degrees(str(tile['dec']))
        fxdbdy.compute(observer)
        try:
            date_rise = observer.next_rising(fxdbdy, start=time.start)            
            date_set = observer.next_setting(fxdbdy, start=time.start)
            if (date_set > sunrise_time): starSet = sunrise_time
            else: starSet = date_set
            if (date_rise < date_set):
                tile['rise_time'], tile['rise_exp'] = date_rise, int((date_rise - time.start)/texpo)+1
                tile['set_time'], tile['set_exp'] = starSet, int((starSet - time.start)/texpo)
            else:
                date_rise = time.start
                tile['rise_time'],tile['rise_exp'] = date_rise, int((date_rise - time.start)/texpo)+1
                tile['set_time'],tile['set_exp'] = starSet, int((starSet - time.start)/texpo)
        except NeverUpError:
            tile['rise_time'], tile['rise_exp'] = Date(1.), 1
            tile['set_time'], tile['set_exp'] = Date(0.), 0
        except AlwaysUpError:
            date_rise = time.start
            date_set = sunrise_time
            tile['rise_time'], tile['rise_exp'] = float(date_rise), 0
            tile['set_time'], tile['set_exp'] = float(date_set), 10
    mask = tileinfo['set_time'] < time.end
    if mask.all(): time.end = max(tileinfo['set_time'])
    return tileinfo


def greedy(tileinfo, time):
    """The Greedy method.
    This method makes the schedule for observation by probability order"""

    #tileinfo = tileinfo[tileinfo['set_time'] > tileinfo['rise_time']]
    tileinfo.sort('prob')
    greedy_tiles = tileinfo[:1]
    for itr in range(int((time.end - time.start)/texpo)):
        mask = tileinfo['set_time'] >= time.start + (itr+1)*texpo
        tileinfo = tileinfo[mask]
        if tileinfo:
            tileinfo[0]['obsrv_time'] = julian_date(time.start + itr*texpo)
            greedy_tiles.add_row(tileinfo[0])
    return greedy_tiles


def settingArray(tileinfo, time):
    """
    Optimizes the observation considering setting of sky into account
    Among the possible candidates,it keeps the best canditates upto 
    each iteration and 'removes' the rest. This algorithm assumes that
    tiles are available for observation.
    """

    #tileinfo = tileinfo[tileinfo['set_time'] > tileinfo['rise_time']]
    observation_time = []
    for itr in range(int((time.end - time.start)/texpo)):
        mask = tileinfo['set_time'] <= time.start + (itr+1)*texpo
        if mask.any():
            setting_tiles = tileinfo[mask]
            tileinfo = tileinfo[mask^True]
            setting_tiles.sort('prob')
            setting_tiles['obsrv_time'] = julian_date(time.start + itr*texpo)
            [tileinfo.add_row(add_tile) for add_tile in setting_tiles[-opts.numb_telescopes:]]
    return tileinfo


def optimal(tile_info, time):
    """The optimal scheduling for the obsevation.
    The tiles are orderd in the optimal way"""

    sear = settingArray(tile_info, time)
    tileinfo = sear.copy()
    for tile in tileinfo:
        for indx in range(len(tileinfo)-1, 0, -1):
            if tileinfo[indx]['prob'] > tileinfo[indx-1]['prob'] and tileinfo[indx-1]['set_time'] >= time.start + (texpo * indx):
                tileinfo.insert_row(indx, vals= tileinfo[indx-1])
                tileinfo.remove_rows([indx+1])
    return tileinfo, sear


def multiepoch(optimaltiles, tileinfo, time):
    """Schedule the tiles for multi epoch observations"""

    factor = 1./(24*3600)
    timegap = opts.timegap * second
    expogap = int(timegap/texpo)
    timegap = expogap * texpo
    rising_mask = tileinfo['rise_time'] <= Date(Time(optimaltiles['obsrv_time'][0], format='jd').datetime)
    tileinfo = tileinfo[rising_mask]
    seting_mask = tileinfo['set_time'] > tileinfo['rise_time']
    tileinfo = tileinfo[seting_mask]
    if tileinfo:
        for i1 in range(len(tileinfo)):
            tileinfo['obsrv_time'][i1] = optimaltiles['obsrv_time'][i1] + texpo/float(opts.epoch) 
            if i1 > expogap: break
        [optimaltiles.add_row(row) for row in tileinfo[:expogap]]    
    obs_time_2nd = Column(data= optimaltiles['obsrv_time'] + (opts.timegap - texpo/opts.epoch) * factor, name='2nd_obsrv_time')
    optimaltiles.add_column(obs_time_2nd)
    optimaltiles.sort('obsrv_time')
    return optimaltiles


def risingPatch(tileinfo, observetory, time):
    """Use the scheduling method on rising patch"""

    optimaltiles = Table(names=('tile_id', 'ra', 'dec', 'prob', 'set_time', 'set_exp', 'rise_time', 'rise_exp', 'obsrv_time'), dtype=('i8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
    set_rise_tiles = setrise_tileinfo(tileinfo, time)
    tileinfo = set_rise_tiles.copy()
    ascii.write(tileinfo, 'setriseTiles.csv', format='csv', fast_writer=False, overwrite=True)

    if opts.epoch > 1: time.end = time.end - (opts.epoch-1) * opts.timegap * second
    # Optimal or SeAr or Greedy method
    while (time.start < time.end):
        time.start = time.start + texpo
        tileinfo_up, tileinfo_down = abovehorizon(tileinfo, time.start)
        if not tileinfo_up and not tileinfo_down: break
        if tileinfo_up:
            if opts.method == 'greedy': scheduled_tile = greedy(tileinfo_up, time)       
            if opts.method == 'sear': scheduled_tile = settingArray(tileinfo_up, time)
            if opts.method == 'optimal': scheduled_tile, _ = optimal(tileinfo_up, time)
            #if opts.method == 'all': 
            if not tileinfo_down:
                [optimaltiles.add_row(row) for row in scheduled_tile]
                break
            else:
                scheduled_tile['obsrv_time'] = julian_date(time.start)
                optimaltiles.add_row(scheduled_tile[0])
                remove_tile_indx = np.argmin(abs(tileinfo['tile_id'] - scheduled_tile[0]['tile_id']))
                tileinfo.remove_row(remove_tile_indx)
    if opts.epoch > 1:
        optimaltiles = multiepoch(optimaltiles, tileinfo, time)
    optimaltiles.remove_columns(['set_exp', 'rise_exp'])

    #print summary
    print "Observation start time is %s" %(str(Time(optimaltiles['obsrv_time'][0], format='jd').datetime))
    print "Observation starting point: ra= %f  dec= %f" %(optimaltiles['ra'][0], optimaltiles['dec'][0])
    print "The covered area is %f square degrees" %(len(optimaltiles) * observetory.fov_w * observetory.fov_h)
    print "Probability %f is covered in %d tiles" %(sum(optimaltiles['prob']), len(optimaltiles))
    ascii.write(optimaltiles, opts.output, format='csv', fast_writer=False, overwrite=True)
    return optimaltiles

risingPatch(tileinfo, observetory, time)
