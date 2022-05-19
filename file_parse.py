
import pathlib
import gzip

import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy as unp
import matplotlib.pyplot as plt


# Path to the database
db_path = '/media/gabriel/USB1TB/GaiaEDR3/datafiles/'

# File that contains the regions delimited by each frame
txt_f = 'frame_ranges.txt'

# TODO: test if it matches Vizier at high DE


def main():
    """
    """
    frame = 'galactic'  # 'equatorial'
    print("Using {frame} frame")

    cl_name = "NGC752"
    # Cluster coordinates: center and box size
    cx, cy = 136.9587510, -23.2884
    # Size of box to query (in arcmin)
    box_s = 360

    all_frames = query(frame, cl_name, cx, cy, box_s)

    all_frames = xyCoords(all_frames)
    all_frames = uncertMags(all_frames)
    all_frames = all_frames.drop(columns=[
        'FG', 'e_FG', 'FBP', 'e_FBP', 'FRP', 'e_FRP'])
    all_frames.to_csv('./out/' + cl_name + '.csv', index=False)


def query(frame, cl_name, cx, cy, box_s):
    """
    """
    print(f"Querying cluster: {cl_name}...")
    # To degrees
    box_s /= 60.

    if frame == 'equatorial':
        # Correct size in RA
        box_s_x = box_s / np.cos(np.deg2rad(cy))
    elif frame == 'galactic':
        box_s_x = box_s
    # Limits of the cluster's region
    xmin_cl, xmax_cl = cx - box_s_x * .5, cx + box_s_x * .5
    ymin_cl, ymax_cl = cy - box_s * .5, cy + box_s * .5

    # Read data about the regions occupied by each frame
    fdata = pd.read_csv(txt_f)

    # These are the points that determine the range of *all* the frames
    ra_min, ra_max = fdata['ra_min'], fdata['ra_max']
    dec_min, dec_max = fdata['dec_min'], fdata['dec_max']

    if frame == 'equatorial':
        xmin_fr, xmax_fr = ra_min, ra_max
        ymin_fr, ymax_fr = dec_min, dec_max
    elif frame == 'galactic':
        gc = SkyCoord(ra=ra_min * u.degree, dec=dec_min * u.degree)
        lb_min = gc.transform_to('galactic')
        xmin_fr, ymin_fr = lb_min.l.value, lb_min.b.value
        gc = SkyCoord(ra=ra_max * u.degree, dec=dec_max * u.degree)
        lb_max = gc.transform_to('galactic')
        xmax_fr, ymax_fr = lb_max.l.value, lb_max.b.value

    # Identify which frames contain the cluster region
    p1_in = (xmin_fr <= xmin_cl) & (xmax_fr >= xmin_cl) &\
        (ymin_fr <= ymin_cl) & (ymax_fr >= ymin_cl)
    p2_in = (xmin_fr <= xmin_cl) & (xmax_fr >= xmin_cl) &\
        (ymin_fr <= ymax_cl) & (ymax_fr >= ymax_cl)
    p3_in = (xmin_fr <= xmax_cl) & (xmax_fr >= xmax_cl) &\
        (ymin_fr <= ymin_cl) & (ymax_fr >= ymin_cl)
    p4_in = (xmin_fr <= xmax_cl) & (xmax_fr >= xmax_cl) &\
        (ymin_fr <= ymax_cl) & (ymax_fr >= ymax_cl)
    # Frames that overlap with cluster's region
    frame_intersec = p1_in | p2_in | p3_in | p4_in

    # plt.scatter(xmin_cl, ymin_cl, c='k')
    # plt.scatter(xmin_cl, ymax_cl, c='k')
    # plt.scatter(xmax_cl, ymin_cl, c='k')
    # plt.scatter(xmax_cl, ymax_cl, c='k')
    # cols, mrks = ('r', 'g', 'b', 'cyan', 'orange'), ('s', 'x', 'v', '^', 'D')
    # j = 0
    # for i, flag in enumerate(frame_intersec):
    #     if flag:
    #         plt.scatter(xmin_fr[i], ymin_fr[i], c=cols[j], marker=mrks[j])
    #         plt.scatter(xmin_fr[i], ymax_fr[i], c=cols[j], marker=mrks[j])
    #         plt.scatter(xmax_fr[i], ymax_fr[i], c=cols[j], marker=mrks[j])
    #         plt.scatter(xmax_fr[i], ymin_fr[i], c=cols[j], marker=mrks[j])
    #         j += 1
    #     if j == 4:
    #         break
    # plt.show()

    data_in_files = list(fdata[frame_intersec]['filename'])
    print(f"Cluster is present in {len(data_in_files)} frames")

    all_frames = []
    for i, file in enumerate(data_in_files):
        with gzip.open(db_path + file) as f:
            data = pd.read_csv(f, index_col=False)

            data_x, data_y = data['ra'], data['dec']
            if frame == 'galactic':
                gc = SkyCoord(ra=data_x * u.degree, dec=data_y * u.degree)
                xy = gc.transform_to('galactic')
                data_x, data_y = xy.l.value, xy.b.value
                data['lon'], data['lat'] = data_x, data_y

            mx = (data_x >= xmin_cl) & (data_x <= xmax_cl)
            my = (data_y >= ymin_cl) & (data_y <= ymax_cl)
            mxy = (mx & my)
            frame_i = data[mxy]
            print(f"Frame {file} contains {len(frame_i)} stars"
                  + " from the cluster region")

            # frame_x, frame_y = frame_i['ra'].values, frame_i['dec'].values
            # if frame == 'galactic':
            #     frame_x, frame_y = frame_i['lon'].values, frame_i['lat'].values
            # if mxy.sum():
            #     if len(frame_i) > 500:
            #         idxs = np.random.choice(len(frame_i), 500, replace=False)
            #         plt.scatter(frame_x[idxs], frame_y[idxs], alpha=.5,
            #                     marker='x')
            #     else:
            #         plt.scatter(frame_x, frame_y, alpha=.5, marker='x')

            all_frames.append(frame_i)
    # plt.show()

    all_frames = pd.concat(all_frames)

    all_frames = all_frames.rename(columns={
        'designation': 'EDR3Name', 'ra': 'RA_ICRS', 'dec': 'DE_ICRS',
        'parallax': 'Plx', 'parallax_error': 'e_Plx',
        'pmra': 'pmRA', 'pmra_error': 'e_pmRA', 'lat': 'GLAT',
        'pmdec': 'pmDE', 'pmdec_error': 'e_pmDE', 'lon': 'GLON',
        'phot_g_mean_flux': 'FG', 'phot_g_mean_flux_error': 'e_FG',
        'phot_bp_mean_flux': 'FBP', 'phot_bp_mean_flux_error': 'e_FBP',
        'phot_rp_mean_flux': 'FRP', 'phot_rp_mean_flux_error': 'e_FRP',
        'dr2_radial_velocity': 'RVDR2', 'dr2_radial_velocity_error': 'e_RVDR2'}
    )

    print(f"Writing {len(all_frames)} stars to file")

    return all_frames


def xyCoords(data):
    """
    Replicate Vizier's '_x, _y' columns
    """
    ra, de = data['RA_ICRS'], data['DE_ICRS']
    # Center of the frame
    ra_0, de_0 = .5 * (np.max(ra) + np.min(ra)), .5 * (np.max(de) + np.min(de))
    coord1 = SkyCoord(ra_0, de_0, frame='icrs', unit='deg')
    coord2 = SkyCoord(ra, de, frame='icrs', unit='deg')
    r = coord1.separation(coord2)
    theta = coord1.position_angle(coord2)
    _x = r * np.cos(theta)
    _y = r * np.sin(theta)

    data['_x'], data['_y'] = _x.value, _y.value
    # data = data.assign(_x=_x)
    # data = data.assign(_y=_y)

    return data


def uncertMags(data):
    """
    # Gaia EDR3 zero points:

    https://www.cosmos.esa.int/web/gaia/edr3-passbands

    """
    # Zero points for the G,BP,RP magnitudes.
    Zp_G = ufloat(25.6873668671, 0.0027553202)
    Zp_BP = ufloat(25.3385422158, 0.0027901700)
    Zp_RP = ufloat(24.7478955012, 0.0037793818)

    # Fluxes
    I_G = unp.uarray(data['FG'], data['e_FG'])
    I_BP = unp.uarray(data['FBP'], data['e_FBP'])
    I_RP = unp.uarray(data['FRP'], data['e_FRP'])
    # Magnitudes
    mag_d = {
        'Gmag': Zp_G + -2.5 * unp.log10(I_G),
        'BPmag': Zp_BP + -2.5 * unp.log10(I_BP),
        'RPmag': Zp_RP + -2.5 * unp.log10(I_RP)}

    e_mag = unp.std_devs(mag_d['Gmag'])
    col_v = mag_d['BPmag'] - mag_d['RPmag']
    e_col = unp.std_devs(col_v)

    data['Gmag'] = unp.nominal_values(mag_d['Gmag'])
    data['BP-RP'] = unp.nominal_values(col_v)
    # Add errors for the magnitude and color
    data['e_Gmag'] = e_mag
    data['e_BP-RP'] = e_col

    return data


if __name__ == '__main__':
    # Create output dir if it does not exist
    path = pathlib.Path('./out')
    path.mkdir(parents=True, exist_ok=True)
    main()
