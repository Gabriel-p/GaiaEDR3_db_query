
import pathlib
import gzip

import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
# from uncertainties import ufloat
# from uncertainties import unumpy as unp
import matplotlib.pyplot as plt


# Path to the database
db_path = '/media/gabriel/backup/gabriel/GaiaEDR3/datafiles/'

# Input file
input_f = "clusters.ini"

# File that contains the regions delimited by each frame
txt_f = 'frame_ranges.txt'

max_mag = 19
# frame = 'equatorial'
frame = 'galactic'

max_box = 3.

# TODO: test if it matches Vizier at high DE


def main():
    """
    box_s is assumed to be in Equatorial coordinates
    """
    print(f"Using {frame} frame")

    clusters = pd.read_csv(input_f)

    # Read data about the regions occupied by each frame
    fdata = pd.read_csv(txt_f)

    # These are the points that determine the range of *all* the frames
    ra_min, ra_max = fdata['ra_min'], fdata['ra_max']
    dec_min, dec_max = fdata['dec_min'], fdata['dec_max']

    if frame == 'equatorial':
        xmin_fr, xmax_fr = ra_min, ra_max
        ymin_fr, ymax_fr = dec_min, dec_max
    elif frame == 'galactic':
        p1 = radec2lonlat(ra_min, dec_min)
        p2 = radec2lonlat(ra_min, dec_max)
        p3 = radec2lonlat(ra_max, dec_min)
        p4 = radec2lonlat(ra_max, dec_max)
        #
        xmin_fr = np.min(np.array([p1[0], p2[0], p3[0], p4[0]]), 0)
        xmax_fr = np.max(np.array([p1[0], p2[0], p3[0], p4[0]]), 0)
        ymin_fr = np.min(np.array([p1[1], p2[1], p3[1], p4[1]]), 0)
        ymax_fr = np.max(np.array([p1[1], p2[1], p3[1], p4[1]]), 0)

    # plt.subplot(121)
    # plt.scatter(ra_min[17], dec_min[17], marker='x')
    # plt.scatter(ra_min[17], dec_max[17], marker='x')
    # plt.scatter(ra_max[17], dec_min[17], marker='x')
    # plt.scatter(ra_max[17], dec_max[17], marker='x')
    # plt.subplot(122)
    # plt.scatter(xmin_fr[17], ymin_fr[17])
    # plt.scatter(xmin_fr[17], ymax_fr[17])
    # plt.scatter(xmax_fr[17], ymin_fr[17])
    # plt.scatter(xmax_fr[17], ymax_fr[17])
    # plt.show()
    # breakpoint()

    #
    for i, cluster in clusters.iterrows():
        cl_name = cluster['Name']
        print(f"\nQuerying cluster: {cl_name}...")

        # Cluster coordinates: center and box size
        c_ra, c_dec = cluster['ra'], cluster['dec']
        # c_lon, c_lat = cluster['GLON'], cluster['GLAT']
        # Size of box to query (in degrees)
        box_s_eq = cluster['box']

        data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl = findFrames(
            frame, max_box, c_ra, c_dec, box_s_eq, fdata, xmin_fr,
            ymin_fr, xmax_fr, ymax_fr)

        all_frames = query(
            frame, data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl)

        print("Adding xy coords and uncertainties")
        all_frames = xyCoords(all_frames)
        all_frames = uncertMags(all_frames)
        all_frames = all_frames.drop(columns=[
            'FG', 'e_FG', 'FBP', 'e_FBP', 'FRP', 'e_FRP'])

        msk = all_frames['Gmag'] < max_mag
        all_frames = all_frames[msk]
        print(f"{len(all_frames)} stars after magnitude cut <{max_mag}")

        print("Save to file")
        all_frames.to_csv('./out/' + cl_name + '.csv', index=False)


def radec2lonlat(ra, dec):
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lb = gc.transform_to('galactic')
    lon, lat = lb.l.value, lb.b.value
    return lon, lat


def findFrames(
    frame, max_box, c_ra, c_dec, box_s_eq, fdata, xmin_fr, ymin_fr,
        xmax_fr, ymax_fr):
    """
    """
    box_s_eq = min(max_box, box_s_eq)

    # Correct size in RA
    box_s_x = box_s_eq / np.cos(np.deg2rad(c_dec))
    xl, yl = box_s_x * .5, box_s_eq * .5

    if frame == 'galactic':
        # Correct length for the rotated frame (approx estimate)
        xl, yl = xl * .8, yl * .8

    # Limits of the cluster's region in Equatorial
    xmin_cl, xmax_cl = c_ra - xl, c_ra + xl
    ymin_cl, ymax_cl = c_dec - yl, c_dec + yl

    if frame == 'galactic':
        p1 = radec2lonlat(xmin_cl, ymin_cl)
        p2 = radec2lonlat(xmin_cl, ymax_cl)
        p3 = radec2lonlat(xmax_cl, ymin_cl)
        p4 = radec2lonlat(xmax_cl, ymax_cl)
        #
        xmin_cl = np.min(np.array([p1[0], p2[0], p3[0], p4[0]]), 0)
        xmax_cl = np.max(np.array([p1[0], p2[0], p3[0], p4[0]]), 0)
        ymin_cl = np.min(np.array([p1[1], p2[1], p3[1], p4[1]]), 0)
        ymax_cl = np.max(np.array([p1[1], p2[1], p3[1], p4[1]]), 0)

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

    data_in_files = list(fdata[frame_intersec]['filename'])
    print(f"Cluster is present in {len(data_in_files)} frames")

    # plt.scatter(xmin_cl, ymin_cl, c='k')
    # plt.scatter(xmin_cl, ymax_cl, c='k')
    # plt.scatter(xmax_cl, ymin_cl, c='k')
    # plt.scatter(xmax_cl, ymax_cl, c='k')
    # from matplotlib.pyplot import cm
    # color = iter(cm.rainbow(np.linspace(0, 1, len(data_in_files))))
    # for i, flag in enumerate(frame_intersec):
    #     if flag:
    #         c = next(color)
    #         plt.scatter(xmin_fr[i], ymin_fr[i], color=c, marker='x')
    #         plt.scatter(xmin_fr[i], ymax_fr[i], color=c, marker='x')
    #         plt.scatter(xmax_fr[i], ymax_fr[i], color=c, marker='x')
    #         plt.scatter(xmax_fr[i], ymin_fr[i], color=c, marker='x')
    # plt.show()

    return data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl


def query(frame, data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl):
    """
    """
    all_frames = []
    for i, file in enumerate(data_in_files):
        with gzip.open(db_path + file) as f:
            data = pd.read_csv(f, index_col=False)

            data_x, data_y = data['ra'], data['dec']
            if frame == 'galactic':
                data_x, data_y = radec2lonlat(data_x.values, data_y.values)
                data['lon'], data['lat'] = data_x, data_y

            mx = (data_x >= xmin_cl) & (data_x <= xmax_cl)
            my = (data_y >= ymin_cl) & (data_y <= ymax_cl)
            msk = (mx & my)
            if msk.sum() == 0:
                continue

            print(f"Frame {file} contains {msk.sum()} cluster stars")
            frame_i = data[msk]

            # frame_x, frame_y = frame_i['ra'].values, frame_i['dec'].values
            # if frame == 'galactic':
            #     frame_x, frame_y = frame_i['lon'].values, frame_i['lat'].values
            # if msk.sum():
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

    print(f"{len(all_frames)} stars queried")
    return all_frames


def xyCoords(data):
    """
    Replicate Vizier's '_x, _y' columns
    """
    ra, de = data['RA_ICRS'], data['DE_ICRS']
    # Center of the frame
    ra_0, de_0 = .5 * (np.max(ra) + np.min(ra)), .5 * (np.max(de) + np.min(de))
    coord1 = SkyCoord(ra_0, de_0, frame='icrs', unit='deg')
    coord2 = SkyCoord(ra.values, de.values, frame='icrs', unit='deg')
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
    # Sigmas are already squared here
    Zp_G, sigma_ZG_2 = 25.6873668671, 0.00000759
    Zp_BP, sigma_ZBP_2 = 25.3385422158, 0.000007785
    Zp_RP, sigma_ZRP_2 = 24.7478955012, 0.00001428

    I_G, e_IG = data['FG'].values, data['e_FG'].values
    I_BP, e_IBP = data['FBP'].values, data['e_FBP'].values
    I_RP, e_IRP = data['FRP'].values, data['e_FRP'].values

    data['Gmag'] = Zp_G + -2.5 * np.log10(I_G)
    BPmag = Zp_BP + -2.5 * np.log10(I_BP)
    RPmag = Zp_RP + -2.5 * np.log10(I_RP)
    data['BP-RP'] = BPmag - RPmag

    e_G = np.sqrt(sigma_ZG_2 + 1.179 * (e_IG / I_G)**2)
    data['e_Gmag'] = e_G
    e_BP = np.sqrt(sigma_ZBP_2 + 1.179 * (e_IBP / I_BP)**2)
    e_RP = np.sqrt(sigma_ZRP_2 + 1.179 * (e_IRP / I_RP)**2)
    data['e_BP-RP'] = np.sqrt(e_BP**2 + e_RP**2)

    # # Zero points for the G,BP,RP magnitudes.
    # Zp_G = ufloat(25.6873668671, 0.0027553202)
    # Zp_BP = ufloat(25.3385422158, 0.0027901700)
    # Zp_RP = ufloat(24.7478955012, 0.0037793818)

    # # Fluxes
    # I_G = unp.uarray(data['FG'], data['e_FG'])
    # I_BP = unp.uarray(data['FBP'], data['e_FBP'])
    # I_RP = unp.uarray(data['FRP'], data['e_FRP'])
    # # Magnitudes
    # mag_d = {
    #     'Gmag': Zp_G + -2.5 * unp.log10(I_G),
    #     'BPmag': Zp_BP + -2.5 * unp.log10(I_BP),
    #     'RPmag': Zp_RP + -2.5 * unp.log10(I_RP)}

    # e_mag = unp.std_devs(mag_d['Gmag'])
    # col_v = mag_d['BPmag'] - mag_d['RPmag']
    # e_col = unp.std_devs(col_v)

    # data['Gmag'] = unp.nominal_values(mag_d['Gmag'])
    # data['BP-RP'] = unp.nominal_values(col_v)
    # # Add errors for the magnitude and color
    # data['e_Gmag'] = e_mag
    # data['e_BP-RP'] = e_col

    return data


if __name__ == '__main__':
    # Create output dir if it does not exist
    path = pathlib.Path('./out')
    path.mkdir(parents=True, exist_ok=True)
    main()
