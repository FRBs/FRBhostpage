""" Methods to generate the FRB and Host galaxy HTML tables"""

import pandas
import os
import numpy as np

from frb import frb
from frb.galaxies import frbgalaxy
from frb.galaxies import utils

from IPython import embed

def get_values(tbl, tbl_units, idx, properties, formats, scale_dict=None,
               units_dict=None, suppress_dict=None):
    values, errors, units = [], [], []
    for key in properties.keys():
        # Proceeed
        value = tbl.iloc[idx][key]
        scale = 1. if scale_dict is None else scale_dict[key]
        if key in formats.keys():
            values.append(format(value / scale, formats[key]))
        else:
            raise IOError("You must add a format for {}!".format(key))
        # Error?
        if hasattr(tbl.iloc[idx], key + '_err'):
            error = tbl.iloc[idx][key + '_err']
            errors.append(format(error / scale, formats[key]))
        else:
            errors.append('--')
        # Unit
        unit = tbl_units[key] if units_dict is None else units_dict[key]
        units.append(unit)

    # Return
    return values, errors, units

def build_frbs(out_path='./html_tables'):
    """
    """
    # Load
    frbs_tbl, tbl_units = frb.build_table_of_frbs()

    # Properties for the Table
    frb_properties = dict(DM='DM_FRB', DMISM='DM_ISM', RM='RM_FRB', fluence='Fluence')
    frb_prop = [item for key, item in frb_properties.items()]

    # Suppress?
    frb_suppress = dict(FRB200430=['DM'])
    for key in frb_suppress:
        idx = np.where(frbs_tbl['FRB'] == key)[0][0]
        for item in frb_suppress[key]:
            frbs_tbl.at[idx, item] =  np.nan
            # Error
            if item+'_err' in frbs_tbl.keys():
                frbs_tbl.at[idx, item+'_err'] = np.nan

    # Formatting
    frb_formats = dict(DM='.1f', DMISM='.1f', RM='.1f', fluence='.1f')

    # Loop me
    for frb_idx in range(len(frbs_tbl)):
        # Quantity to get us started
        frb_tbl = pandas.DataFrame(dict(Quantity=frb_prop))

        # The rest
        values, errors, units = get_values(frbs_tbl, tbl_units, frb_idx,
                                           frb_properties, frb_formats)
        # Add em in
        frb_tbl['Measured Value'] = values
        frb_tbl['Measured Error'] = errors
        frb_tbl['Units'] = units

        # To HTML and disk
        filename = os.path.join(out_path, frbs_tbl.iloc[frb_idx].FRB+'.html')
        frb_tbl.to_html(open(filename, 'w'))

def build_hosts(out_path='./html_tables'):
    hosts_tbl, host_tbl_units = utils.build_table_of_hosts()

    host_formats = dict(Mstar='0.2f', SFR_nebular='0.2f', EBV_photom='0.2f', M_r='0.2f', age_mass='0.2f', reff_kpc='0.2f', physical='0.2f',
                        Halpha='0.2f', Hbeta='0.2f', [OIII] 5007='0.2f', '[NII] 6584'='0.2f')
    host_scale = dict(Mstar=1e9, SFR_nebular=1., EBV_photom=1., M_r=1., age_mass=1.e-3, reff_kpc=1., physical=1.,
                      Halpha=1e-16, Hbeta=1e-16, '[OIII] 5007'=1e-16, '[NII] 6584'=1e-16)
    host_units = dict(Mstar='Msun', SFR_nebular='Msun/yr', EBV_photom='mag', M_r='mag', age_mass='Gyr', reff_kpc='kpc', physical='kpc',
                      Halpha='erg/s/cm^2', Hbeta='erg/s/cm^2','[OIII] 5007'='erg/s/cm^2', '[NII] 6584'='erg/s/cm^2')

    # Derived
    host_derived_properties = dict(Mstar='Stellar Mass', SFR_nebular='SFR', EBV_photom='<i>E(B-V)</i>', M_r='Absolute r-band mag.', u-r='rest-frame u-r color',
                                    age_mass='Mass-weighted age', reff_kpc='Half-light radius', physical='Offset from FRB to host galaxy')
    host_emission_properties = dict(Halpha='H&alpha;', Hbeta='H&beta;','[OIII] 5007'='[OIII] 5007', '[NII] 6584'='[NII] 6584')

    # Photometry -- *way* easier with JSON files
    suffix = 'photom'
    for host_idx in range(len(hosts_tbl)):
        # Load JSON
        host = frbgalaxy.FRBHost.by_name(hosts_tbl.iloc[host_idx].Host[2:])
        telescopes, filters, values, errors = [], [], [], []
        for key in host.photom.keys():
            if 'err' in key:
                continue
            if key in ['EBV']:
                continue
            # Parse
            items = key.split('_')
            telescopes.append('-'.join(items[:-1]))
            filters.append(items[-1])
            values.append(format(host.photom[key], '0.2f'))
            errors.append(format(host.photom[key+'_err'], '0.2f'))

        host_tbl = pandas.DataFrame(dict(Telescope=telescopes,
                                         Filter=filters,
                                         Value=values,
                                         Error=errors))
        # To HTML and disk
        filename = os.path.join(out_path, hosts_tbl.iloc[host_idx].Host + '_{}.html'.format(suffix))
        host_tbl.to_html(open(filename, 'w'))

    # Loop me
    for properties, suffix in zip([host_emission_properties, host_derived_properties],
                          ['emission', 'derived']):
        host_prop = [item for key, item in properties.items()]
        # Loop on Table
        for host_idx in range(len(hosts_tbl)):
            # Quantity to get us started
            host_tbl = pandas.DataFrame(dict(Quantity=host_prop))
            # The rest
            values, errors, units = get_values(hosts_tbl, host_tbl_units,
                                               host_idx, properties,
                                               host_formats,
                                               scale_dict=host_scale,
                                               units_dict=host_units)
            # Add em in
            host_tbl['Measured Value'] = values
            host_tbl['Measured Error'] = errors
            host_tbl['Units'] = units

            # To HTML and disk
            filename = os.path.join(out_path, hosts_tbl.iloc[host_idx].Host+'_{}.html'.format(suffix))
            host_tbl.to_html(open(filename, 'w'))



def main(inflg='all', options=None):

    if inflg == 'all':
        flg = np.sum(np.array( [2**ii for ii in range(25)]))
    else:
        flg = int(inflg)

    # FRBs
    if flg & (2**0):
        build_frbs()

    # Hosts
    if flg & (2**1):
        build_hosts()

# Command line execution
if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1:
        flg = 0
        #flg += 2**0   # FRB tables
        flg += 2**1   # Host tables
        #flg_fig += 2**1   # FRB 190102
        #flg_fig += 2**16   # Check impacts
    else:
        flg = sys.argv[1]

    main(inflg=flg)
