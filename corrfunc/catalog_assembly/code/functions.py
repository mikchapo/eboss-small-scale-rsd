import numpy as np
import random

from astropy.table import Table
from scipy.integrate import quad


# H0 = 67.6
H0 = 100.
c = 299792.458
Omega_M = 0.31
Omega_L = 0.69

# Omega_M = 0.307115
# Omega_L = 1. - Omega_M


def comoving_dist_integral(z):
    return 1. / np.sqrt(Omega_M * (1. + z) * (1. + z) * (1. + z) + Omega_L)


def calc_position(ra, dec, zs):
    dH = c / H0
    r = dH * quad(comoving_dist_integral, 0, zs)[0]
    x = r * np.cos(dec*np.pi/180.) * np.cos(ra*np.pi/180.)
    y = r * np.cos(dec*np.pi/180.) * np.sin(ra*np.pi/180.)
    z = r * np.sin(dec*np.pi/180.)
    return x, y, z


def assemble_data(input_path, output_path, mock=False):
    cat = Table.read(input_path, hdu=1)

    file = open(output_path, 'w')
    file.write("#LRG_ID\t RA\t\t DEC\t\t ZS\t\t X\t\t Y\t\t Z\t\t WEIGHT_SYSTOT\t WEIGHT_CP\t WEIGHT_NOZ\t WEIGHT_FKP\n")

    if mock:
        id_field = "EBOSS_TARGET_ID"

        for i in range(len(cat['RA'])):
             if cat['Z'][i]>=0.6 and cat['Z'][i]<=1.:
                 x, y, z = calc_position(cat['RA'][i], cat['DEC'][i], cat['Z'][i])

                 file.write("%i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n" % 
                   (cat[id_field][i],cat['RA'][i],cat['DEC'][i],cat['Z'][i],x,y,z,1.,1.,1.,cat['WEIGHT_FKP'][i]))

    else:
        id_field = "LRG_ID"

        for i in range(len(cat['RA'])):
             x, y, z = calc_position(cat['RA'][i], cat['DEC'][i], cat['Z'][i])

             file.write("%i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n" % 
               (cat[id_field][i],cat['RA'][i],cat['DEC'][i],cat['Z'][i],x,y,z,cat['WEIGHT_SYSTOT'][i],cat['WEIGHT_CP'][i],cat['WEIGHT_NOZ'][i],cat['WEIGHT_FKP'][i]))

    file.close()


def assemble_rand(input_path, output_path, dilute=False, seed=1, threshold=1., mock=False):
    cat = Table.read(input_path, hdu=1)

    file = open(output_path, 'w')
    file.write("#RA\t\t DEC\t\t ZS\t\t X\t\t Y\t\t Z\t\t WEIGHT_SYSTOT\t WEIGHT_CP\t WEIGHT_NOZ\t WEIGHT_FKP\n")

    if dilute:
        random.seed(a=seed)
        if mock:
            for i in range(len(cat['RA'])):
                if random.random() < threshold and cat['Z'][i]>=0.6 and cat['Z'][i]<=1.:
                    x, y, z = calc_position(cat['RA'][i], cat['DEC'][i], cat['Z'][i])
                    file.write("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n" % (cat['RA'][i],cat['DEC'][i],cat['Z'][i],x,y,z,1.,1.,1.,cat['WEIGHT_FKP'][i]))

        else:
            for i in range(len(cat['RA'])):
                if random.random() < threshold:
                    x, y, z = calc_position(cat['RA'][i], cat['DEC'][i], cat['Z'][i])
                    file.write("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n" % (cat['RA'][i],cat['DEC'][i],cat['Z'][i],x,y,z,cat['WEIGHT_SYSTOT'][i],cat['WEIGHT_CP'][i],cat['WEIGHT_NOZ'][i],cat['WEIGHT_FKP'][i]))

    else:
        if mock:
            for i in range(len(cat['RA'])):
                if cat['Z'][i]>=0.6 and cat['Z'][i]<=1.:
                    x, y, z = calc_position(cat['RA'][i], cat['DEC'][i], cat['Z'][i])
                    file.write("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n" % (cat['RA'][i],cat['DEC'][i],cat['Z'][i],x,y,z,1.,1.,1.,cat['WEIGHT_FKP'][i]))

        else:
            for i in range(len(cat['RA'])):
                x, y, z = calc_position(cat['RA'][i], cat['DEC'][i], cat['Z'][i])
                file.write("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n" % (cat['RA'][i],cat['DEC'][i],cat['Z'][i],x,y,z,cat['WEIGHT_SYSTOT'][i],cat['WEIGHT_CP'][i],cat['WEIGHT_NOZ'][i],cat['WEIGHT_FKP'][i]))

    file.close()


def assemble_data_pip(input_path, output_path, bw_path, id_field="EBOSS_TARGET_ID"):
    cat = Table.read(input_path, hdu=1)

    file = open(output_path, 'w')
    file.write("#LRG_ID\t RA\t\t DEC\t\t ZS\t\t X\t\t Y\t\t Z\t\t WEIGHT_SYSTOT\t WEIGHT_CP\t WEIGHT_NOZ\t WEIGHT_FKP\t FIBER\t CLUSTERING\n")

    bw_file = open(bw_path, 'w')
    bw_file.write("#LRG_ID\t\t RA\t\t DEC\t\t ZS\t\t WEIGHT_BW\n")

    for i in range(len(cat['RA'])):
         x, y, z = calc_position(cat['RA'][i], cat['DEC'][i], cat['Z'][i])

         file.write("%i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %i\t %i\n" % 
           (cat[id_field][i],cat['RA'][i],cat['DEC'][i],cat['Z'][i],x,y,z,cat['WEIGHT_SYSTOT'][i],cat['WEIGHT_CP'][i],cat['WEIGHT_NOZ'][i],cat['WEIGHT_FKP'][i], cat['FIBER'][i], cat['CLUSTERING'][i]))

         bw_file.write("%i\t %f\t %f\t %f\t %s\n" % 
           (cat[id_field][i],cat['RA'][i],cat['DEC'][i],cat['Z'][i]," ".join(map(str, cat['WEIGHT_BW'][i]))))

    file.close()


def assemble_rand_mock(input_path, output_path, dilute=False, seed=1, threshold=1.):
    cat = np.loadtxt(input_path)

    file = open(output_path, 'w')
    file.write("#RA\t\t DEC\t\t ZS\t\t X\t\t Y\t\t Z\t\t WEIGHT_SYSTOT\t WEIGHT_CP\t WEIGHT_NOZ\t WEIGHT_FKP\n")

    # if dilute:
    #     random.seed(a=seed)
    #     for i in range(len(cat['RA'])):
    #         if random.random() < threshold:
    #             x, y, z = calc_position(cat['RA'][i], cat['DEC'][i], cat['Z'][i])
    #            file.write("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n" % (cat['RA'][i],cat['DEC'][i],cat['Z'][i],x,y,z,cat['WEIGHT_SYSTOT'][i],cat['WEIGHT_CP'][i],cat['WEIGHT_NOZ'][i],cat['WEIGHT_FKP'][i]))

    # else:
    for i in range(cat.shape[0]):
        if cat[i,2]>=0.6 and cat[i,2]<=1.:
            x, y, z = calc_position(cat[i, 0], cat[i, 1], cat[i, 2])
            file.write("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n" % (cat[i, 0],cat[i, 1],cat[i, 2],x,y,z,cat[i,3],cat[i,4],cat[i,5],cat[i,6]))

    file.close()
