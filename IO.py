__author__ = 'ank'

import numpy as np
from csv import reader, writer

source_path = '/home/ank/Documents/External_Predictions/Tamara_RNASeq/rpkm.csv'

output_path = '/home/ank/Documents/External_Predictions/Tamara_RNASeq/rpkm_cr.csv'
output_file = open(output_path, 'wb')


def import_data(anot_size, data_size):
    """
    """
    names = np.zeros((anot_size,1)).T
    table = np.zeros((data_size, 1)).T
    with open(source_path, 'rb') as source_file:
        rdr = reader(source_file, 'excel-tab')
        rdr.next()
        for row in rdr:
            rrw = [float(f) for f in row[anot_size : data_size + anot_size] ]
            table = np.concatenate( (table, np.array([rrw])))
            names = np.concatenate( (names, np.array([row[0 : anot_size]])))
    print table.shape
    return table[1: , :], names[1: , :]


def prima_filter(data, groups, filtration_fraction):
    """
    filters out samples where there are too much intra_sample noise
    """
    noise_list = []
    for group in groups:
        std = np.std(data[:, group], axis = 1)
        mean = np.mean(data[:, group], axis = 1)
        noise_list.append(np.reshape(std/mean,(data.shape[0],1)))

    noise = np.concatenate(tuple(noise_list), axis = 1)
    ret_set =  np.max(noise, axis = 1) < filtration_fraction
    print 'prima_filter_passed:', ret_set.nonzero()[0].__len__()
    return ret_set

def reshape(data, ms_array):
    """ """
    return np.reshape(ms_array, (data.shape[0],1))

def noise_removed_difference(data, group_1, group_2):
    """ """
    var_1 = np.var(data[:, group_1], axis = 1)
    mean_1 = np.mean(data[:, group_1], axis = 1)
    var_1, mean_1 = (reshape(data, var_1), reshape(data, mean_1))
    print var_1
    print mean_1

    var_2 = np.var(data[:, group_2], axis = 1)
    mean_2 = np.mean(data[:, group_2], axis = 1)
    var_2, mean_2 = (reshape(data, var_2), reshape(data, mean_2))
    print var_2
    print mean_2

    intra_v = np.mean(np.concatenate((var_1,var_2), axis=1), axis=1)
    print intra_v, intra_v.shape
    inter_v = np.var(np.concatenate((var_1,var_2), axis=1), axis=1)
    print inter_v, inter_v.shape
    return inter_v / intra_v

if __name__ == "__main__":
    data, names = import_data(1,9)

    groups = [[0,1,2],[3,4,5],[6,7,8]]
    prima_filter(data, groups, 0.05)

    print noise_removed_difference(data,groups[0],groups[1])

    tests = [(0,(1, 2)),(1, 2)]


