__author__ = 'ank'

import numpy as np
from csv import reader, writer

source_path = '/home/ank/Documents/External_Predictions/Tamara_RNASeq/rpkm.csv'

output_path = '/home/ank/Documents/External_Predictions/Tamara_RNASeq/rpkm_cr.csv'
output_file = open(output_path,'wb')

def import_data(anot_size, data_size):
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

    noise = np.concatenate(tuple(noise_list), axis=1)
    ret_set =  np.max(noise, axis=1) < filtration_fraction
    print 'prima_filter_passed:', ret_set.nonzero()[0].__len__()
    return ret_set

def noise_removed_difference(group_1,group_2):
    std_1 = np.std(data[:, group_1], axis = 1)
    mean_1 = np.mean(data[:, group_1], axis = 1)

    std_2 = np.std(data[:, group_2], axis = 1)
    mean_2 = np.mean(data[:, group_2], axis = 1)

    


if __name__ == "__main__":
    data, names = import_data(1,9)

    groups = [[0,1,2],[3,4,5],[6,7,8]]
    prima_filter(data, groups, 0.05)

    tests = [(0,(1, 2)),(1, 2)]
