__author__ = 'ank'

from csv import reader
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
from pprint import PrettyPrinter
from itertools import izip, chain

root_path = '/n/projects/ank/2014/Chromo_Motility/'

Nup2_source_path = root_path+'HapN2_peaks.tsv'
Nup60_source_path = root_path+'HapN60_peaks.tsv'


def parse_source(source):
    """
    Parses a nup-containing file with first 24 lines removed and returns chromosome-specific parse

    :param source: source file
    :return: dict mapping chromosome to a list containing in the order: start, end, abs summit
    :rtype : dict
    """
    nup_pos_dict = defaultdict(list)
    with open(source, 'rb') as source_file:
        readr = reader(source_file, 'excel-tab')
        readr.next()
        for row in readr:
            nup_pos_dict[row[0]].append([int(row[1]), int(row[2]), int(row[3]), int(row[4]), (int(row[1])+int(row[2]))/2, float(row[7])])
    return dict(nup_pos_dict)


def get_distances(parse_dict):
    dist_dict = defaultdict(list)
    for chr, call_list in parse_dict.iteritems():
        for call1, call2 in izip(call_list[1:], call_list[:-1]):
            dist_dict[chr].append([call1[1]-call2[0], call1[0]-call2[1], call1[3]-call2[3], call1[4]-call2[4]])
    return dict(dist_dict)


def get_min_cross_distances(parse_dict1, parse_dict2):
    distance_dict = defaultdict(list)
    chr2overload = {}
    for chrom, call_list in parse_dict1.iteritems():
        overloading = defaultdict(list)
        re_overloading = defaultdict(list)
        search_list = parse_dict2[chrom]
        for idx_, call1 in enumerate(call_list):
            centroid_position = call1[4]
            idx, val = min(enumerate(search_list), key=lambda  x: abs(x[1][4]-centroid_position))
            call2 = search_list[idx]
            distance_dict[chr].append([abs(call1[1]-call2[0]), abs(call1[0]-call2[1]), abs(call1[3]-call2[3]), abs(call1[4]-call2[4])])
            overloading[idx].append((idx_, idx))
            re_overloading[idx_].append((idx_, idx))
        chr2overload[chrom] = (len([ value for value in overloading.itervalues() if len(value) > 1]),
                                len([ value for value in re_overloading.itervalues() if len(value) > 1]))

    return dict(distance_dict), chr2overload


def collapse_and_show_distances(dist_dict, limiter=None, log=True):
    array_hist_0 = np.array([ call[0] for call in chain(*dist_dict.itervalues())])
    array_hist_1 = np.array([ call[1] for call in chain(*dist_dict.itervalues())])
    array_hist_2 = np.array([ call[2] for call in chain(*dist_dict.itervalues())])
    array_hist_3 = np.array([ call[3] for call in chain(*dist_dict.itervalues())])

    plt.figure()

    plt.subplot(221)
    plt.title('closest side to side')
    plt.hist(array_hist_0, bins=50, log=log, range=limiter)

    plt.subplot(222)
    plt.title('furthest side to side')
    plt.hist(array_hist_1, bins=50, log=log, range=limiter)

    plt.subplot(223)
    plt.title('cemtroid to centroid')
    plt.hist(array_hist_2, bins=50, log=log, range=limiter)

    plt.subplot(224)
    plt.title('pic to peak')
    plt.hist(array_hist_3, bins=50, log=log, range=limiter)

    plt.show()

    means = np.mean(array_hist_0), np.mean(array_hist_1), np.mean(array_hist_2), np.mean(array_hist_3)
    print means[0]-means[1], means[2], means[3], means[3]-means[2]
    means = [m/1.6 for m in means]
    print means[0]-means[1], means[2], means[3], means[3]-means[2]


if __name__ == "__main__":
    pprinter = PrettyPrinter(indent=4)

    Nup2_parse = parse_source(Nup2_source_path)
    Nup60_parse = parse_source(Nup60_source_path)

    Nup2_distances = get_distances(Nup2_parse)
    Nup60_distance = get_distances(Nup60_parse)

    crossdists, overload_dict = get_min_cross_distances(Nup2_parse, Nup60_parse)

    collapse_and_show_distances(Nup2_distances)
    collapse_and_show_distances(Nup60_distance)
    collapse_and_show_distances(crossdists)
    collapse_and_show_distances(crossdists, limiter=(0,10000), log=False)
    collapse_and_show_distances(crossdists, limiter=(0,2000), log=False)
    pprinter.pprint(overload_dict)

    # pprinter.pprint()
    # pprinter.pprint()