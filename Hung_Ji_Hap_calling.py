__author__ = 'ank'

from csv import reader, writer
from collections import defaultdict, OrderedDict
import numpy as np
from matplotlib import pyplot as plt
from pprint import PrettyPrinter
from itertools import izip, chain

root_path = '/n/projects/ank/2014/Chromo_Motility/'

Nup2_source_path = root_path+'AN2_peaks.tsv'
Nup60_source_path = root_path+'AN60_peaks.tsv'

output_path = root_path+'AN2_vs_N60.tsv'
filtered_output = root_path+'AN2_vs_N60_filtered.tsv'
better_output = root_path+'AN2_vs_N60_better.tsv'


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
        namerow = readr.next()
        namemap = [ namerow[1], namerow[2], namerow[3], namerow[4], 'centroid', namerow[7], namerow[9], namerow[6]]
        for row in readr:
            nup_pos_dict[row[0]].append([int(row[1]), int(row[2]), int(row[3]), int(row[4]), (int(row[1])+int(row[2]))/2, float(row[7]), row[9], float(row[6])])
    return dict(nup_pos_dict), namemap


def composite_writer(out_path, Closest_dict, descript_dict1, descript_dict2, namemap):
    with open(out_path, 'wb') as out_file:
        out_writer = writer(out_file, 'excel-tab')
        out_writer.writerow(['chr','centroid_dist','peak_dist'] + namemap)
        for chrom, closeness_dict in Closest_dict.iteritems():
            out_writer.writerow([chrom])
            for dict2_idx, dict1_idx2dists in closeness_dict.iteritems():
                out_writer.writerow(['','----','----'] + descript_dict2[chrom][dict2_idx])
                for dict1_idx, dists in dict1_idx2dists.iteritems():
                    out_writer.writerow(['']+dists + descript_dict1[chrom][dict1_idx])
    return None


def b_dist_writer(out_path, better_dict, namemap):
    with open(out_path, 'wb') as out_file:
        out_writer = writer(out_file, 'excel-tab')
        out_writer.writerow(['chr','centroid_dist','peak_dist'] + namemap[2:]+namemap[2:])
        for line in better_dict:
            out_writer.writerow(line)
    return None


def get_distances(parse_dict):
    dist_dict = defaultdict(list)
    for chr, call_list in parse_dict.iteritems():
        for call1, call2 in izip(call_list[1:], call_list[:-1]):
            dist_dict[chr].append([call1[1]-call2[0], call1[0]-call2[1], call1[3]-call2[3], call1[4]-call2[4]])
    return dict(dist_dict)


def get_min_cross_distances(parse_dict1, parse_dict2):
    distance_dict = defaultdict(list)
    chr2overload = {}
    chr2memory = {}
    for chrom, call_list in parse_dict1.iteritems():
        overloading = defaultdict(list)
        re_overloading = defaultdict(list)
        overload_memory = defaultdict(OrderedDict)
        search_list = parse_dict2[chrom]
        for idx_, call1 in enumerate(call_list):
            centroid_position = call1[4]
            idx, val = min(enumerate(search_list), key=lambda  x: abs(x[1][4]-centroid_position))
            call2 = search_list[idx]
            distmap  = [abs(call1[1]-call2[0]), abs(call1[0]-call2[1]), abs(call1[3]-call2[3]), abs(call1[4]-call2[4])]
            distance_dict[chr].append(distmap)
            overloading[idx].append((idx_, idx))
            re_overloading[idx_].append((idx_, idx))
            overload_memory[idx][idx_] = distmap[2:]
        chr2overload[chrom] = (len([ value for value in overloading.itervalues() if len(value) > 1]),
                                len([ value for value in re_overloading.itervalues() if len(value) > 1]))
        chr2memory[chrom] = OrderedDict(sorted(overload_memory.items(), key=lambda t: t[0]))

    return dict(distance_dict), chr2overload,  OrderedDict(sorted(chr2memory.items(), key=lambda t: t[0]))


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


def filter_dist_dict(range, Closest_dict, descript_dict1, descript_dict2,):
    """
    Takes a ch2memory dict and outputs three dictionaries:
        one with the chrom to disclist of the cross-dist to memory,
        tuple of dictionaries: the parse dicts of the elements that failed the tests
        tuple of dictionaries: the parse dicts of the elements that suceeded in the tests

    Filtering procedure: filters on intercentroid distance, then on the interpeak distance and finally on the significant
        intercoverage of the probabilities of presence of the peaks.

    """
    Filtered_closest_dict =  {}
    plus_descript_dict1 = defaultdict(list)
    plus_descript_dict2 = defaultdict(list)
    minus_descript_dict1 = defaultdict(list)
    minus_descript_dict2 = defaultdict(list)
    better_dist_dict = []
    for chrom, closeness_dict in Closest_dict.iteritems():
        overload_memory = defaultdict(OrderedDict)
        for dict2_idx, dict1_idx2dists in closeness_dict.iteritems():
            for dict1_idx, dists in dict1_idx2dists.iteritems():
                if dists[0]<range or dists[1]<range or abs(dists[0]-range) < 0.25*(descript_dict1[chrom][dict1_idx][2] +
                                                                                  descript_dict2[chrom][dict2_idx][2])/2.0:
                    overload_memory[dict2_idx][dict1_idx] = dists
                    plus_descript_dict2[chrom].append(descript_dict2[chrom][dict2_idx])
                    plus_descript_dict1[chrom].append(descript_dict1[chrom][dict1_idx])
                    better_dist_dict.append([chrom]+dists+
                                            descript_dict1[chrom][dict1_idx][2:]+descript_dict2[chrom][dict2_idx][2:])
                else:
                    minus_descript_dict2[chrom].append(descript_dict2[chrom][dict2_idx])
                    minus_descript_dict1[chrom].append(descript_dict1[chrom][dict1_idx])

        Filtered_closest_dict[chrom] = OrderedDict(sorted(overload_memory.items(), key=lambda t: t[0]))

    return better_dist_dict, Filtered_closest_dict, plus_descript_dict1, plus_descript_dict2, plus_descript_dict1, plus_descript_dict2



if __name__ == "__main__":
    pprinter = PrettyPrinter(indent=4)

    Nup2_parse, namemap = parse_source(Nup2_source_path)
    Nup60_parse, _ = parse_source(Nup60_source_path)

    Nup2_distances = get_distances(Nup2_parse)
    Nup60_distance = get_distances(Nup60_parse)

    crossdists, overload_dict, chr2mem = get_min_cross_distances(Nup2_parse, Nup60_parse)

    collapse_and_show_distances(Nup2_distances)
    collapse_and_show_distances(Nup60_distance)
    collapse_and_show_distances(crossdists)
    collapse_and_show_distances(crossdists, limiter=(0,10000), log=False)
    collapse_and_show_distances(crossdists, limiter=(0,2000), log=False)


    better_dist_dict, filtered_dists, Nup2_pos_dists, Nup60_pos_dists, Nup2_neg_dists, Nup60_neg_dists = filter_dist_dict(500, chr2mem, Nup2_parse, Nup60_parse)

    collapse_and_show_distances(Nup2_pos_dists)
    collapse_and_show_distances(Nup60_pos_dists)

    collapse_and_show_distances(Nup2_neg_dists)
    collapse_and_show_distances(Nup60_neg_dists)

    pprinter.pprint(overload_dict)

    composite_writer(output_path, chr2mem, Nup2_parse, Nup60_parse, namemap)
    composite_writer(filtered_output, filtered_dists, Nup2_parse, Nup60_parse, namemap)
    b_dist_writer(better_output,better_dist_dict, namemap)

    # pprinter.pprint()
    # pprinter.pprint()