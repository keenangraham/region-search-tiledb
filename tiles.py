import gzip
import io
import json
import numba as nb
import numpy as np
import os
import pandas as pd
import pybedtools
import requests
import tiledb

from collections import defaultdict
from pybedtools import BedTool


FILE_INDEX_START = int(os.environ.get('FILE_INDEX_START', 0))

GAP = 1000

ATTRIBUTE_NAME = 'region'

CHROMOSOME_TO_NUMBER = {
    'X': 23,
    'Y': 24,
    'M': 25
}


def create_region_array(array_name):
    dom = tiledb.Domain(
        tiledb.Dim(
            name='chr',
            domain=(1, 25),
            tile=25,
            dtype=np.int32,
        ),
        tiledb.Dim(
            name="genomic_position",
            domain=(0, 250000000),
            tile=1000,
            dtype=np.int32,
        ),
    )
    schema = tiledb.ArraySchema(
        domain=dom,
        sparse=True,
        attrs=[
            tiledb.Attr(
                name="region",
                dtype=np.dtype(
                    [
                        ('chrom', np.int32),
                        ('start', np.int32),
                        ('end', np.int32)
                    ]
                ),
                filters=tiledb.FilterList(
                    [
                        tiledb.GzipFilter()
                    ]
                )
            ),
        ],
        allows_duplicates=True,
    )
    tiledb.SparseArray.create(array_name, schema)


def get_remote_file(url):
    print('getting', url)
    return requests.get(url).content


def get_remote_bed_file(url):
    fh = gzip.open(
        io.BytesIO(
            get_remote_file(url)
        ),
        'rt'
    )
    return BedTool(fh)


def get_local_bed_file(filepath):
    return BedTool(filepath)


def parse_chromosome(chromosome):
    chromosome = chromosome.replace('chr', '')
    return CHROMOSOME_TO_NUMBER.get(
        chromosome,
        chromosome
    )


def get_regions_from_bed(bed):
    for feature in bed:
        try:
            parsed_chromosome = int(parse_chromosome(feature.chrom))
        except ValueError:
            print('skipping', feature.chrom)
            continue
        yield (
            parsed_chromosome,
            feature.start,
            feature.end
        )


def save_map_to_json_file(map_, name):
    print('Saving', name)
    with open(name, 'w') as f:
        json.dump(map_, f)


def load_map_from_json_file(name, key_to_int=False):
    print('loading', name)
    with open(name, 'r') as f:
        map_ = json.load(f)
    if key_to_int:
        return {int(k): v for k, v in map_.items()}
    return map_


def get_find_index_and_find_file_maps():
    current_index = FILE_INDEX_START
    accession_to_index = {}
    index_to_accession = {}

    def find_index(accession):
        nonlocal current_index
        nonlocal accession_to_index
        nonlocal index_to_accession
        if accession not in accession_to_index:
            accession_to_index[accession] = current_index
            index_to_accession[current_index] = accession
            current_index += 1
        return accession_to_index[accession]

    def find_file(index):
        nonlocal index_to_accession
        return index_to_accession.get(index)

    def save_maps(database):
        nonlocal accession_to_index
        nonlocal index_to_accession
        save_map_to_json_file(accession_to_index, f'{database}_a_to_i.json')
        save_map_to_json_file(index_to_accession, f'{database}_i_to_a.json')

    def load_maps(database):
        nonlocal accession_to_index
        nonlocal index_to_accession
        accession_to_index = load_map_from_json_file(f'{database}_a_to_i.json')
        index_to_accession = load_map_from_json_file(f'{database}_i_to_a.json', key_to_int=True)

    def get_maps():
        nonlocal accession_to_index
        nonlocal index_to_accession
        return accession_to_index, index_to_accession

    def set_maps(a_to_i, i_to_a):
        nonlocal accession_to_index
        nonlocal index_to_accession
        accession_to_index = a_to_i
        index_to_accession = i_to_a

    return find_index, find_file, save_maps, load_maps, get_maps, set_maps


def get_positions_for_region(start, end, gap):
    return (
        position
        for position in range(start, end, gap)
    )


def get_fill_value_for_region(filename, start, end):
    return (filename, start, end)


def get_data_for_region(file_index, region, gap):
   chrom, start, end = region
   fill_value = get_fill_value_for_region(file_index, start, end)
   return (
       (chrom, position, fill_value)
       for position in get_positions_for_region(start, end, gap)
   )


def generate_positions_from_regions(file_index, file_, gap):
    return (
        data
        for region in get_regions_from_bed(file_)
        for data in get_data_for_region(file_index, region, gap)
    )


def write_data_for_file(open_array, file_index, file_):
    print('generating positions')
    chroms, positions, data = zip(
        *generate_positions_from_regions(
            file_index,
            file_,
            GAP
        )
    )
    print('writing data')
    open_array[chroms, positions] = list(data)
    print('done')


@nb.njit
def query_intersects_result(result, start, end):
    # f1 -> start, f2 -> end.
    return (result['f2'] >= start) and (result['f1'] <= end)


# https://stackoverflow.com/a/58422691
@nb.njit
def filter_results(results, start, end, limit):
    length = 0
    filtered_length = 0
    while filtered_length < limit and length < results.size:
        if query_intersects_result(results[length], start, end):
            filtered_length += 1
        length += 1
    filtered_results = np.empty(filtered_length, dtype=results.dtype)
    filtered_length = 0
    for i in range(length):
        if query_intersects_result(results[i], start, end):
            filtered_results[filtered_length] = results[i]
            filtered_length += 1
    return filtered_results


def recover_file_from_result(result):
    idx, start, end = result
    return (find_file(idx), start, end)


def query_region(database, chrom, start, end, limit=25):
    with tiledb.SparseArray(database, 'r') as A:
        results = pd.unique(
            A[chrom, max(0, (start - GAP)):(end + 1)][ATTRIBUTE_NAME]
        )
        filtered_results = [
            recover_file_from_result(result)
            for result in filter_results(results, start, end, limit)
        ]
    return filtered_results


def get_filtered_file_results(results, start, end, limit):
    for result in filter_results(results, start, end, limit):
        filename = recover_file_from_result(result)[0]
        if filename:
            yield filename


def query_file(database, chrom, start, end, limit=25):
    with tiledb.SparseArray(database, 'r') as A:
        results = pd.unique(
            A[chrom, max(0, (start - GAP)):(end + 1)][ATTRIBUTE_NAME]
        )
        filtered_results = pd.unique(
            list(
                get_filtered_file_results(
                    results,
                    start,
                    end,
                    limit,
                )
            )
        )
    return filtered_results


def load(database, urls):
    with tiledb.SparseArray(database, 'w') as A:
        for i, url in enumerate(urls):
            print(i)
            file_ = get_remote_bed_file(url)
            accession = url.split('/')[-1].replace('.bed.gz', '')
            file_index = find_index(accession)
            print(accession, file_index)
            write_data_for_file(A, file_index, file_)


def load_local(database, files):
    with tiledb.SparseArray(database, mode='w') as A:
        for i, filepath in enumerate(files):
            print(i)
            file_ = get_local_bed_file(filepath)
            accession = filepath.split('/')[-1].replace('.bed.gz', '')
            file_index = find_index(filepath.split('/')[-1])
            print(accession, file_index)
            write_data_for_file(A, file_index, file_)


def clean_beds():
    pybedtools.cleanup(remove_all=True)


def clean(database):
    clean_beds()
    config = tiledb.Config(
        {
            "sm.consolidation.steps": 3,
            "sm.consolidation.mode": "fragment_meta"
        }
    )
    ctx = tiledb.Ctx(config)
    tiledb.consolidate(database, ctx=ctx)
    tiledb.vacuum(database)


find_index, find_file, save_maps, load_maps, get_maps, set_maps = get_find_index_and_find_file_maps()
