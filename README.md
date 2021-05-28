# Region search with tileDB
Ingest and query genomic intervals from multiple BED files.


## Ingest
```python
>>> import glob
>>> import tiles
# Get paths to local BED files.
>>> files = sorted(glob.glob('/path/to/beds/*.bed.gz'))
>>> files[:5]
[
    'ENCFF001SOF.bed.gz',
    'ENCFF001SOH.bed.gz',
    'ENCFF001SOM.bed.gz',
    'ENCFF001SON.bed.gz',
    'ENCFF001SOO.bed.gz',
]
# Create chromosome by genomic position database.
>>> tiles.create_region_array('regions')
# Load first five BEDs into tileDB.
>>> tiles.load_local('regions', files[:5])
...
generating positions
writing data
done
...
# Save accession to file index maps.
>>> tiles.save_maps('regions')
# Consolidate and vacuum database.
>>> tiles.clean('regions')
```

## Query
```python
>>> import tiles
# Load accession to file index maps.
>> tiles.load_maps('regions')
# Get first five intervals around POMC gene on chromosome two.
>>> tiles.query_region('regions', 2, 25132492, 25192278, limit=5)
[
    ('ENCFF001SOF.bed.gz', 25139665, 25139815),
    ('ENCFF001SOF.bed.gz', 25141485, 25141635),
    ('ENCFF001SOF.bed.gz', 25141625, 25141775),
    ('ENCFF001SOF.bed.gz', 25142065, 25142215),
    ('ENCFF001SOF.bed.gz', 25142665, 25142815),
]
# Get total number of intervals in that region.
>>> len(tiles.query_region('regions', 2, 25132492, 25192278, limit=1000))
52
# Get files with intervals in that region.
>>> tiles.query_file('regions', 2, 25132492, 25192278)
[
    'ENCFF001SOF.bed.gz',
    'ENCFF001SOH.bed.gz',
    'ENCFF001SOM.bed.gz',
]
```
