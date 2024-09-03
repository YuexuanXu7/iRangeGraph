# iRangeGraph



## Quick Start

### Build

```bash
mkdir build && cd build && cmake .. && make
```

### Construct Index

#### parameters:

**`--data_path`**: The input data over which to build an index, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time. 
The data points should be already sorted in ascending order by the attribute.

**`--index_file`**: The constructed index will be saved to this file, in .bin format.

**`--M`**: The degree of the graph index.

**`--ef_construction`**: The size of result set during index building.

**`--threads`**: The number of threads for index building.


#### command:
```bash
./tests/buildindex --data_path [path to data points] --index_file [file path to save index] --M [integer] --ef_construction [integer] --threads [integer]
```


### Search For Single-Attribute

#### parameters:

**`--data_path`**: The data points over which the index is built, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time.
The data points should be already sorted in ascending order by the attribute.

**`--query_path`**: The query vectors, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the query one query point in a time.

**`--range_saveprefix`**: The path of folder where query range files will be saved.

**`--groundtruth_saveprefix`**: The path of folder where groundtruth files will be saved.

**`--index_file`**: The file path where the constructed index is saved, in .bin format. 

**`--result_saveprefix`**: The path of folder where result files will be saved.

**`--M`**: The degree of the graph index. It should equal the 'M' used for constructing index.

#### command:
```bash
./tests/search --data_path [path to data points] --query_path [path to query points] --range_saveprefix [folder path to save query ranges] --groundtruth_saveprefix [folder path to save groundtruth] --index_file [path of the index file] --result_saveprefix [folder path to save results] --M [integer]
```


### Search For Multi-Attribute

#### parameter:
**`--data_path`**:  The data points over which the index is built, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time.
The data points should not be sorted by the attribute.

**`--query_path`**: The query vectors, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the query one query point in a time.

**`--range_saveprefix`**: The path of folder where query range files will be saved.

**`--groundtruth_saveprefix`**: The path of folder where groundtruth files will be saved.

**`--index_file`**: The file path where the constructed index is saved, which is built with data points sorted by the first attribute.

**`--result_saveprefix`**: The path of folder where result files will be saved.

**`--attribute1_file`**: The path of the first attribute file, in .bin format. `n*sizeof(int)` bytes contain the first attributes of the data for one data point in a time.

**`--attribute2_file`**: The path of the second attribute file, in .bin format. `n*sizeof(int)` bytes contain the second attributes of the data for one data point in a time.

**`--M`**: The degree of the graph index. It should equal the 'M' used for constructing index by the first attribute.


#### command:
```bash
./tests/search_multi --data_path [path to data points] --query_path [path to query points] --range_saveprefix [folder path to save query ranges] --groundtruth_saveprefix [folder path to save groundtruth] --index_file [path of the index file] --result_saveprefix [folder path to save results] --attribute1 [path to first attributes] --attribute2 [path to second attributes] --M [integer]
```