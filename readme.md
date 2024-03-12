# Preference Shapley


## This is a dummy Repo for the paper of Computing Shapley Values in Preference Queries as the paper is still under review.

## Compile

```
source compile.sh
```

If you only want to compile some certain files

```
export LD_LIBRARY_PATH=$(pwd)/lib:$LD_LIBRARY_PATH
g++ main_compute_sv_2d.cpp compute_sv_2d.cpp local_partition.cpp intersect_lp_solve.cpp intersect_sweep_line.cpp global_partition_2d.cpp -o main_compute_sv_2d -L ./lib -lcdd -llpsolve55 -pthread -lClp -lCoinUtils
g++ main_compute_partition_hd.cpp local_partition.cpp intersect_lp_solve.cpp global_partition.cpp dominance.cpp coalition.cpp  global_partition_2d.cpp -o main_compute_partition_hd -L ./lib -lcdd -llpsolve55 -pthread -lClp -lCoinUtils
g++ main_hyperplane_arrangement.cpp local_partition.cpp global_partition.cpp intersect_lp_solve.cpp  -o main_hyperplane_arrangement -L ./lib -lcdd -llpsolve55 -pthread -lClp -lCoinUtils
g++ main_merge_partition.cpp local_partition.cpp global_partition.cpp intersect_lp_solve.cpp -o main_merge_partition -L ./lib -lcdd -llpsolve55 -pthread -lClp -lCoinUtils
g++ main_permutation.cpp local_partition.cpp intersect_lp_solve.cpp global_partition.cpp permutation.cpp -o main_permutation -L ./lib -lcdd -llpsolve55 -pthread -lClp -lCoinUtils
g++ main_compute_skyline.cpp local_partition.cpp intersect_lp_solve.cpp -o main_compute_skyline -L ./lib -lcdd -llpsolve55 -pthread -lClp -lCoinUtils
g++ gen_data.cpp -o gen_data
```

### requirements:
Already included in this repo

[lpsolve](https://lpsolve.sourceforge.net/5.5/)

[cddlib](https://github.com/cddlib/cddlib/tree/master)

[clp](https://github.com/coin-or/Clp)

### Virtual Environment

Install sage 10.0 using conda [Instruction](https://doc.sagemath.org/html/en/installation/conda.html#sec-installation-conda-develop)

```
conda activate sage
pip install vegas tqdm polytope
```

## Integrity Test
```
python test_function.py
```

## Dataset Preparation

### On Generated Dataset
We provide several datasets at ```./data/``` to run

If you want to generate some new datasets

#### Generate New Dataset

To generate new dataset, use the following command:

```
./gen_data [total number of points] [dimensionality] [data distribution] [partition type] [number of data owners] [number of datasets] [-seeds xxxx]
```

- [data distribution]: 1. uniform 2. corr 3. anti

- [partition type]: 1: random 4. layer

- [number of data owners]: [number of data owners] should be smaller than [total number of points]

- [number of datasets]: the number of datasets generated of this type

- [-seeds xxx xxx xxx]: a list of seeds used to generate datasets. the first seed for the first file, the second seed for the second file... total need [number of datasets] seeds. Some seed examples:
3838846863, 
3814455767, 
4294778037, 
1523572004, 
2197961198, 
3953912073, 
385400468, 
4089372537, 
2298599557, 
3878630449

Example:
```
./gen_data 100 2 1 1 10 3 -seeds 3838846863 3814455767 4294778037
```

#### Modify Your Own Datasets
If you want to run the model on your own dataset, please follow the following format

##### Input Format

The raw data file uses the following format:
```
Number of data owners (n)
Dimension of the data points (d)
Data sets from n data owner
```

Then, for each data owner, data set uses the following format:
```
Number of data points (m) for the ith data owner
Data points, each line containing the point ID followed by the d attributes, separated by spaces
```

Example:
```
2                # Specifies 2 data owners
3                # Specifies a 3-dimensional data point
2                # Specifies 2 data points for the first data owner
1 51 35 14       # First data point for data owner 1
2 49 30 14       # Second data point for data owner 1
1                # Specifies 1 data point for the second data owner
3 47 32 13       # Data point for data owner 2
```


## Calculate Shapley Values

To calculate the Shapley values, use the following command:

```
python calculate_sv.py [OPTIONS]
```
### Options

- -i, --input Path to the input data
- -t, Specifies the model type.
    - Candidates: wssv, wssv2d, pssv, pssv2d, enum, enum2d, perm, perm2d, parr, parr2d
- -p Parallelism level. Default is 16.
- -s Number of Permutations. Only used for perm and perm2d models. Default 100.
- --niter Number of iterations. Default is 5.
- --neval Number of evaluations. Default is 100.
- --timeout Set the timeout for the calculation. Default is 3600 second.

Example
```
python calculate_sv.py -i ${PWD}/data/data_for_test/2d10000n_10p_anti_layer_0 -t perm2d -s 100 -p 16 --niter 5 --neval 100 --timeout 3600
```

## Find the Number of Samples, neval

To find the number of samples needed for the approximated Integration, use the following command:

```
python find_sample.py [OPTIONS]
```
### Options

- -i, input Path to the input data
- -o, output directory. Default ./output 
- -t, Specifies the model type.
    - Candidates: wssv, wssv2d, pssv, pssv2d, enum, perm, parr
- -p Parallelism level. Default is 16.
- -s Number of Samples. Only used for perm model. Default 100.
- --niters Number of iterations to be tested. Default is 5.
- --nevals Number of evaluations to be tests. Default is 10 20 50 100 200.
- --times Number of times to run a certain configuration. Default is 10.
- --th Set the threshold for the volatility of SV, the finding process ceases when the calculated volatility is less than the threshold. Default is 0.01.

Example
```
python find_samples.py -t wssv --niters 5 --times 3 --nevals 10 20 50 -p 16 -i ${PWD}/data/data_for_test/2d10000n_10p_anti_layer_0
```

For more examples, please refer to ```./examples.sh```
