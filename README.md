# Amerced DTW demonstration application

Amerced DTW (ADTW) is a variant of DTW, "amercing" ("penalizing") warping steps by a fixed penalty.
This is a NN1 classification demonstration application for our paper (submitted)
> Amercing:  An Intuitive, Elegant and Effective Constraint for Dynamic Time Warping

This program can load `ts` files from the UCR archive (see [timeseriesclassification.com](timeseriesclassification.com)).
In the submitted paper, the penalty is find through Leave-One-Out cross validation; see the `results` folder.

## How to compile and run
As a cmake project, most IDE should be able to open/compile/run the application.
If using the command line (tested under Linux), use the following steps:

```sh
mkdir cmake-build-release
cd cmake-build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make ADTWNN1
./ADTWNN1
```

This last line run the application; it will complain that no argument are provided:
```sh
Error: argument parsing
ADTW NN1 classifier.
Usage:
  ./path/exec <ucr folder> <dataset name> <penalty> <threads>
```

The error message tells us what are the expected arguments:
* `<ucr folder>` path to a folder containing the UCR datsets, in the `ts` format.
* `<dataset name>` name of the dataset in the folder. The program will access the corresponding `_TRAIN.ts` and `_TEST.ts` files.
* `<penalty>` The additive penalty, should be tune per dataset
* `<thread>`  Number of thread to use. Must be given, if only using `1` thread.

## Early Abandoned and Pruned Implementation (EAP)
The ADTW distance is efficiently implementated following
> Early abandoning and pruning for elastic distances including dynamic time warping.
> 
> Data Min Knowl Disc (2021). https://doi.org/10.1007/s10618-021-00782-4
> 
> Herrmann, M., Webb, G.I.

