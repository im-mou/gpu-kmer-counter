# gpu-kmer-counter
Speeding up the algorithm to count K-mers in a genetic sequence using GPUs.

## Abstract
K-mer counting is a process with the goal of creating a histogram of all possible combinations of length k for an input string S. From an algorithmic point of view, counting k-mers in a string seems like a very simple task but with recent advances in sequencing technology, more and more sequencing machines are generating a large amount of data in a very short time and makes the simple task of generating a histogram a challenge. In recent years, the performance of k-mer counting algorithms has improved significantly, and there has been much interest in using graphics processing units (GPUs) to accomplish the task of counting k-mers. The fundamental purpose of this research is to analyze different algorithms to count the number of occurrences in a sequence with different k-mer settings and subsequently to optimize and speed up one of the algorithms by using GPUs.

## Source repository
[Source code repository: lh3/kmer-cnt](https://github.com/lh3/kmer-cnt)
<br/>
This repository contains all the the source code of the diferent implementations that has be used and experimented with for this research.

## Instructions to use and test the implementations

```sh
git clone https://github.com/im-mou/gpu-kmer-counter
cd gpu-kmer-counter
make
```

### Download and parse the dataset for different scripts

```sh
wget https://github.com/lh3/kmer-cnt/releases/download/v0.1/M_abscessus_HiSeq_10M.fa.gz
./parse-data ./M_abscessus_HiSeq_10M.fa.gz
```

## Excute implementations

### Secuential: kc-c1-fast.c
```sh
sbatch ./slurm_scripts/slurm-kc-c1-fast.sub
```

### Parallel: cuda-fast.c - Best Implementation
```sh
sbatch ./slurm_scripts/slurm-cuda-fast.sub
```

### Parallel: cuda-dumb.c - Pretty dumb. Non-atomic, experimental purpose only.
```sh
sbatch ./slurm_scripts/slurm-cuda-dumb.sub
```

## Scripts
The two properly working final scripts with the correct outputs are the following:
- kc-c1-fast.c
- cuda-fast.cu
