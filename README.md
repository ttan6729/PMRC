# PMRC: Path graph based Multiple reads set Clustering and Compression

PMRC is a novel clustering approach for the compression of multiple FASTQ files which construct features space based on the freuqency of minimizers in each file and divided the files into subgroups based on path graph.

## Download & Install

	git clone git@github.com:ttan6729/PMRC.git
	cd MRC
	sh install.sh

## Usage
Usage:
Compression - compresses FASTQ datasets. Output written to '*.PMRC' file
```
./PMRC.sh -a m -r test.txt (compress with minicom, file contain name of to be compressed files)
./PMRC.sh -a p -r test.txt (compress with PgRc, file contain name of to be compressed files)
```
Decompression - decompress compressed files, decompressed fastq files are written to '.PMRC' folder
```
./PMRC.sh -d -r file.PMRC
```
Options:
        -r      compression mode
        -h              print help message
        -t              number of threads, default: 12
        -k              length of k-mer, k <= 10, default: 8
        -e              threshold percentage, default: 2.0
        -d              a compressed file .MRC [only for decompression]
        -t              number of threads, default: 12
        
##Example
```
./PMRC.sh -a p -r test.txt
./PMRc.sh -d test_p.PMRC
```
