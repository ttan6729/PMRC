#!/bin/bash
set -e

#to be added: 
usage()
{
cat << EOF
MRC is a clsutering-bsaed tool for selecting high-simialrity groups from multiple FASTQ datasets, then apply minicom/PgRC for compression.

Usage: 
Compression - compresses FASTQ datasets. Output written to '*.MRC' file
./MRC.sh -a m -r list1.txt (compress with minicom, file contain name of to be compressed files)
./MRC.sh -a p -r file.txt (compress with PgRc, file contain name of to be compressed files)
./MRC.sh -d file.MRC
Options:
	-r      compression mode
	-h 		print help message
	-t 		number of threads, default: 12
	-k 		length of k-mer, k <= 10, default: 8
	-e 		threshold percentage, default: 2
Decompression - decompresses reads. Output written to 'dec' folder
./minicom -d file.minicom 
	-d 		a compressed file .MRC [only for decompression]
 	-t 		number of threads, default: 24
#See README and more supplementary information at:
EOF
# exit 0
}

compress()
{
	output=${filename%.*}_${alg}${mark}_PMRC
	echo "file:${filename} folder:${output}, alg: ${alg}, mark: ${mark}"
	rm -rf $output ${filename%.*}_${alg}${mark}.MRC
	mkdir $output
	echo ${alg} > ${output}/info

	echo "./PMRC -r ${filename} -o ${output}/cluster"
	./PMRC -r ${filename} -o ${output}/cluster #-t ${num_thr} -k ${k} -e ${threshold_per} -o ${output}

	listVar=( )
	while read a b
	do
	        listVar+=($a)
	done < "${filename}"
	while IFS= read -r line
	do
	        fp=""
	        for a in $line
	        do
	                fp="${fp}${a}_"
	        done
  
			rm -rf ${output}/buff.fastq
        	for a in $line
        	do
            	cat ${listVar[a]} >> ${output}/buff.fastq
        	done
        	mv ${output}/buff.fastq ${output}/${fp}.fastq            	

            if [[ $alg = "p" ]]; then
            	./PgRC -o -i ${output}/${fp}.fastq ${output}/${fp}.pgrc
            elif [[ $alg = "m" ]]; then
            	cd minicom 
            	./minicom -r ../${output}/${fp}.fastq -p
            	cd ../
            elif [[ $alg = "s" ]]; then
                        ./spring -c -i ${output}/${fp}.fastq -o ${output}/${fp}.spring
            elif [[ $alg = "s2" ]]; then
                        ./spring -c -i ${output}/${fp}.fastq --no-quality --no-ids -o ${output}/${fp}.spring
            elif [[ $alg = "f" ]]; then
                    echo ${fp}.fastq
                    _cwd="$PWD"
                    cd FaStore
                    sh ./fastore_compress.sh --lossless --in ../${output}/${fp}.fastq --out ../$
                    cd ../
       		fi
			rm -rf ${output}/${fp}.fastq
	done < "${output}/cluster"
	tar -cf ${filename%.*}_${alg}${mark}.PMRC ${output}
	echo "result write to ${filename%.*}_${alg}${mark}.PMRC"
	#rm -rf ${output}
}

decompress()
{	
	dir=${filename%.PMRC*}_PMRC
	output=${filename%.PMRC*}_decompress
	echo "dir is ${dir}"
	rm -rf ${dir}
	echo "decompress, file: ${filename}"
	echo "overwrite -xvkf ${filename}"
	tar --overwrite -xvkf  ${filename}
	echo "finish tar"

	rm -rf $output
	mkdir $output
	
	file_list=( )
	length_list=( )

	alg=$(head -n 1 ${dir}/info)
	echo "alg:${alg}"

	{
	read
    while read a b
	do
	        file_list+=($a)
	        length_list+=($b)
	done
	} < "${dir}/info"

	while IFS= read -r line
	do
		fp=""
		IFS=', ' read -r -a array <<< $line
		for element in "${array[@]}"
		do
    		fp="${fp}${element}_"
		done

        if [[ $alg = "p" ]]; then
        	./PgRC -d ${dir}/${fp}.pgrc
        	start=1
        	end=0
        	for element in "${array[@]}"
			do
				end=$(($end + ${length_list[$element]}))
				sed -n -e "${start},${end} p" -e "${end} q" ${dir}/${fp}.pgrc_out > ${output}/${file_list[$element]}
				#sed -n \'\' ${dir}/${fp}.pgrc_out > ${output}/${file_list[$element]}
				start=$(($end+1))
			done

        elif [[ $alg = "m" ]]; then
        	cd minicom 
        	./minicom -r ../${dir}/${fp}.fastq -p
        	cd ../
        elif [[ $alg = "s" ]]; then
        	fp="${fp}.spring"
        	output=""
        	for element in "${array[@]}"
			do
		    	output="${output} ${file_list[$element]}"
			done
        	./spring -d -i ${fp} -o ${output}
		fi		
	done < "${dir}/cluster"
	#rm -rf $dir

	echo "finished, decompressed file write to ${output}"

}
#Initialize variables to default values.
num_thr=12
threshold_per=2
m_dict=0
k=7
alg=""
filename=""
#Check the number of arguments. If none are passed, print help and exit.
argnum=$#
if [[ $argnum -eq 0 || $1 == "-h" ]]; then
 usage
 exit 1
fi
mark=""
mode="c"
if [[ $1 == "-d" ]]; then
	mode="d"
fi

while getopts ":a:r:t:k:e:m:" opt; do
	case "$opt" in
		a) alg=$OPTARG;;
		r) filename=$OPTARG;;
		t) num_thr=$OPTARG;;
		k) k=$OPTARG;; #length of k
		e) threshold=$OPTARG;; #different threshold
		m) mark=$OPTARG;; #additional suffix for test
		#\?) usage; echo -e "\033[31m Error parameters. \033[0m"; exit 0;;
		#*) usage; echo -e "\033[31m Error parameters. \033[0m"; exit 0;;
	esac
done

if [[ $mode == "c" ]]; then
	compress
elif [[ $mode == "d" ]]; then
	decompress
fi
