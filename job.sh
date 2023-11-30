#!/bin/sh

#module load mpi/intelmpi/2019u3
#module load mathlib/lapack/intel/3.4.2

for((i=1;i<=3;i++)); do
	U=$(awk -v U=$i 'BEGIN{printf "%.1f", (U-1)*0.20}')
	for((j=1;j<=11;j++)); do
		V=$(awk -v V=$j 'BEGIN{printf "%.1f", (V-1)*0.20}')
		for((k=1;k<=1;k++)); do
			# mu=$(awk -v mu=$k 'BEGIN{printf "%.3f", (mu-1)*0.001}')
            # dir=u${U}v${V}m${mu}

			dir=${U}_${V}
            if [ ! -d "$dir" ]; then
				mkdir ${dir}
	        fi
			# mkdir ${dir}
			# cp -r ../code/. ./${dir}

            cd ${dir}
			cp -r ../code/. .
			
			sed -i "s/U=0;V1=0.00;V2=0.00;V3=0;/U=0;V1=$U;V2=$V;V3=0;/g" main.m
			sed -i "s/#SBATCH -J matlab/#SBATCH -J ${U}_${V}/g" matlab.sh

			#   Create a local working directory on scratch
			# mkdir -p $SCRATCH/$SLURM_JOB_ID
			# mkdir -p /public/home/jiangkun/zj/matlabtmp/$SLURM_JOB_ID
			# sbatch matlab.sh
			# Cleanup local working directory
			# rm -rf $SCRATCH/$SLURM_JOB_ID
			# rm -rf /public/home/jiangkun/zj/matlabtmp/$SLURM_JOB_ID

	        # if [ ! -e "flow.dat" ]; then
		    #     # cp ../smfrg_seq.out .
			# 	# cp -r ../code/. ./${dir}
	        #     # echo ${U} > U.input
        	# 	# echo ${V} > Vnn.input
        	# 	# echo ${mu} > mu.input
        	# 	# ./smfrg_seq.out > output
			# fi

			# if [ -e "flow.dat" ]; then
			# 	#echo ${mu} >> ../mu.dat
	    	# 	#sed -n '1,1p' cdwform.dat >> ../f0.dat
	    	# 	tail -1 flow.dat >> ../phd.dat
			# 	#echo $mu
			# 	#tail -1 flow.dat
    		# fi
	
			cd ..
		done
    done
done
