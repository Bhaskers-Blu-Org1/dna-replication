#!/usr/bin/env bash

total_duration=0
current_iteration=0
for n_particles in 10 50 100 150 200 400 600 800 1000
do
	for d_coef in "0.5" "2.0" "3.5"
	do
		for r_attr in "0.015" "0.025" "0.035"
		do
			for p_bind in 1
			do
				for p_act in 1
				do
					t_start=$(date +%s.%N)
					current_iteration=$((current_iteration + 1))
					out_key="SIM_N${n_particles}_D${d_coef}_R${r_attr}_P${p_bind}"
					mpirun /u/amp/jonas/dnarepl-batch -o /u/amp/simulations/setup/data/oripos_heichinger.csv -c /u/amp/simulations/setup/data/contig_wood50.csv -k $out_key -s /u/amp/simulations/setup/data/G1_constraint_100 -w /u/amp/simulations/runA_ran/out -i 5 -g 0.02 -r 1.33 -x 1.54 -f 50 -n $n_particles -d $d_coef -p $p_bind -b $r_attr -a $p_act >> /u/amp/simulations/runA_ran/error
					iteration_duration=$(echo "$(date +%s.%N) - ${t_start}" | bc)
					total_duration=$(echo "${total_duration} + ${iteration_duration}" | bc)
				done
			done
		done
	done
done
