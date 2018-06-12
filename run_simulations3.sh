#! /bin/bash
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#////////////// Only block of code that should need to be edited //////////////
repmax=4              #..number of replicas
tmin=310              #..effective temperature minimum
tmax=430              #..430..#..effective temperature maximum
repstep=1             #..starting step number
repstepmax=2          #..ending step number
run_length=4          #..it's always 2000
spring_constant=3.8   #..spring constant for constraints.tcl
batch=1               #..batch number
run=1                 #..run number
CUDA=no               #..CUDA? "yes" or "no"
number_cpus=1         #..number of processors to use per replica
cold_lipid=yes        #..scaling lipids? "yes" or "no"
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////
repmax_1=$(($repmax - 1))    #..number of replicas - 1
repstep_1=$(($repstep + 1))  #..starting step number + 1
R=0.001986                   #..kcal/molK
num_atoms=$(grep '^ATOM' ionized.pdb | wc -l)
RT_for_tmin=$(echo "$R*$tmin" | bc -l | awk '{printf"%2.3f",$0}')
beta_for_tmin=$(echo "1/($R*$tmin)" | bc -l | awk '{printf"%2.5f",$0}')

#edit exchange.f90, compile, and set up the directory
cp random_template.f90 random.f90
a=$(<random.f90);  b=${a//numberreplicas/$repmax};
echo "$b" > random.f90
gfortran -static-libgfortran -o random.exe random.f90

cp exchange_rest2_npt_simple_template.f90 exchange_rest2_npt_simple.f90
a=$(<exchange_rest2_npt_simple.f90); b=${a//numberatoms/$num_atoms}
c=${b//numberreplicas/$repmax};  d=${c//mtemp_RT/$RT_for_tmin}
e=${d//min_temp/$tmin};          f=${e//mintemp_beta/$beta_for_tmin}
echo "$f" > exchange_rest2_npt_simple.f90
gfortran -static-libgfortran -o exchange_rest2_npt_simple.exe exchange_rest2_npt_simple.f90

for i in $(seq 1 $repmax); {
    cp mem_equil.coor abf_quench${i}_0.coor;
    cp mem_equil.vel abf_quench${i}_0.vel;
    cp mem_equil.xsc abf_quench${i}_0.xsc;
}

cp abf_quench_template.namd abf_quench_template_.namd
cp abf_quench_template_simple.namd abf_quench_template_simple_.namd

#..preprocessing files (happens once)..////////////////////////////////////////
while [ $repstep -lt $repstep_1 ];do
    list=$(awk -v n=$repmax -v tmin=$tmin -v tmax=$tmax \
    'BEGIN{for(i=0;i<n;i++){
        t=tmin*exp(i*log(tmax/tmin)/(n-1));printf(t);if(i<n-1)printf(" ");}}')
		
    #...obtaining prefactors for scaling /////////////////////////////////////
    p=0
    while [ $p -lt $repmax_1 ];do
        for a in $(echo $list);{
            beta[p]=$(echo "(1/($R*$a))" | bc -l | awk '{printf"%2.5f",$0}')
            test1[p]=$(echo "(sqrt(${beta[0]}/${beta[$p]}))" | bc -l | awk '{printf"%2.5f",$0}')
            test2[p]=$(echo "(${beta[0]}/${beta[$p]})" | bc -l | awk '{printf"%2.5f",$0}')
            p=$(($p+1))
        }
    done
	
    # scale psf ///////////////////////////////////////////////////////
    o=1
    while [ $o -le $repmax ];do
        for a in ${test1[@]};do
            head -21 ionized.psf > temp1
            #..lipid
            awk '{if ($2=="MEMB") print $0;}' ionized.psf > temp2
            if [ $cold_lipid == yes ];then
                awk '{printf "%.6f\n", $7*a}' a=$a temp2 > temp2_1
                awk 'FNR==NR{a[NR]=$0;next} {sub($7, a[FNR])}1' temp2_1 temp2 > temp2_2
            fi
            #..water
            awk '{if ($2=="TIP3") print $0;}' ionized.psf > temp3
            awk '{printf "%.6f\n", $7*a}' a=$a temp3 > temp3_1
            awk 'FNR==NR{a[NR]=$0;next} {sub($7, a[FNR])}1' temp3_1 temp3 > temp3_2
            #..water
            awk '{if ($2=="WT1") print $0;}' ionized.psf > temp4
            awk '{printf "%.6f\n", $7*a}' a=$a temp4 > temp4_1
            awk 'FNR==NR{a[NR]=$0;next} {sub($7, a[FNR])}1' temp4_1 temp4 > temp4_2
            #..water
            awk '{if ($2=="WT2") print $0;}' ionized.psf > temp5
            awk '{printf "%.6f\n", $7*a}' a=$a temp5 > temp5_1
            awk 'FNR==NR{a[NR]=$0;next} {sub($7, a[FNR])}1' temp5_1 temp5 > temp5_2
            #..peptides and ions
            awk '{if ($2=="P1") print $0;}' ionized.psf > temp6
            awk '{if ($2=="R2") print $0;}' ionized.psf > temp7
            awk '{if ($2=="ION") print $0;}' ionized.psf > temp8
            sed -n '27061,65237p' ionized.psf > temp9
            if [ $cold_lipid == yes ];then
                cat temp1 temp2_2 temp3_2 temp4_2 temp5_2 temp6 temp7 temp8 temp9 > ionized$o.psf
            elif [ $cold_lipid == no ];then
                cat temp1 temp2 temp3_2 temp4_2 temp5_2 temp6 temp7 temp8 temp9 > ionized$o.psf
            fi
            rm temp*
            o=$(($o+1))
        done
    done
	
    #scale prm for water ////////////////////////////////////////////////////
    for ((i=1;i<=$repmax;i++)) do
        for y in ${test2[@]};do
            scaled_OT2=$(echo "$y*-0.1521" | bc -l | awk '{printf "%f", $0}')
            scaled_HT2=$(echo "$y*-0.046" | bc -l | awk '{printf "%f", $0}')
            cp par_all36_prot.prm par_all36_prot$y.prm
            sed '3461s/-0.1521/'$scaled_OT2'/g' par_all36_prot$y.prm > temp1
            sed '3460s/-0.046/'$scaled_HT2'/g' temp1 > temp2
            mv temp2 par_all36_prot$i.prm
            rm temp1 par_all36_prot$y.prm
			i=$(($i+1))
        done
    done
	
    #scale lipid prm ////////////////////////////////////////////////////////
    for ((i=1;i<=$repmax;i++)) do
        for y in ${test2[@]};do
            if [ $cold_lipid == yes ];then
                head -58 par_all36_lipid.prm > temp1
                sed -n '59,108p' par_all36_lipid.prm > temp2
                awk '{printf "%.2f\n", $3*y}' y=$y temp2 > temp2_1
                awk 'FNR==NR{a[NR]=$0;next} {sub($3, a[FNR])}1' temp2_1 temp2 > temp2_2
                sed -n '109,123p' par_all36_lipid.prm > temp3  ### okay
				
                sed -n '124,254p' par_all36_lipid.prm > temp4
                awk '{printf "%.2f\n", $4*y}' y=$y temp4 > temp4_1
                awk 'FNR==NR{a[NR]=$0;next} {sub($4, a[FNR])}1' temp4_1 temp4 > temp4_2
                awk '{if ($6!="!"){printf "%.2f\n", $6*y}else{print "!"}}' y=$y temp4_2 > temp4_3
                awk 'FNR==NR{a[NR]=$0;next} {sub($6, a[FNR])}1' temp4_3 temp4_2 > temp4_4  
                sed -n '255,265p' par_all36_lipid.prm > temp5  ### okay
				
                sed -n '266,458p' par_all36_lipid.prm > temp6
                awk '{if ($5!="!"){printf "%.4f\n", $5*y}else{print "!"}}' y=$y temp6 > temp6_1
                awk 'FNR==NR{a[NR]=$0;next} {sub($5, a[FNR])}1' temp6_1 temp6 > temp6_2 
                sed -n '459,470p' par_all36_lipid.prm > temp7  ### okay
				
                sed -n '471,474p' par_all36_lipid.prm > temp8
                awk '{if ($5!="!"){printf "%.2f\n", $5*y}else{print "!"}}' y=$y temp8 > temp8_1
                awk 'FNR==NR{a[NR]=$0;next} {sub($5, a[FNR])}1' temp8_1 temp8 > temp8_2
                sed -n '475,486p' par_all36_lipid.prm > temp9  ### okay
				
                sed -n '487,521p' par_all36_lipid.prm > temp10
				awk '{if ($3!="!"){printf "%.4f\n", $3*y}else{print "!"}}' y=$y temp10 > temp10_1
                awk 'FNR==NR{a[NR]=$0;next} {sub($3, a[FNR])}1' temp10_1 temp10 > temp10_2
				awk '{if ($6!="!"){printf "%.2f\n", $6*y}else{print "!"}}' y=$y temp10_2 > temp10_3
                awk 'FNR==NR{a[NR]=$0;next} {sub($6, a[FNR])}1' temp10_3 temp10_2 > temp10_4
                sed -n '522,533p' par_all36_lipid.prm > temp11
                cat temp1 temp2_2 temp3 temp4_4 temp5 temp6_2 temp7 temp8_2 temp9 temp10_4 temp11 > par_all36_lipid$i.prm
                rm temp*
            elif [ $cold_lipid == no ];then
                cp par_all36_lipid.prm par_all36_lipid$i.prm
            fi
			i=$(($i+1))
        done
    done
	
    # scale water and ion ////////////////////////////////////////////////////
    for ((i=1;i<=$repmax;i++)) do
        for y in ${test2[@]};do
            scaled_k=$(echo "$y*$spring_constant" | bc -l | awk '{printf "%f", $0}')
            cp constraints_template.tcl constraints$i.tcl
            a=$(<constraints$i.tcl)
            b=${a//insert_k_here/$scaled_k}
            echo "$b" > constraints$i.tcl
            i=$(($i+1))
        done
    done
    printf "%s\n" ${list[@]} > RepTemp.dat

    if [ -r RunStatus${batch} ];then
        rm RunStatus${batch}
    fi
    rm abf_quench*.out xsc*.dat #abf_quench_template_*.namd
	
    #..preparing config files for replicas..///////////////////////////////////
	for ((rep=1;rep<=repmax;rep++)) do
		reptemp=`head -$rep RepTemp.dat | tail -1`
	  	cp abf_quench_template_.namd abf_quench_template_${rep}.namd
		sed 's/FinalFile/abf_quench'$rep'/g' abf_quench_template_${rep}.namd > temp1
	  	sed 's/RepTemp/'$reptemp'/g' temp1 > temp2
	  	sed 's/RunLength/'$run_length'/g' temp2 > temp3
		sed 's/RestartFile/abf_quench'$rep'_i/g' temp3 > temp4
	  	sed 's/OutputFrequency/'$run_length'/g' temp4 > temp5
	  	sed 's/Param1/par_all36_prot'$rep'.prm/g' temp5 > temp6
	  	sed 's/Param2/par_all36_lipid_edit'$rep'.prm/g' temp6 > temp7
	  	sed 's/PSF/ionized'$rep'.psf/g' temp7 > temp8
		sed 's/constraints.tcl/constraints'$rep'.tcl/g' temp8 > temp9
	  	mv temp9 abf_quench_template_${rep}.namd
		rm temp*
		if [ $repstep -eq 1 ];then
	    	rm abf_quench${rep}.* abf_quench${rep}_i.*
	 	fi
	done
	
    #..Start of simulations ///////////////////////////////////////////////
    for ((repstep=repstep;repstep<=repstepmax;repstep++)) do
        echo 'Replica step: '$repstep' ' >> RunStatus${batch}
        for ((rep=1;rep<=repmax;rep++)) do
            reptemp=$(head -$rep RepTemp.dat | tail -1)
            echo $repstep $rep $reptemp
            echo $repstep $run > status.dat
            if [ $repstep -lt 10 ]; then
                repstepi=000$repstep
            fi
            if [ $repstep -ge 10 ]; then
                if [ $repstep -lt 100 ]; then
                    repstepi=00$repstep
                fi
            fi
            if [ $repstep -ge 100 ]; then
                if [ $repstep -lt 1000 ]; then
                    repstepi=0$repstep
                fi
            fi
            if [ $repstep -ge 1000 ]; then
                repstepi=$repstep
            fi
            if [ $repstep -ne 1 ]; then
                mv abf_quench${rep}.coor abf_quench${rep}_i.coor
                mv abf_quench${rep}.vel abf_quench${rep}_i.vel
                mv abf_quench${rep}.xsc abf_quench${rep}_i.xsc
            else
                cp abf_quench${rep}_0.coor abf_quench${rep}_i.coor
                cp abf_quench${rep}_0.vel abf_quench${rep}_i.vel
                cp abf_quench${rep}_0.xsc abf_quench${rep}_i.xsc
            fi
            tail -1 abf_quench${rep}_i.xsc > temp1000
            read dum x dum dum dum y dum dum dum z dum dum dum < temp1000
            rm temp1000
            sed 's/RandSeed/201'${run}${repstep}${rep}'/g' abf_quench_template_${rep}.namd > temp1
            sed 's/DCDFile/abf_quench'${run}'_'${repstepi}'_'${rep}'/g' temp1 > temp2
            mv temp2 abf_quench${rep}.namd
            rm temp*
			
            ready=0
            while [ $ready -eq 0 ]; do
                if [ -r abf_quench${rep}.namd ]; then
                    nlines=$(grep 'wrapAll' abf_quench${rep}.namd | wc -l)
                    if [ $nlines -eq 1 ]; then
                        ready=1
                    fi
                fi
                sleep 1
            done
			
            ######################## if using CUDA compatible namd then call it namd3
            if [ $CUDA == no ];then
                ./namd2 +p$number_cpus abf_quench${rep}.namd > abf_quench${rep}_${repstep}.out &
            elif [ $CUDA == yes ];then
                ./namd3 +idlepoll +p$number_cpus abf_quench${rep}.namd > abf_quench${rep}_${repstep}.out &
            fi
            echo '   sending replica '$rep' '$reptemp' ' >>  RunStatus${batch}
        done
		
		complete=0
		while [ $complete -ne $repmax ]; do  #...checking if sims are done
			sleep 1
			complete=0
			for ((rep=1;rep<=repmax;rep++)) do
				if [ -r abf_quench${rep}.coor ]; then
					nlines=$(grep '^ATOM' abf_quench${rep}.coor | wc -l)
					if [ $nlines -eq $num_atoms ]; then
						complete=$(($complete + 1))
					fi
				fi
			done
		done
		echo 'Replica step: '$repstep' completed' >>  RunStatus${batch}
		#cp abf_quench1_${repstep}.out normal1_${repstep}.out ########### for saving normal.out
		#.. get potential energy from simulations
		
		for ((rep=1;rep<=repmax;rep++)) do
			grep '^ENERGY' abf_quench${rep}_${repstep}.out | cut -c 192-205 > potential_${rep}.dat
			potential_normal[rep]=$(tail -1 potential_${rep}.dat | cut -c1-20)
		done
		
		rm RepBox_normal.dat
		for ((rep=1;rep<=repmax;rep++)) do
			tail -1 abf_quench${rep}.xsc > temp1
			read dum x dum dum dum y dum dum dum z dum dum dum < temp1
			echo $x $y $z >> RepBox_normal.dat
			cat temp1 >> xsc_normal${rep}.dat
			rm temp*
		done
		
		#.. all simulatons are completed, now namd_energy
		./random.exe < status.dat > namd_energy_out
		#namd_energy_out1=$(cat namd_energy_out)
		#echo "$namd_energy_out1" >> is_it_random.dat
		npairs=$(wc -l < namd_energy_out)
		for ((np=1;np<=npairs;np++)) do
			echo $repstep $run $np > status1.dat
			head -$np namd_energy_out|tail -1 > namd_energy_out_temp
			read rep1 rep2 < namd_energy_out_temp
			
			rep1temp=($(awk -v a=$rep1 'NR==a' RepTemp.dat))
			rep2temp=($(awk -v a=$rep2 'NR==a' RepTemp.dat))
              
		  	cp abf_quench_template_simple_.namd abf_quench_simple_${rep1}.namd
		  	sed 's/RestartFile/abf_quench'$rep2'_i/g' abf_quench_simple_${rep1}.namd > temp1  ######RestartFile/abf_quench'$rep2'
		  	sed 's/OutputFrequency/'$run_length'/g' temp1 > temp2
		  	sed 's/Param1/par_all36_prot'$rep1'.prm/g' temp2 > temp3
		  	sed 's/Param2/par_all36_lipid_edit'$rep1'.prm/g' temp3 > temp4
		  	sed 's/PSF/ionized'$rep1'.psf/g' temp4 > temp5
			sed 's/DCDFile/abf_quench'${run}'_'${repstepi}'_'${rep2}'/g' temp5 > temp6    ######'${rep2}'/g'
			sed 's/constraints.tcl/constraints'$rep1'.tcl/g' temp6 > temp7
			mv temp7 abf_quench_simple_${rep1}.namd
			rm temp*
			
		  	cp abf_quench_template_simple_.namd abf_quench_simple_${rep2}.namd
		  	sed 's/RestartFile/abf_quench'$rep1'_i/g' abf_quench_simple_${rep2}.namd > temp1  ######RestartFile/abf_quench'$rep1'
		  	sed 's/OutputFrequency/'$run_length'/g' temp1 > temp2
		  	sed 's/Param1/par_all36_prot'$rep2'.prm/g' temp2 > temp3
		  	sed 's/Param2/par_all36_lipid_edit'$rep2'.prm/g' temp3 > temp4
		  	sed 's/PSF/ionized'$rep2'.psf/g' temp4 > temp5
			sed 's/DCDFile/abf_quench'${run}'_'${repstepi}'_'${rep1}'/g' temp5 > temp6    ######'${rep1}'/g'
			sed 's/constraints.tcl/constraints'$rep2'.tcl/g' temp6 > temp7
		  	mv temp7 abf_quench_simple_${rep2}.namd
			rm temp*

			./namd2 +p$number_cpus abf_quench_simple_${rep1}.namd > abf_quench_simple_${rep1}_${repstep}.out &
			./namd2 +p$number_cpus abf_quench_simple_${rep2}.namd > abf_quench_simple_${rep2}_${repstep}.out &
			
	        complete=0
	        while [ $complete -lt 1 ]; do
	            sleep 1
	            complete=0
	            for ((rep=1;rep<=repmax;rep++)) do  #...checking to see if done
	                if [ -r abf_quench_simple_${rep1}_${repstep}.out ] && [ -r abf_quench_simple_${rep2}_${repstep}.out ]; then
	                    nlines=$(cat abf_quench_simple_${rep1}_${repstep}.out abf_quench_simple_${rep2}_${repstep}.out | grep '^Program'|wc -l)
	                    if [ $nlines -eq 2 ]; then
	                        complete=$(($complete + 1))
	                    fi
	                fi
	            done
	        done
			
	        echo 'Replica step: '$repstep' completed rep1 rep2' >>  RunStatus${batch}
			#cp abf_quench_simple_1_${repstep}.out simple_1_${repstep}.out ########### for saving simple.out
			rm RepBox.dat RepEnergy.dat RepTemp1.dat
            grep '^ENERGY' abf_quench_simple_${rep1}_${repstep}.out | cut -c 192-205 > potential_simple_${rep1}.dat
			grep '^ENERGY' abf_quench_simple_${rep2}_${repstep}.out | cut -c 192-205 > potential_simple_${rep2}.dat
			tail -1 abf_quench${rep1}.xsc > temp1
            read dum x dum dum dum y dum dum dum z dum dum dum < temp1
            echo $x $y $z >> RepBox.dat
            cat temp1 >> xsc${rep1}.dat
			tail -1 abf_quench${rep2}.xsc > temp2
            read dum x dum dum dum y dum dum dum z dum dum dum < temp2
            echo $x $y $z >> RepBox.dat
            cat temp2 >> xsc${rep2}.dat
            rm temp*
			
			for ((p=1;p<=repmax;p++)) do
				for a in $(echo $list); do
					beta[p]=$(echo "(1/($R*$a))" | bc -l | awk '{printf "%2.5f", $0}')
					p=$(($p+1))
				done
	        done
			
			potential1=$(tail -1 potential_simple_${rep1}.dat | cut -c1-20)
			potential2=$(tail -1 potential_simple_${rep2}.dat | cut -c1-20)
			echo "${beta[$rep1]} ${potential_normal[$rep1]} $potential1 $rep1" >> RepEnergy.dat
			echo "${beta[$rep2]} ${potential_normal[$rep2]} $potential2 $rep2" >> RepEnergy.dat
			awk -v a=$rep1 'NR==a' RepTemp.dat >> RepTemp1.dat
			awk -v a=$rep2 'NR==a' RepTemp.dat >> RepTemp1.dat
		
	        ./exchange_rest2_npt_simple.exe < status1.dat > exchange_out
            read dum1 dum2 iexchange < exchange_out
            if [ $iexchange -eq 1 ]; then
                mv abf_quench${rep1}.coor restart_coor1
                mv abf_quench${rep2}.coor restart_coor2
                mv restart_coor1 abf_quench${rep2}.coor
                mv restart_coor2 abf_quench${rep1}.coor
                mv abf_quench1s.vel abf_quench${rep1}.vel
                mv abf_quench2s.vel abf_quench${rep2}.vel
                mv abf_quench${rep1}.xsc restart_xsc1
                mv abf_quench${rep2}.xsc restart_xsc2
                mv restart_xsc1 abf_quench${rep2}.xsc
                mv restart_xsc2 abf_quench${rep1}.xsc
            fi
            echo $np $rep1 $rep2 'exchange:' $iexchange >>  RunStatus${batch}
		done
        if [ $repstep -lt $repstepmax ]; then
            rm abf_quench*.out
        fi
	done
done
exit
