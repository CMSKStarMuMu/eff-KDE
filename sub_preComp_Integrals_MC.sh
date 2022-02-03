#!/bin/bash

par=${1}
[ -z "${par}" ] && par=1	# set default
year=${2}
[ -z "${year}" ] && year=2016	# set default

# Compile file
if make preComp_Integrals_MC; then

    tag=1

    while [ ${tag} -ge 0 ]; do
	while read -a line; do

	    bin=${line[0]}
	    vers=${line[7]}

	    # Create directory for log files
	    if [ ! -d logs_preComp_v${vers} ]; then mkdir logs_preComp_v${vers}; fi
	    if [ ! -d tmpint_v${vers} ]; then mkdir tmpint_v${vers}; fi

	    if [ ! -d tmpdir-precom_v${vers} ]; then mkdir tmpdir-precom_v${vers}; fi
	    cd tmpdir-precom_v${vers}

	    echo "Submitting tag ${tag}, bin ${bin}, vers ${vers}"

	    sbatch -a 0-49 \
		   --export=bin=${bin},par=${par},tag=${tag},cnt_hit=100000000,year=${year},vers=${vers} \
		   ../run_preComp_Integrals_MC.sh

	    # sbatch -a 0-499 \
	    # 	   --export=bin=${bin},par=${par},tag=${tag},cnt_hit=10000000,year=${year},vers=${vers} \
	    # 	   ../run_preComp_Integrals_MC.sh

	    cd ..
	    
	done < ../confSF/KDE_SF.list
	tag=$(( ${tag} - 1 ))
    done
fi
