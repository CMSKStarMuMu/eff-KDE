#!/bin/bash

par=1
tag=1
year=${1}

# Create directory for log files
if [ ! -d logs_preComp ]; then mkdir logs_preComp; fi

# Compile file
make preComp_Integrals_MC_Swave

while [ ${tag} -ge 0 ]; do
    while read -a line; do

	bin=${line[0]}

	cat << EOF > temp_sub_preComp_Integrals_MC_Swave.sub
Executable  = run_preComp_Integrals_MC_Swave.sh
bin         = ${bin}
par	    = ${par}
tag         = ${tag}
cnt_hit     = 100000000
seed        = \$(ProcId) + 51
year        = ${year}
Arguments   = \$INT(bin) \$INT(par) \$INT(tag) \$INT(cnt_hit) \$INT(seed) \$INT(year)
Log         = logs_preComp/sub_\$(ClusterId).log
Output      = logs_preComp/preComp_Integrals_MC_Swave_\$INT(bin)_\$INT(par)_\$INT(tag)_\$INT(cnt_hit)_\$INT(seed)_\$INT(year).out
Error       = logs_preComp/preComp_Integrals_MC_Swave_\$INT(bin)_\$INT(par)_\$INT(tag)_\$INT(cnt_hit)_\$INT(seed)_\$INT(year).err
+JobFlavour = "tomorrow"
EOF

        if [ "${USER}" == "fiorendi" ]; then
            echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_preComp_Integrals_MC_Swave.sub
        fi

        echo 'Queue 50'>>temp_sub_preComp_Integrals_MC_Swave.sub

        echo "Submit par: "${par}" tag:"${tag}" bin:"${bin} " year:"${year}
	condor_submit temp_sub_preComp_Integrals_MC_Swave.sub
	rm temp_sub_preComp_Integrals_MC_Swave.sub
	    
    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
