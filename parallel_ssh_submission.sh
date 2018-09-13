#!/bin/bash

connect_ip () {
	ip=$1
	if [ $ip -eq 181 -o $ip -eq 183 -o $ip -eq 184 -o $ip -eq 197 ]
	then
		return 1 
	else 
		sshhost="10.0.2."$ip
		ssh -o "StrictHostKeyChecking no" bineet.dash@$sshhost "cd ~/Documents/ED/; ./run_run_parallel.sh $2 $3 $4 $5 "
		echo $sshhost done.
	fi	
}

r1=41
r2=$(($r1+1))
r3=$(($r1+2))
r4=$(($r1+3))

for ip in {169..170}
	do
	 connect_ip $ip $r1 $r2 $r3 $r4 &
	 r1=$(($r1+4)); r2=$(($r2+4)); r3=$(($r3+4)); r4=$(($r4+4))
	done

#skip_list=181,183,184,197

# echo $sshhost $2 $3 $4 $5 >> check_oui