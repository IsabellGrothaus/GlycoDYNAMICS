#!/bin/bash
#!!!IMPORTANT!!!
#This is a bash script that generates a REST-RECT input for GROMACS-PLUMED use, having replica0&1 as an unbiased replica at ground temperature and replica2 at a higher temperature. 

# !!! IMPORTANT!!! Uses an updated local version (July 2021) of partial_tempering that can also handle the GLYCAM06 FF with its different "pairs" section! Modified by Giovanni Bussi directly. 

ml chem/GROMACS

input=$1 #(Your structure file in .gro format)
nrep=$2 #(Number of replicas)
# "effective" temperature range
tmin=$3 #(Temp. of ground replica)
tmax=$4 #(Temp. of highest replica)
# build geometric progression
mdpfile=$5 #(mdp file like in classical MD at ground temperature)
index=$6

list=$(
awk -v n=$nrep \
        -v tmin=$tmin \
        -v tmax=$tmax \
        'BEGIN{for(i=0;i<n;i++){
        t=tmin*exp(i*log(tmax/tmin)/(n-1));
        printf(t); if(i<n-1)printf(",");
                }
        }'
       )

echo 'temp-list', $list

for((i=0;i<$nrep;i++))
do
  # choose lambda as T[0]/T[i]
  # remember that high temperature is equivalent to low lambda
  lambda=$(echo $list | awk 'BEGIN{FS=",";}{print $1/$'$((i))';}')
echo 'lambda', $lambda
  # process topology
  # (if you are curious, try "diff topol0.top topol1.top" to see the changes)
  partial_tempering_pairs.sh $lambda < processed_m.top > topol$i.top
  # prepare tpr file
  # -maxwarn is often needed because box could be charged
  gmx grompp -maxwarn 2 -o topol$i.tpr -f $mdpfile -p topol$i.top -n $index -c $input
done

for((i=0;i<$nrep;i++))
do
mkdir TOPO$i
mv topol$i.top TOPO$i
mv topol$i.tpr TOPO$i
cp plumed_allangles.dat TOPO$i
cp reference.pdb TOPO$i
cd TOPO$i
mv topol$i.tpr topol.tpr
cd ..
done                                 
