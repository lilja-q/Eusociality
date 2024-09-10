PATH=$PATH:/revbayes/projects/installer:/revbayes/projects/installer/rb


#export R_HOME="$R_HOME:/opt/R/4.2.3/lib/R" 
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/R/4.2.3/lib"
#export HOME=“/julia”
#export JULIA_DEPOT_PATH="/storage1/fs1/michael.landis/Active/RFBS_RIS/julia_packages"
#export JULIA_CPU_THREADS=1
#export JULIA_DEPOT_PATH="/root/.julia"
#export JULIA_DEPOT_PATH="/julia/.julia"


IFS="__"
read -ra arr <<< $LSB_JOBNAME
IFS=" "

RUN_NUM=${arr[0]}


rb_command="RUN_NUM=\"$RUN_NUM\";source(\"/storage1/fs1/michael.landis/Active/Eusociality/code/Eusoc_rjMuHiSSE.Rev\");"

echo $rb_command | rb

