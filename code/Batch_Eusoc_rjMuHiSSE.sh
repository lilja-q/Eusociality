#!/bin/bash

PATH=$PATH:/revbayes/projects/installer:/revbayes/projects/installer/rb

label_LIST=$(seq 0 0)
label=("RIS")



rep_nums=$(seq 0 9)
RUN_NUM=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")


# Initialize a counter for job submissions
counter=0


for t in ${label_LIST[@]}
do	



# Loop through each rep_num
for r in ${rep_nums[@]}; do

    # Construct the job name
    NAME="${RUN_NUM[$r]}_${label[$t]}"
    
      # Submit the job
      bsub -G compute-michael.landis \
      -cwd /storage1/fs1/michael.landis/Active/Eusociality/ \
      -o /storage1/fs1/michael.landis/Active/Eusociality/output_MuHiSSE/stdout/$NAME  \
      -J $NAME \
      -q general \
      -g /m.seanwmchugh/Eusociality \
      -n 1 -M 3GB -R "rusage [mem=3GB] span[hosts=1]" \
      -a 'docker(sswiston/rb_tp:7)' /bin/bash /storage1/fs1/michael.landis/Active/Eusociality/code/Eusoc_rjMuHiSSE.sh

    # Increment the counter
    ((counter++))

    # Check if counter has reached 20
    if [ $counter -eq 20 ]; then
        # Reset counter
        counter=0
        # Wait for 10 seconds
        sleep 10
    fi

done
done