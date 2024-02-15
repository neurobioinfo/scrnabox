#!/bin/bash

function call_parameter () {
   a=$1
   my_string=${a#*=}
   my_array=($(echo $my_string | tr "," "\n"))
   for i in "${my_array[@]}"
   do
      eval  $i
   done
}

function remove_argument () {
    X=(ACCOUNT MODULEUSE= THREADS MEM_ARRAY WALLTIME_ARRAY CELLRANGER R_VERSIO step integrate)
    for item in ${X[@]}; do
        sed -i $1 -e "/${item}/d"
    done
}

