#!/bin/bash

unset OUTPUT EXPORT


TIMESTAMP=`date +%FT%H.%M.%S`


# create function to handle error messages
# ===============================================
Usage() {
	echo
	echo -e "Usage:\t$0 [arguments]"
	echo -e "\tmandatory arguments:\n" \
          "\t\t-e  (--export)  = Working directory (where all the outputs will be printed) (give full path)\n" \
          "\t\t-o  (--output)  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)\n" 
	echo -e "\toptional arguments:\n " \
          "\t\t-h  (--help)  = See helps regarding the pipeline options. \n" 
echo 
}
          # "\t\t-v  (--verbose)  = set verbosity level [CURRENT \"$VERBOSE\"]\n" 


# ===============================================
# PARSING ARGUMENTS
# ===============================================
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:e:o: --longoptions output:,export:,help -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    exit 42
fi



set -- $options

while [ $# -gt 0 ]

do
    case $1 in
    -x| --extra) 
      EXTRA_CONF="$2" ;
      if [ -f $EXTRA_CONF ]; then
        echo "* LOADING EXTRA CONFIG FILE $EXTRA_CONF";
        . $EXTRA_CONF
      else
        echo "ERROR: invalid EXTRA CONFIG file: $EXTRA_CONF";
        echo "Please check options and try again"; exit 42;
      fi
    esac
    shift
done

# ===============================================
# LOAD ALL OTHER OPTIONS
# ===============================================
set -- $options

while [ $# -gt 0 ]
do
    case $1 in
    -h| --help) Usage; exit 0;;
    -o| --output) OUTPUT="$2" ; shift ;;
    -e| --export) EXPORT="$2" ; shift ;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 42;;
    (*) break;;
    esac
    shift
done

