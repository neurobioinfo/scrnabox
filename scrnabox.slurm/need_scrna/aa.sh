#/bin/bash
unset EXPORT OUTPUT

############
# FUNCTION TO PRINT HELP MESSAGE
############
Usage() {
        echo
        echo -e "Usage:\t$0 [arguments]"
        echo -e "\tmandatory arguments:\n" \
          "\t\t-e  (--export)  = choose export. \n" \
          "\t\t-o  (--output)   = choose output \n" 
	echo -e "\toptional arguments:\n " \
          "\t\t-h  (--help)  = See helps regarding the pipeline options. \n" 
    echo
}


############
# PARSE ARGUMENTS
############
if ! options=$(getopt --name $(basename $0) --alternative --unquoted --options hr:e:o: --longoptions export:,output:,help -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    exit 42
fi

set -- $options
while [ $# -gt 0 ]
do
    case $1 in
    -h| --help) Usage; exit 0;;
    -e| --export) EXPORT="$2"; shift ;;
    -o| --output) OUTPUT="$2"; shift ;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; failure=1;;
    (*) break;;
    esac
    shift
done


my_export=${EXPORT}  
my_export2=($(echo $my_export | tr "," "\n"))

for i in "${my_export2[@]}"
do
    export $i
done

###############################################
echo "-------------------------------------------"
echo "* step2 submitted at `date +%FT%H.%M.%S`"
echo "-------------------------------------------"
echo "* PIPELINE_HOME:        "
echo "* OUTPUT DIR:            "
echo "-------------------------------------------"

Rscript ./test.R 

exit 0
