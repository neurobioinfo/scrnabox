#!/bin/bash

source $3
OUTFILE=$1/$2
if [ -f $OUTFILE ]; then
    rm $OUTFILE
fi

cat <<EOF > $OUTFILE
#!/bin/bash
unset RUN_NAME LIBRARY EXPECT_CELLS SLURM_TEMPLATE failure REF_DIR

Usage() {
        echo
        echo -e "Usage:\t\$0 [arguments]"
        echo -e "\tmandatory arguments:\n" \\
          "\t\t-r  (--run_name)         run name, serves as prefix for output files\n"
        echo -e "\toptional arguments:\n" \\
          "\t\t-l  (--library_csv)      /path/to/library.csv file       (default: ./library.csv)\n" \\
          "\t\t-c  (--expect_cells)     expected # of cells             (default: 6000)\n" \\
          "\t\t-t  (--slurm_template)   /path/to/slurm.template file    (default: ./slurm.template)\n" \\
          "\t\t-h  (--help)             run this help message and exit\n"
        echo
}

if ! options=\$(getopt --name \$(basename \$0) --alternative --unquoted --options hr:l:f:c:t: --longoptions run_name:,library_csv:,expect_cells:,slurm_template:,help -- "\$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    exit 42
fi

set -- \$options
while [ \$# -gt 0 ]
do
    case \$1 in
    -h| --help) Usage; exit 0;;
    -r| --run_name) RUN_NAME="\$2"; shift ;;
    -l| --library_csv) LIBRARY="\$2"; shift ;;
    -c| --expect_cells) EXPECT_CELLS="\$2"; shift ;;
    -t| --slurm_template) SLURM_TEMPLATE="\$2"; shift ;;
    (--) shift; break;;
    (-*) echo "\$0: error - unrecognized option \$1" 1>&2; failure=1;;
    (*) break;;
    esac
    shift
done

[[ -z \${RUN_NAME} ]] && echo "error: run name must be specified with -r or --run_name" && failure=1
[[ \$failure -eq 1 ]] && echo "ERRORS FOUND. Exiting" && exit 42

source \$OUTPUT_DIR/job_info/.tmp/step1_par.txt
source \$OUTPUT_DIR/job_info/.tmp/temp_config.ini

LIBRARY=\${LIBRARY:-./library.csv}
SLURM_TEMPLATE=\${SLURM_TEMPLATE:-./slurm.template}
# EXPECT_CELLS=${par_expect_cells:-6000}

# module use \$MODULEUSE
if [[ \$MODULEUSE ]]; then module use \$MODULEUSE ; fi
# module load \$MODULECELLRANGER
module load \${CELLRANGER}/\${CELLRANGER_VERSION}

pwd_dir=\$(pwd)
TEMPLOG=\$OUTPUT_DIR/job_info/logs/step_1_\$(basename \${pwd_dir}).log
echo "CELL RANGER is currently running on \$(basename \${pwd_dir}). Please leave it undisturbed until it finishes. "

EOF

cat <<EOF >> $OUTFILE
cellranger count \\
    --id=\${RUN_NAME} \\
    --libraries=\${LIBRARY} \\
    --jobmode=\${SLURM_TEMPLATE} \\
EOF

# par_ref_dir_grch
if [[ -n "${par_ref_dir_grch}" ]]; then
cat <<EOF >> $OUTFILE
    --transcriptome=\${par_ref_dir_grch} \\
EOF
fi

# par_include_introns=TRUE
if [[ -n "${par_include_introns}" ]]; then
if  [[  ${par_include_introns}  =~  yes  ]]  ; then
cat <<EOF >> $OUTFILE
    --include-introns \\
EOF
fi
fi

# par_mempercode=20
if [[ -n "${par_mempercode}" ]]; then
cat <<EOF >> $OUTFILE
    --mempercore=\${par_mempercode} \\
EOF
fi

# par_r1_length=20
if [[ -n "${par_r1_length}" ]]; then
cat <<EOF >> $OUTFILE
    --r1-length ${par_r1_length} \\
EOF
fi

# par_r2_length=20
if [[ -n "${par_r2_length}" ]]; then
cat <<EOF >> $OUTFILE
    --r2-length ${par_r2_length} \\
EOF
fi

# par_no_target_umi_filter
if [[ -n "${par_no_target_umi_filter}" ]]; then
if  [[  ${par_no_target_umi_filter}  =~  yes  ]]  ; then
cat <<EOF >> $OUTFILE
    --no-target-umi-filter \\
EOF
fi
fi

# par_expect_cells
if [[ -n "${par_expect_cells}" ]]; then
cat <<EOF >> $OUTFILE
    --expect-cells ${par_expect_cells} \\
EOF
fi

# par_force_cells
if [[ -n "${par_force_cells}" ]]; then
cat <<EOF >> $OUTFILE
    --force-cells ${par_force_cells} \\
EOF
fi

# par_no_bam
if [[ -n "${par_no_bam}" ]]; then
if  [[  ${par_no_bam}  =~  yes  ]]  ; then
cat <<EOF >> $OUTFILE
    --no-bam \\
EOF
fi
fi

# par_no_libraries
if [[ -n "${par_no_libraries}" ]]; then
if  [[  ${par_no_libraries}  =~  yes  ]]  ; then
cat <<EOF >> $OUTFILE
    --no-libraries \\
EOF
fi
fi

########## END 
cat <<EOF >> $OUTFILE
    2>&1|tee -a \${RUN_NAME}.\$(date +%Y%m%d_%H%M).log   >  \${TEMPLOG}
EOF

cat <<EOF >> $OUTFILE
echo -e "The computation on  \$(basename \${pwd_dir}) is done. "
awk '/Pipestance/'  \${TEMPLOG}
EOF
