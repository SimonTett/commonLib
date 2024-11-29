#!/bin/bash -l

# next two lines only needed for reticulate
# module load scitools
# conda activate sjb-iris-env

shopt -s extglob

# prog expects
# cargs  <- commandArgs(trailingOnly = TRUE)
# flist  <- cargs[1] # "/project/hires_rcm/UKCP18/cpm_rerun_output/TS1/halfhourly/r001i1p00000/pr/r001i1p00000_19801201-19801230_pr.nc"
# stseas <- cargs[2] # 'djf'

f2run="batch.filter_hourly.R"

### HOURLY  ########################################################
stseas='jja'
flist0=/scratch/hadsx/cpm/hourly/filtered/r001i1p00000/pr/XXX/r001i1p00000_????@(06|07|08)01-????@(06|07|08)30_pr.nc

# stseas='djf'
# flist0=/project/hires_rcm/UKCP18/cpm_rerun_output/XXX/hourly/r001i1p00000/pr/r001i1p00000_????@(12|01|02)01-????@(12|01|02)30_pr.nc
# flist0=/scratch/hadsx/cpm/hourly/filtered/r001i1p00000/pr/XXX/r001i1p00000_????@(12|01|02)01-????@(12|01|02)30_pr.nc

TSall='TS1 TS2 TS3 TS4 TS5' #  'TS1' #'TS1 TS2' #

####################################################################
spicecount() {
    squeue -u hadsx |awk '{printf("%s\n",$3)}'|grep wrap|wc -l
}


for TS0 in $TSall; do

    flist=${flist0/XXX/$TS0}

    out1="rout/batch.filter.$TS0.$stseas"
    # echo $out1

    count=0

    for f in $flist; do

        let count=$count+1
        echo $f

        fout=$out1$count'.rout'
        # echo $fout

        # echo $PWD/$f2run $exptid $TS0 $stseas

        f3run=$PWD/$f2run

        nspice="$(spicecount)"    # cnt="$(cmd)"
        while [ $nspice -gt 300 ]; do
            # echo "Sleeping "
            sleep 60
            nspice="$(spicecount)"
            # echo "count is" $count
        done

        sbatch --time=240 --mem=8000 -o $PWD/$fout --wrap="Rscript $f3run $f $stseas"
        # sbatch --time=30 -o $PWD/$fout --wrap="Rscript $f3run $f $stseas"

        # echo $PWD/$fout
        # echo

    done
done

# sbatch --time=15 --mem=8000 -o test_reticulate.rout --wrap="Rscript /home/h03/hadsx/extremes/tawn/jordan/code/cpm/test_reticulate.R"
