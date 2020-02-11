#!/usr/bin/env bash
# Run CellRanger
# Peter Hickey
# 2019-09-18

# NOTE: This was sample was sequenced in 2 runs. The first containing 5' gene
#       expression (GEX) and hashtag oligos (HTOs) and the second containing
#       VDJ sequencing (VDJ).

# Setup ------------------------------------------------------------------------

module load cellranger/3.1.0
module load bcl2fastq/2.19.1

# Project specific variables ---------------------------------------------------

PROJECT="C057_Cooney"
# NOTE: eval needed for path expansion of '~' character.
eval PROJECT_ROOT="~/SCORE/${PROJECT}"
SEQDIR=${PROJECT_ROOT}/extdata/jamesC_10X_040919

RUNDIR_GEX_HTO=${SEQDIR}/190903_NB500916_0620_AHHW3NBGXC_jamesC_10X5p_CITE
RUNDIR_VDJ=${SEQDIR}/190904_NB500916_0621_AHNYGGAFXY_jamesC_10XVDJ
SAMPLESHEET_GEX_HTO=${PROJECT_ROOT}/data/sample_sheets/SampleSheet_GEX_HTO.csv
SAMPLESHEET_VDJ=${PROJECT_ROOT}/data/sample_sheets/SampleSheet_VDJ.csv
LIBRARIESCSV=${PROJECT_ROOT}/data/sample_sheets/libraries.csv
FEATUREREFCSV=${PROJECT_ROOT}/data/sample_sheets/feature_ref.csv
TRANSCRIPTOME_GEX_HTO="/wehisan/home/allstaff/h/hickey/grpu_mritchie_0/tools/cellranger/refdata-cellranger-GRCh38-3.0.0"
TRANSCRIPTOME_VDJ="/wehisan/home/allstaff/h/hickey/grpu_mritchie_0/tools/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0"

SAMPLES=( HTLV )

# General variables ------------------------------------------------------------

OUTDIR=${PROJECT_ROOT}/output/CellRanger
mkdir -p ${OUTDIR}
CELLRANGERDIR_GEX_HTO=${SEQDIR}/CellRanger_GEX_HTO
mkdir -p ${CELLRANGERDIR_GEX_HTO}
CELLRANGERDIR_VDJ=${SEQDIR}/CellRanger_VDJ
mkdir -p ${CELLRANGERDIR_VDJ}
DATADIR=${PROJECT_ROOT}/data/CellRanger
mkdir -p ${DATADIR}

# Generate FASTQs with cellranger mkfastq  -------------------------------------

# NOTE: Not specifying --id or --output-dir to avoid GH and VDJ outputs
#       clobbering each other.
for SAMPLE in "${SAMPLES[@]}"
do
  echo ${SAMPLE}

  # NOTE: CellRanger writes a bunch of output to $(pwd) so have to move to the
  #       directory where you want the output to go. Sigh.
  cd ${CELLRANGERDIR_GEX_HTO} || exit
  cellranger mkfastq --run=${RUNDIR_GEX_HTO} \
                     --id=${PROJECT}_GEX_HTO \
                     --simple-csv=${SAMPLESHEET_GEX_HTO} \
                     --qc \
                     --jobmode=local \
                     --localcores=30 \
                     --localmem=200 \
                     --output-dir=${CELLRANGERDIR_GEX_HTO}

  # NOTE: CellRanger writes a bunch of output to $(pwd) so have to move to the
  #       directory where you want the output to go. Sigh.
  cd ${CELLRANGERDIR_VDJ} || exit
  cellranger mkfastq --run=${RUNDIR_VDJ} \
                     --id=${PROJECT}_VDJ \
                     --simple-csv=${SAMPLESHEET_VDJ} \
                     --qc \
                     --jobmode=local \
                     --localcores=30 \
                     --localmem=200 \
                     --output-dir=${CELLRANGERDIR_VDJ}
done

# Construct libraries.csv ------------------------------------------------------

# NOTE: This file is only needed for feature barcoding experiments
#       (e.g., CITE-seq).

echo "fastqs,sample,library_type" > ${LIBRARIESCSV}
for SAMPLE in "${SAMPLES[@]}"
do
  echo "${CELLRANGERDIR_GEX_HTO}/HHW3NBGXC,5GEX,Gene Expression" >> ${LIBRARIESCSV}
  echo "${CELLRANGERDIR_GEX_HTO}/HHW3NBGXC,HashTag,Antibody Capture" >> ${LIBRARIESCSV}
done

# Single-sample analyses with cellranger count / vdj ---------------------------

for SAMPLE in "${SAMPLES[@]}"
do
  echo ${SAMPLE}

  # NOTE: CellRanger writes a bunch of output to $(pwd) so have to move to the
  #       directory where you want the output to go. Sigh.
  cd ${CELLRANGERDIR_GEX_HTO} || exit

  # NOTE: Normally would do extract lines, but samples aren't labelled properly.
  # head -n 1 ${LIBRARIESCSV} > ${SAMPLE}.libraries.csv
  # grep "${SAMPLE}" ${LIBRARIESCSV} >> ${SAMPLE}.libraries.csv
  cp ${LIBRARIESCSV} ${SAMPLE}_GEX_HTO.libraries.csv

  cellranger count --id=${SAMPLE}_GEX_HTO \
                   --transcriptome=${TRANSCRIPTOME_GEX_HTO} \
                   --libraries=${SAMPLE}_GEX_HTO.libraries.csv \
                   --feature-ref=${FEATUREREFCSV} \
                   --localcores=30 \
                   --localmem=200

  cp ${SAMPLE}_GEX_HTO/outs/web_summary.html \
     ${OUTDIR}/${SAMPLE}_GEX_HTO.web_summary.html
  cp ${SAMPLE}_GEX_HTO/outs/cloupe.cloupe ${OUTDIR}/${SAMPLE}_GEX_HTO.cloupe

  # NOTE: CellRanger writes a bunch of output to $(pwd) so have to move to the
  #       directory where you want the output to go. Sigh.
  cd ${CELLRANGERDIR_VDJ} || exit
  cellranger vdj --id=${SAMPLE}_VDJ \
                 --reference=${TRANSCRIPTOME_VDJ} \
                 --fastqs=${CELLRANGERDIR_VDJ}/HNYGGAFXY \
                 --sample=VDJ \
                 --localcores=30 \
                 --localmem=200

  cp ${SAMPLE}_VDJ/outs/web_summary.html \
     ${OUTDIR}/${SAMPLE}_VDJ.web_summary.html
  cp ${SAMPLE}_VDJ/outs/vloupe.vloupe ${OUTDIR}/${SAMPLE}_VDJ.vloupe
  cp ${SAMPLE}_VDJ/outs/filtered_contig_annotations.csv \
     ${DATADIR}/${SAMPLE}.filtered_contig_annotations.csv
done
