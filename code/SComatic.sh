# Activate conda environment if needed
# conda activate SComatic

# Create an output folder and go to the main SComatic folder
SCOMATIC=/scratch/devel/pnieto/projects/CSF/data/SComatic/SComatic-main
project=CSF
subproject=CSF_01
sample=4839
output_dir=/scratch/devel/pnieto/projects/$project/data/$subproject/output/05_SComatic/$sample
mkdir -p $output_dir

# Step 1: Splitting alignment file in cell type specific bams
output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1
bam_file=/scratch/devel/pnieto/projects/$project/data/$subproject/jobs/$sample/$sample/outs/per_sample_outs/$sample/count/sample_alignments.bam

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py --bam $bam_file \
        --meta $output_dir/barcode_annotations.tsv \
        --id ${sample} \
        --min_MQ=0 \
        --outdir $output_dir1


# Step 2: Collecting base count information

REF=/scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa

output_dir2=$output_dir/Step2_BaseCellCounts
mkdir -p $output_dir2

for bam in $(ls -d $output_dir1/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 12

  rm -rf $temp
done

# Step 3: Merging base count matrices
output_dir3=$output_dir/Step3_BaseCellCountsMerged
mkdir -p $output_dir3

python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv


# Step 4: Detection of somatic mutations

# Step 4.1
output_dir4=$output_dir/Step4_VariantCalling
mkdir -p $output_dir4

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/${sample}.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/${sample} \
          --ref $REF

# Step 4.2
editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/${sample}.calling.step1.tsv \
          --outfile ${output_dir4}/${sample} \
          --editing $editing \
          --pon $PON

bedtools intersect -header -a ${output_dir4}/${sample}.calling.step2.tsv -b $SCOMATIC/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == "PASS"' > ${output_dir4}/${sample}.calling.step2.pass.tsv

# Computing the number of callable sites per cell type
output_dir5=$output_dir/CellTypeCallableSites
mkdir -p $output_dir5

python $SCOMATIC/scripts/GetCallableSites/GetAllCallableSites.py --infile $output_dir4/${sample}.calling.step1.tsv  \
   --outfile $output_dir5/${sample} \
   --max_cov 150 --min_cell_types 2

# Computing the number of callable sites per cell
STEP4_1=$output_dir4/${sample}.calling.step1.tsv

output_dir6=$output_dir/UniqueCellCallableSites
mkdir -p $output_dir6

for bam in $(ls -d $output_dir1/*bam);do  
    cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
    echo $cell_type
    
    temp=$output_dir6/temp_${cell_type}
    mkdir -p $temp

    python $SCOMATIC/scripts/SitesPerCell/SitesPerCell.py --bam $bam    \
       --infile $output_dir4/${sample}.calling.step1.tsv   \
       --ref $REF \
       --out_folder $output_dir6 --tmp_dir $temp --nprocs 1
    echo
done

# Computing the genotype for each cell at the variant sites
META=$output_dir/barcode_annotations.tsv
STEP4_2_pass=${output_dir4}/${sample}.calling.step2.pass.tsv

output_dir7=$output_dir/SingleCellAlleles
mkdir -p $output_dir7

for bam in $(ls -d $output_dir1/*bam);do  
    cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
    
    temp=$output_dir7/temp_${cell_type}
    mkdir -p $temp

    python $SCOMATIC/scripts/SingleCellGenotype/SingleCellGenotype.py --bam $bam  \
        --infile ${STEP4_2_pass}   \
        --nprocs 1   \
        --meta $META   \
        --outfile ${output_dir7}/${cell_type}.single_cell_genotype.tsv  \
        --tmp_dir $temp  \
        --ref $REF

    rm -rf $temp
done