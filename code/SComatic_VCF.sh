# Prepare SComatic input from an external vcf for the script SingleCellGenotype.py

SCOMATIC=/scratch/devel/pnieto/projects/CSF/data/SComatic/SComatic-main
project=CSF
subproject=CSF_01
sample=4839
output_dir=/scratch/devel/pnieto/projects/$project/data/$subproject/output/05_SComatic/$sample
VCF=/scratch/devel/pnieto/projects/$project/data/$subproject/jobs/$sample/4839_mutation.vcf
grep -v '^#' $VCF | awk -F'\t' -v OFS='\t' '{print $1,$2,$2,$4,$5,".","Cell_type",".",".",".",".",".",".","."}' > /scratch/devel/pnieto/projects/$project/data/$subproject/jobs/$sample/Mutations.other_caller.tsv

cell_type_bam=/path/to/cell_type.bam
META=/path/to/cell_type_annotations.tsv
REF=/path/to/reference_genome.fasta

python $SCOMATIC/scripts/SingleCellGenotype/SingleCellGenotype.py --bam $cell_type_bam  \
    --infile Mutations.other_caller.tsv \
    --nprocs 16  \
    --meta $META   \
    --outfile Mutations.other_caller.single_cell_genotype.tsv  \
    --tmp_dir temp  \
    --ref $REF