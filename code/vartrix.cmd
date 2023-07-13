project=CSF
subproject=CSF_01
sample=4839
vcf_file=/scratch/devel/pnieto/projects/$project/data/$subproject/jobs/$sample/4839_mutation.vcf
bam_file=/scratch/devel/pnieto/projects/$project/data/$subproject/jobs/$sample/$sample/outs/per_sample_outs/$sample/count/sample_alignments.bam
bc_file=/scratch/devel/pnieto/projects/$project/data/$subproject/jobs/$sample/$sample/outs/per_sample_outs/$sample/count/sample_filtered_feature_bc_matrix/barcodes.tsv
REF=/scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa
output_dir=/scratch/devel/pnieto/projects/$project/data/$subproject/output/05_Vartrix/$sample

vartrix -v $vcf_file \
        -b $bam_file \
        -f $REF \
        -c $bc_file \
        --threads=6 \
        -o $output_dir