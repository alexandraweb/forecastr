
rule all:
    input:
        expand("{outdir}/prediction_snp_array_coords.bed", outdir=config['outdir'])
    


rule array_coord_bed:
    input:
        config["predict_snp_vcf"]
    output:
        "{outdir}/prediction_snp_array_coords.bed"
    shell:
        """
        module load Bioinformatics
        module load bcftools

        bcftools view -H {input} | awk 'BEGIN {{OFS="\t"}} {{print $1,$2-1,$2}}' \
          > {output}
        """
