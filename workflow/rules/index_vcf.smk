rule zip_vcf:
    input:
        "{file}.vcf"
    output:
        vcf="{file}.vcf.gz",
    shell:
        """
        module load Bioinformatics
        module load tabix

        bgzip {input} 
        """
        
rule index_vcf:
    input:
        "{file}.vcf.gz"
    output:
        tbi="{file}.vcf.gz.tbi"
    shell:
        """
        module load Bioinformatics
        module load tabix

        tabix -p vcf {input}
        """
