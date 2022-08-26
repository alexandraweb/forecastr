
rule repeat_by_chr:
    input:
        config["repeat_file"]
    output:
        temp("{outdir}/chr{chr}_repeat_file.txt")
    params:
        chr="{chr}"
    shell:
        """
        awk '$1=="chr{params.chr}"' < {input} | uniq > {output}
        """

rule trainingSnpVcf_subsetVariantsAndSamples:
    input:
        vcf=config["training_snp_vcf"],
        array_coords=config["array_coords"]
    output:
        temp("{outdir}/chr{chr}_training_snps.vcf")
    params:
        chr="{chr}",
    resources:
        time="48:00:00"
    container:
        "docker://aleweb/forecastr:latest"
    shell:
        """
        mkdir -p {TMPDIR}/t_snp_{params.chr}

        ##############
        # Merge VCFs if needed
        ##############

        file_count=0
        for i in {input.vcf}; do file_count=`expr $file_count + 1`; done
        if [ $file_count -gt 1 ]; then
          bcftools merge -m all -o {TMPDIR}/t_snp_{params.chr}/training_snps_{params.chr}.vcf {input.vcf}
          bgzip {TMPDIR}/t_snp_{params.chr}/training_snps_{params.chr}.vcf 
          tabix -p vcf {TMPDIR}/t_snp_{params.chr}/training_snps_{params.chr}.vcf.gz
        else
          cp {input.vcf} {TMPDIR}/t_snp_{params.chr}/training_snps_{params.chr}.vcf.gz
        fi

        ##############
        # Subset variants
        ##############
 
        ## Subset the SNP VCF to region surrounding the STR and remove indels 
        bcftools view -f 'PASS' {TMPDIR}/t_snp_{params.chr}/training_snps_{params.chr}.vcf.gz | \
        bcftools view -v 'snps' | \
        bcftools view -V 'indels' \
        > {TMPDIR}/t_snp_{params.chr}/chr{params.chr}_hq_snps.vcf
        
        ## Subset the prediction SNPs to those on the SNP array
        bedtools intersect -header -wa -a {TMPDIR}/t_snp_{params.chr}/chr{params.chr}_hq_snps.vcf -b {input.array_coords} \
           > {output}

        rm -r {TMPDIR}/t_snp_{params.chr}/
        """

rule predictionSnpVcf_subsetVariantsAndSamples:
    input:
        vcf=config["predict_snp_vcf"],
        population=config["predict_sample_populations"],
        array_coords=config["array_coords"]
    output:
       temp("{outdir}/chr{chr}_prediction_snps.vcf")
    params:
        chr="{chr}",
        sample_file=config["predict_samples"]
    resources:
        time="4:00:00"
    container:
        "docker://aleweb/forecastr:latest"
    shell:
        """
        mkdir -p {TMPDIR}/p_snp_{params.chr}


        ##############
        # Subset Variants 
        ##############

        ## Subset the prediction SNPs to those on the SNP array
        #** If the cohort is the prediction dataset - this step shouldn't change anything
        #** If a WGS dataset is the prediction dataset - this step will simulate an array 
        bedtools intersect -header -wa -a {input.vcf} -b {input.array_coords} \
            > {TMPDIR}/p_snp_{params.chr}/chr{params.chr}_array_snps.vcf

        ##############
        # Subset samples 
        ##############

        ## Subset samples based on the file if supplied in the config file
        if [ -f {params.sample_file} ]; then
            bcftools view -S {params.sample_file} {TMPDIR}/p_snp_{params.chr}/chr{params.chr}_array_snps.vcf > {TMPDIR}/p_snp_{params.chr}/chr{params.chr}_array_configSamples_snps.vcf
            grep -f {params.sample_file} {input.population} > {TMPDIR}/p_snp_{params.chr}/population_info.txt
        else
            cp {TMPDIR}/p_snp_{params.chr}/chr{params.chr}_array_snps.vcf {TMPDIR}/p_snp_{params.chr}/chr{params.chr}_array_configSamples_snps.vcf
            cp {input.population} > {TMPDIR}/p_snp_{params.chr}/population_info.txt
        fi

        ## Subset the cohort array to samples with population information
        #cut -f 1 {input.population} > {TMPDIR}/p_snp_{params.chr}/pop_samples.txt
        cut -f 1 {TMPDIR}/p_snp_{params.chr}/population_info.txt > {TMPDIR}/p_snp_{params.chr}/pop_samples.txt
        bcftools view \
            --force-samples -S {TMPDIR}/p_snp_{params.chr}/pop_samples.txt \
            {TMPDIR}/p_snp_{params.chr}/chr{params.chr}_array_configSamples_snps.vcf \
            > {output}

        rm -r {TMPDIR}/p_snp_{params.chr}/
        """


rule combine_files:
    input:
        training_str_vcf=expand(config["training_str_vcf"], method=method, allow_missing=True),
        training_pop=config["training_sample_populations"],
        predict_pop=config["predict_sample_populations"]
    output:
        vcf=temp("{outdir}/training_strs_{chr}.vcf.gz"),
        population=temp("{outdir}/chr{chr}_pop.txt") 
    params:
        chr="{chr}",
        out="{outdir}"
    container:
        "docker://aleweb/forecastr:latest"
    shell:
        """
        mkdir -p {TMPDIR}/comb_{params.chr}

        ##############
        # Merge VCFs if needed
        ##############
        file_count=0
        for i in {input.training_str_vcf}; do file_count=`expr $file_count + 1`; done
        if [ $file_count -gt 1 ]; then
          bcftools merge -m all -o {params.out}/training_strs_{params.chr}.vcf {input.training_str_vcf}
          bgzip {params.out}/training_strs_{params.chr}.vcf 
          tabix -p vcf {params.out}/training_strs_{params.chr}.vcf.gz
        else
          cp {input.training_str_vcf} {params.out}/training_strs_{params.chr}.vcf.gz
        fi

        ##############
        # Combine all population files
        ##############
        cat {input.training_pop} {input.predict_pop} | sort | uniq > {params.out}/chr{params.chr}_pop.txt


        rm -r {TMPDIR}/comb_{params.chr}/
        """

rule predictStrVcf:
    input:
        training_str_vcf="{outdir}/training_strs_{chr}.vcf.gz",
        training_snp_vcf="{outdir}/chr{chr}_training_snps.vcf", 
        predict_snp_vcf="{outdir}/chr{chr}_prediction_snps.vcf", 
        population="{outdir}/chr{chr}_pop.txt", 
        repeat_file="{outdir}/chr{chr}_repeat_file.txt"
    output:
        temp("{outdir}/chr{chr}_forecastr.vcf")
    params:
        chr="{chr}",
        out="{outdir}",
        prediction_method=config['prediction_method']
    conda:
        "../envs/new_forecastr.yml"
    threads: 16
    resources:
        mem_mb=100000,
        time="100:00:00"
    shell:
        """
        ##############
        # Run imputation script
        ##############
        echo "submitting python command"
        python workflow/scripts/predict_repeat_lengths.py \
            {input.training_str_vcf} \
            {input.training_snp_vcf} \
            {input.predict_snp_vcf} \
            {input.repeat_file}\
            {input.population} \
            {output} \
            -m {params.prediction_method}

        """
