import numpy as np
import pandas as pd

class vcf_data:
    def __init__(self, variants_info, samples, genotype_info):

        self.variants_info = variants_info
        self.samples = list(samples)
        self.genotype_info = genotype_info

    def get_variant_genotypes_in_range(self, chrom, start, end, ru=None):
        if ru:
            ## Get the index of the variants that fit the criteria (to subset the genotype data)
            variant_idx = self.variants_info.index[(self.variants_info['CHROM'] == chrom) & (self.variants_info['POS'] >= start) & (self.variants_info['POS'] <= end) & (self.variants_info['RU'] == ru)] 
            ## Subset the variant into to the ones that fit the criteria
            variant_in_range_variant_info = self.variants_info[(self.variants_info['CHROM'] == chrom) & (self.variants_info['POS'] >= start) & (self.variants_info['POS'] <= end) & (self.variants_info['RU'] == ru)]
        else:
            ## Get the index of the variants that fit the criteria (to subset the genotype data)
            variant_idx = self.variants_info.index[(self.variants_info['CHROM'] == chrom) & (self.variants_info['POS'] >= start) & (self.variants_info['POS'] <= end)] 
            ## Subset the variant into to the ones that fit the criteria
            variant_in_range_variant_info = self.variants_info[(self.variants_info['CHROM'] == chrom) & (self.variants_info['POS'] >= start) & (self.variants_info['POS'] <= end)] 

        ## Subset the genotype data based on the indecies of the varints that fit the criteria
        variant_in_range_genotypes = self.genotype_info[variant_idx]

        ## Create the new (subset) vcf_data object
        new_vcf_data = vcf_data(variant_in_range_variant_info, self.samples, variant_in_range_genotypes)
        return(new_vcf_data)

    def subset_samples(self, sample_list):
        ## Determine the indecies of the samples to keep 
        keep_idx = [idx for idx, sample in enumerate(self.samples) if sample in sample_list]

        ## Subset the genotype data with the samples to keep
        genotype_data = self.genotype_info[:, keep_idx]

        ## Subset the sample data with the samples to keep
        sample_data = [sample for sample in self.samples if sample in sample_list]

        ## Create the new (subset) vcf_data object
        new_vcf_data = vcf_data(self.variants_info, sample_data, genotype_data)
        return(new_vcf_data)

    def format_for_repeat_data(self):
        var_id_chr = list(self.variants_info['CHROM'])
        var_id_pos = list(self.variants_info['POS'])
        var_id_ref = list(self.variants_info['REF'])
        var_id_alt = list(self.variants_info['ALT_1'])
        var_id_data = [f'{chrom}:{pos}:{ref}:{alt}' for chrom,pos,ref,alt in zip(var_id_chr, var_id_pos, var_id_ref, var_id_alt)] 

        ## allel will create a 3d structure if the genotypes a split by a ',' as they are for the SNP
        ## If they are split by a '/', it treats them as a single string, here I am creating the 3d structure
        if len(self.genotype_info.shape) < 3:
            genotype_data = []
            for variant in self.genotype_info:
                variant_genotype_data = []
                for sample_geno in variant:
                    hap1, hap2 = sample_geno.split('/')
                    if hap1 != '.':
                        hap1 = int(round(float(hap1)))
                    if hap2 != '.':
                        hap2 = int(round(float(hap2)))
                    variant_genotype_data.append([hap1, hap2])
                genotype_data.append(variant_genotype_data)
            genotype_data = np.array(genotype_data)
            #genotype_data = np.array([[list(map(int, list(map(float, variant.split('/'))))) if variant!="./." else ['.','.'] for variant in sample] for sample in self.genotype_info])
        else:
            genotype_data = self.genotype_info
        hap_1s = [genotype_data[:, i, 0]  for i in range(genotype_data.shape[1])] 
        hap_2s = [genotype_data[:, i, 1]  for i in range(genotype_data.shape[1])] 

        ## Create data frames from the 2D lists with variant IDs as column names and sample id as row names
        hap1_df = pd.DataFrame(hap_1s, columns=var_id_data, index=[f'{sample}+1' for sample in self.samples])
        hap2_df = pd.DataFrame(hap_2s, columns=var_id_data, index=[f'{sample}+2' for sample in self.samples])
        hap_df = pd.concat([hap1_df, hap2_df])
        return(hap_df)
