import pandas as pd

class repeat_data:
    def __init__(self, repeat_info, sample_data, hom_sample_count, hom_allele_count, num_snps, snp_data=None, str_data=None):

        self.repeat_info = repeat_info
        self.chrom = repeat_info['chrom']
        self.pos = repeat_info['pos']
        self.end = repeat_info['end']
        self.ref = repeat_info['ref']
        self.repeat_unit = repeat_info['repeat_unit']
        self.sample_data = sample_data
        self.hom_sample_count= hom_sample_count
        self.hom_allele_count= hom_allele_count
        self.num_snps = num_snps

        self.snp_data = snp_data

        if str_data is None:
            ## If no SNP data is given, the STRs can't be predicted and will be set to '.'
            if snp_data is None:
                hap1_df = pd.DataFrame(['.' for sample in self.sample_data.index], columns=['str_length'], index=[f'{sample}+1' for sample in self.sample_data.index])
                hap2_df = pd.DataFrame(['.' for sample in self.sample_data.index], columns=['str_length'], index=[f'{sample}+2' for sample in self.sample_data.index])
                self.str_data = pd.concat([hap1_df, hap2_df])
            ## If SNP data is provided, the STRs will be predicted later
            else:
                self.str_data = pd.DataFrame(columns=['str_length'], dtype='float')
        else:
            self.str_data = str_data

    def add_str_data(self, str_data):
        self.str_data = pd.concat([self.str_data, str_data])

    def match_snps(self, other):
        self_snps = set(self.snp_data.columns)
        other_snps = set(other.snp_data.columns)
        
        shared_snps = list(other_snps.intersection(self_snps))
        self.snp_data = self.snp_data[shared_snps]
        other.snp_data = other.snp_data[shared_snps]

    def subset_by_population(self, population):
        snp_data = self.snp_data.loc[self.sample_data['population'] == population]
        sample_data = self.sample_data.loc[self.sample_data['population'] == population]
        str_data = self.str_data.loc[self.sample_data['population'] == population]

        repeat_data_subset = repeat_data(self.repeat_info, sample_data, self.hom_sample_count, self.hom_allele_count, self.num_snp, snp_data=snp_data, str_data=str_data)

        return(repeat_data_subset)

    def order_str_data(self):
        self.str_data = self.str_data.reindex(self.snp_data.index)

    def remove_constants(self):
        '''
        remove any SNP columns that have a single allele across the population
        '''
        self.snp_data = self.snp_data.loc[:, (self.snp_data != self.snp_data.iloc[0]).any()]

    def create_plotting_df(self):
        return(self.sample_data.merge(self.str_data, right_index=True, left_index=True))       

    def create_vcf_entry(self, method = 'eh'):

        ## Determine the alternate alleles present in the samples
        if method == 'gangstr':
            alt_alleles = [(self.repeat_unit * int(repeat_count)).upper() for repeat_count in set(self.str_data['str_length']) if not pd.isna(repeat_count) and repeat_count != '.']
        if method == 'eh':
            alt_alleles = [f'<STR{int(repeat_count)}>' for repeat_count in set(self.str_data['str_length']) if not pd.isna(repeat_count) and repeat_count != '.']

        if len(alt_alleles) == 0:
            alt_alleles = ['.']
           
        ## Remove the reference allele so it can be added to the front of the list and the
        ## indecies can be used to create the genotypes for the VCF
        try:
            if method == 'gangstr':
                alt_alleles.remove(self.ref)
        except ValueError:
            pass

        all_alleles = [self.ref] + alt_alleles

        all_samples = list(set([sample.split('+')[0] for sample in self.sample_data.index]))
        all_samples.sort()

        sample_info_columns = []
        for sample in all_samples:
            repeat_count1 = self.str_data.loc[f'{sample}+1','str_length']
            if not pd.isna(repeat_count1) and repeat_count1 != '.':
                if method == 'gangstr':
                    allele_idx1 = all_alleles.index((self.repeat_unit * int(repeat_count1)).upper())
                if method == 'eh':
                    allele_idx1 = all_alleles.index(f'<STR{int(repeat_count1)}>')
            else: 
                allele_idx1 = '.'

            repeat_count2 = self.str_data.loc[f'{sample}+2','str_length'] 
            if not pd.isna(repeat_count2) and repeat_count2 != '.':
                if method == 'gangstr':
                    allele_idx2 = all_alleles.index((self.repeat_unit * int(repeat_count2)).upper()) 
                if method == 'eh':
                    allele_idx2 = all_alleles.index(f'<STR{int(repeat_count2)}>') 
            else:
                allele_idx2 = '.'
            
            if pd.isna(repeat_count1):
                repeat_count1 = '.'
            if pd.isna(repeat_count2):
                repeat_count2 = '.'

            if method == 'gangstr':
                sample_info_string = f'{allele_idx1}/{allele_idx2}:{repeat_count1},{repeat_count2}'
            if method == 'eh':
                sample_info_string = f'{allele_idx1}/{allele_idx2}:{repeat_count1}/{repeat_count2}'
 
            sample_info_columns.append(sample_info_string)

        variant_info_string = f'{self.chrom}\t{self.pos}\t{self.chrom}.{self.pos}\t{self.ref}\t{",".join(alt_alleles)}\t.\t.\tEND={self.end};RU={self.repeat_unit};NTH={self.hom_sample_count};NTA={self.hom_allele_count};NTS={self.num_snps}\tGT:REPCN\t'
        sample_info_string = '\t'.join(sample_info_columns)
        return(variant_info_string + sample_info_string + "\n")

