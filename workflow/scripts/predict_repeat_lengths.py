import argparse
import sys
import allel

import datetime
import numpy as np
import pandas as pd

from repeat_data import repeat_data
from vcf_data import vcf_data
from exceptions import noMatchingStrError,multipleMatchingStrsError,notEnoughHomSamplesError,notEnoughSnpsError

from sklearn.linear_model import LinearRegression
import statsmodels.miscmodels.ordinal_model as ordinal_model

from concurrent.futures import ProcessPoolExecutor as Pool
from collections import namedtuple
import random
#from mpi4py.futures import MPIPoolExecutor


def extract_data(vcf, str_vcf=False):
 
    if str_vcf:
        variant_info = allel.vcf_to_dataframe(vcf, numbers={'ALT': 30}, fields=['samples', 'variants/CHROM', 'variants/RU', 'variants/POS','variants/REF', 'variants/ALT', 'variants/END'])
        genotype_info = allel.read_vcf(vcf, fields=['calldata/REPCN', 'calldata/GT', 'samples'])
        vcf_obj = vcf_data(variant_info, genotype_info['samples'], genotype_info['calldata/REPCN'])
        
    else:
        variant_info = allel.vcf_to_dataframe(vcf, numbers={'ALT': 30})
        genotype_info = allel.read_vcf(vcf, fields=['calldata/GT', 'samples'])
        vcf_obj = vcf_data(variant_info, genotype_info['samples'], genotype_info['calldata/GT'])
    return(vcf_obj)

def create_file_dict(dict_file):
    '''
    create a dictionary from a file with the first column as the keys and 
    second column as the values
    (This is used for the population dictioanry)
    '''
    file_dict = {}
    with open(dict_file, 'r') as f:
        for line in f.readlines():
            file_dict[line.split('\t')[0]] = line.split('\t')[1].strip()
    return(file_dict)

def ordinal_linear_model_estimation(ref_snp_data, ref_str_data, cohort_snp_data):
    '''
    Predict the STR length of the cohort based on an ordinal linear model trained on the 
    reference data
    '''
    ## Creating a data ordered data type for the ordinal model
    repeat_length_dt = pd.CategoricalDtype(categories=range(501), ordered=True)

    ## Convert the reference STR lengths to the categorical data type
    ref_str_data = pd.Series(ref_str_data['str_length'], dtype=repeat_length_dt)


    ## This dictionary is used to translate the cutoff group results to the
    ## correct repeat count range (all results will start and 0 and need to start
    ## at the minimum repeat length used to train the model)
    convert_dict = {}
    values = sorted(set(ref_str_data.values))
    for i in range(len(values)):
        convert_dict[i] = values[i]

    ## Set up the model variables
    X = ref_snp_data
    y = ref_str_data.to_frame()

    ## Check if there is only one STR allele in the training set
    constant_y = (y['str_length'] == y['str_length'][0]).all()

    ## If there is a single STR allele - assign that allele to all prediction
    ## samples without creating the model
    if constant_y:
        cohort_str_data = pd.DataFrame([y['str_length'][0]] * len(cohort_snp_data.index), columns=['str_length'], index=cohort_snp_data.index)
    else:
        ## Train the model
        model_log = ordinal_model.OrderedModel(y, X, distr='logit', hasconst=False)
        res_log = model_log.fit(method='bfgs', disp=False)

        ## Use the model to predict the cohort STR length
        try:
            cohort_str_predictions = res_log.model.predict(res_log.params, exog=cohort_snp_data)
            cohort_str_data = cohort_str_predictions.argmax(1)

            # Convert values to the right range
            cohort_str_data = [convert_dict[value] for value in cohort_str_data]

        except ValueError:
            cohort_str_data = pd.DataFrame(['.'] * len(cohort_snp_data.index), columns=['str_length'], index=cohort_snp_data.index)

        ## Create a dataframe with the results
        cohort_str_data = pd.DataFrame(cohort_str_data, columns=['str_length'], index=cohort_snp_data.index)

    return(cohort_str_data)

def random_selection(ref_snp_data, ref_str_data, cohort_snp_data):

    ## Get a list of the STR lengths present in the reference set
    ref_str_data = list(ref_str_data['str_length'])
 
    ## Select an STR length randomly from the reference set for each haplotype in the cohort set
    cohort_str_data = [random.choice(ref_str_data) for value in cohort_snp_data.index]
    cohort_str_data = pd.DataFrame(cohort_str_data, columns=['str_length'], index=cohort_snp_data.index)
    return(cohort_str_data)

def predict_str_length2(ref, cohort, pop_dict, prediction_method="ordinal", by_pop=False):

    ## Find common SNPs between datasets
    cohort.match_snps(ref)
   
    if by_pop:
        ## Estiamte STR lengths in the cohort (done by population)
        for population in set(pop_dict.values()):

            ## Subset the samples by population
            pop_cohort = cohort.subset_by_population(population)
            pop_the_ref = ref.subset_by_population(population)
            
            cohort_str_data = predict_str_length(pop_ref, pop_cohort, prediction_method)

            ## Add the predicted STRs to the cohort object
            cohort.add_str_data(cohort_str_data)
    else:
        cohort_str_data = predict_str_length(ref, cohort, prediction_method)
        cohort.add_str_data(cohort_str_data)

    ## Order the STR data to match the sample order of the SNP data
    cohort.order_str_data()

    

def predict_str_length(ref, cohort, prediction_method):
    '''
    predict the cohort STR data population and add to the data to the cohort object
    '''
    print(prediction_method)
    ref.order_str_data()

    ## Continue as long as there are samples for the given population
    if cohort.snp_data.shape[0] > 0 and ref.snp_data.shape[0] > 0:

        ## Remove any SNPs without variability across the training sub-population
        ## (also remove these from the prediction set)
        ref.remove_constants()
        ref.match_snps(cohort)

        ## Predict cohort STR lengths
        if cohort.snp_data.shape[1] > 0 and ref.snp_data.shape[1] > 0:
            if prediction_method == "ordinal":
                cohort_str_data = ordinal_linear_model_estimation(ref.snp_data, ref.str_data, cohort.snp_data)
            elif prediction_method == "random":
                cohort_str_data = random_selection(ref.snp_data, ref.str_data, cohort.snp_data)
            ## TODO: need to handle the case where nothing is currently getting returned
            return(cohort_str_data)
        else:
            print('not enough SNPs')
    else:
        print('not enough samples')


#global_lock = threading.Lock()
vcf_entries = []


def mapped_create_vcf_entry(named_tuple_input):
    vcf_entry = create_vcf_entry(named_tuple_input.repeat, named_tuple_input.train_str_vcf_obj, named_tuple_input.train_snp_vcf_obj, named_tuple_input.predict_snp_vcf_obj, named_tuple_input.pop_dict, named_tuple_input.samples, named_tuple_input.prediction_method)
    return(vcf_entry)

def create_vcf_entry(repeat, train_str_vcf_obj, train_snp_vcf_obj, predict_snp_vcf_obj, pop_dict, samples, prediction_method):
    ## Save the STR information
    chrom, start, stop, ru = repeat.split()
    start = int(start)
    stop = int(stop)

    repeat_info = {'repeat_unit': ru, 'chrom': chrom, 'pos': start, 'end': stop, 'ref': '.'}

    try:
        #########
        ## Train STR
        #########
        ## Get the STR variant
        variant_train_str_vcf_obj = train_str_vcf_obj.get_variant_genotypes_in_range(chrom, start, stop, ru)
        num_var, num_samples = variant_train_str_vcf_obj.genotype_info.shape
        if num_var == 0:
            print('No matching STR was found in the training set')
            raise noMatchingStrError
            #continue 
        if num_var > 1:
            print('Multiple STR matches were found')
            raise multipleMatchingStrsError
            #continue
        
        ## Determine the homozyougs samples
        ##   Note: there should only be one VCF entry for the current repeat but it is still an ndarray
        geno_list = [genotype.split('/') for genotype in variant_train_str_vcf_obj.genotype_info[0]]
        hom_sample_idx = [idx for idx, geno in enumerate(geno_list) if geno[0] == geno[1] and geno[0] != '.']
        num_hom_alleles = len(set([geno[0] for idx, geno in enumerate(geno_list) if geno[0] == geno[1] and geno[0] != '.']))
        hom_samples = [train_str_vcf_obj.samples[i] for i in hom_sample_idx]
        num_hom_samples = len(hom_samples)
        if num_hom_samples < 3:
            print('Not enough homozygous samples were found in the training set')
            raise notEnoughHomSamplesError
            #continue 
        ## Subset the data to the homozygous samples
        hom_train_str_vcf_obj = variant_train_str_vcf_obj.subset_samples(hom_samples)
        
        #########
        ## Train SNP
        #########
        ## Get the variants surrounding the STR
        variants_train_snp_vcf_obj = train_snp_vcf_obj.get_variant_genotypes_in_range(chrom, start-100000, stop + 100000)
        num_var, num_samples, ploidy = variants_train_snp_vcf_obj.genotype_info.shape
        num_training_snps = num_var
        if num_var < 3:
            print('Not enough SNPs were found surrounding the STR variant in the training set')
            raise notEnoughSnpsError
            #continue
        
        ## Subset the data to the homozygous samples
        hom_train_snp_vcf_obj = train_snp_vcf_obj.subset_samples(hom_samples)
        
        #########
        ## Prediction SNP
        #########
        ## Get the variants surrounding the STR
        variants_predict_snp_vcf_obj = predict_snp_vcf_obj.get_variant_genotypes_in_range(chrom, start - 100000, stop + 100000)
        num_var, num_samples, ploidy = variants_predict_snp_vcf_obj.genotype_info.shape
        num_prediction_snps = num_var
        if num_var < 3:
            print('Not enough SNPs were found surrounding the STR variant in the prediction set')
            raise notEnoughSnpsError
            #continue 
        
        ## Convert formats
        train_str_data = hom_train_str_vcf_obj.format_for_repeat_data()
        train_str_data.columns = ['str_length']

        train_snp_data = hom_train_snp_vcf_obj.format_for_repeat_data()

        train_sample_data = pd.DataFrame([[pop_dict[sample.split('+')[0]]] for sample in train_snp_data.index], columns=['population'], index=train_snp_data.index)
        
        predict_snp_data = variants_predict_snp_vcf_obj.format_for_repeat_data()
        predict_sample_data = pd.DataFrame([[pop_dict[sample.split('+')[0]]] for sample in predict_snp_data.index], columns=['population'], index=predict_snp_data.index)
        
        ref = repeat_data(repeat_info, train_sample_data, num_hom_samples, num_hom_alleles, num_training_snps, snp_data=train_snp_data, str_data=train_str_data)
        cohort = repeat_data(repeat_info, predict_sample_data, num_hom_samples, num_hom_alleles, num_prediction_snps, snp_data=predict_snp_data)
        
        predict_str_length2(ref, cohort, pop_dict, prediction_method)

    except (noMatchingStrError, multipleMatchingStrsError, notEnoughHomSamplesError, notEnoughSnpsError) as e:
        predict_sample_data = pd.DataFrame([[pop_dict[sample.split('+')[0]]] for sample in predict_snp_vcf_obj.samples], columns=['population'], index=predict_snp_vcf_obj.samples)
        cohort = repeat_data(repeat_info, predict_sample_data, -1, -1, -1)

    vcf_entry = cohort.create_vcf_entry()
    return(vcf_entry)

def write_vcf_file(samples, vcf_entries, out_file):
    with open(out_file,'w') as o:
        o.write('##fileformat=VCFv4.2\n')
        o.write(f'##fileDate={datetime.datetime.now().strftime("%Y%m%d")}\n')
        o.write('##source=forecastr\n')
        o.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">\n')
        o.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat Unit">\n')
        o.write('##INFO=<ID=NTH,Number=1,Type=String,Description="Number of training haplotypes">\n')
        o.write('##INFO=<ID=NTA,Number=1,Type=String,Description="Number of training alleles">\n')
        o.write('##INFO=<ID=NTS,Number=1,Type=String,Description="Number of training SNPs">\n')
        o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        o.write('##FORMAT=<ID=REPCN,Number=1,Type=String,Description="Repeat_count">\n')
        for i in range(1,23):
            o.write(f'##contig=<ID=chr{i}>\n')
        o.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n')
        o.write(''.join(vcf_entries))
    

if __name__ == '__main__':
    # Configure the argument parser
    parser = argparse.ArgumentParser(
      description='estimate STR lengths from surrounding SNP haplotypes')

    parser.add_argument('train_str_vcf', help='/path/to/train_str.vcf')
    parser.add_argument('train_snp_vcf', help='/path/to/train_snp.vcf')
    parser.add_argument('predict_snp_vcf', help='/path/to/predict_snp.vcf')
    parser.add_argument('repeat_file', help='/path/to/repeat_file.txt (foramt CHR\tSTART\tEND\tREPEAT_UNIT)')
    parser.add_argument('pop_file', help='/path/to/population_file.txt (format SAMPLE_NAME\tPOP_NAME)')
    parser.add_argument('out_path', help='/path/to/output_prefix')
    parser.add_argument('--prediction_method', '-m', default="ordinal", help='prediction method type (ordinal, random, linear; Default: ordinal)')

    args = parser.parse_args()

    train_str_vcf = args.train_str_vcf
    train_snp_vcf = args.train_snp_vcf
    predict_snp_vcf = args.predict_snp_vcf
    repeat_file = args.repeat_file
    pop_file = args.pop_file
    out = args.out_path
    prediction_method = args.prediction_method


    ## Read in the prediction STR vcf
    train_str_vcf_obj = extract_data(train_str_vcf, str_vcf=True)
    # Read in the training SNP vcf
    train_snp_vcf_obj = extract_data(train_snp_vcf)
    ## Read in the prediction SNP vcf
    predict_snp_vcf_obj = extract_data(predict_snp_vcf)
    
    print('Finished reading in VCF files')

    pop_dict = create_file_dict(pop_file)
    
    samples = predict_snp_vcf_obj.samples
    samples.sort()
    samples = "\t".join(samples)

    print('Finished reading info files')
        
    vcf_entry_input = namedtuple('vcf_entry_input', ['repeat',
        'train_str_vcf_obj', 'train_snp_vcf_obj', 'predict_snp_vcf_obj',
        'pop_dict', 'samples', 'prediction_method'])
        
    repeat_vcf_entry_inputs = [vcf_entry_input(repeat=repeat, \
          train_str_vcf_obj=train_str_vcf_obj, train_snp_vcf_obj=train_snp_vcf_obj, \
          predict_snp_vcf_obj=predict_snp_vcf_obj, pop_dict=pop_dict, samples=samples, prediction_method=prediction_method) \
          for repeat in open(repeat_file, 'r')]

    vcf_entries = []
    count = 0
    for repeat in repeat_vcf_entry_inputs:
        count += 1
        vcf_entries.append(mapped_create_vcf_entry(repeat))

#    executor = Pool()
#    vcf_entries = list(executor.map(mapped_create_vcf_entry, repeat_vcf_entry_inputs))
#    executor.shutdown()
    
    write_vcf_file(samples, vcf_entries, out)
