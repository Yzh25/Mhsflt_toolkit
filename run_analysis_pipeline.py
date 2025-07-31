#!/usr/bin/env python3
#-*-coding:utf-8-*-

import os
import argparse
import yaml
import re,copy
import os,sys
import subprocess
import shutil
from datetime import datetime

def load_yaml(yaml_file):
    if not os.path.exists(yaml_file):
        raise FileNotFoundError(f'YAML file not exist: {yaml_file}')
    
    with open(yaml_file, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    return data


def check_nested_empty(data, path=''):
    if isinstance(data, dict):
        for k, v in data.items():
            check_nested_empty(v, f"{path}.{k}" if path else k)
    elif isinstance(data, (list, tuple)):
        for i, v in enumerate(data):
            check_nested_empty(v, f"{path}[{i}]")
    elif data is None:
        if path not in ['reference.additional_files', 'barcodes.barcode_num', 'reference.additional_STAR_params', 'barcodes.barcode_file']:
            raise ValueError(f'Find empty value: {path}')


def check_species(data):
    species=data['sample']['sample_species'].lower()
    if species not in ['human', 'mouse']:
        raise ValueError(f'Species now only support "human" and "mouse".')
    

def make_dir(data):
    out_path=data['out_dir']
    if os.path.exists(out_path):
        raise OSError('Out dirctory in yaml is exist, please delete it')
    
    for dirs in ['data','config','analysis','results']:
        os.makedirs(f'{out_path}/{dirs}') 


def validate_type_id(data):
    if data['sample']['sample_type'].lower() != 'manual' and data['sample']['sample_type'].lower() != 'auto':
        raise ValueError('Wrong sample type in yaml file, must be "manual" or "auto".')
    
    if data['sample']['sample_type'] == '' or data['sample']['sample_id'] == '':      
        raise ValueError('sample type and sample id must set in yaml file')

    if data['sample']['sample_type'].lower() == 'manual':
        s_id=str(data['sample']['sample_id']).split(',')
        for ind in s_id:
            if int(ind) < 1 or int(ind) > 24:
                raise ValueError('Manual version sample id must be set within the range of 1-24, with multiple ids separated by commas.')                
    else:
        try: 
            ind=int(data['sample']['sample_id'])
            if ind < 1 or ind > 12:
                raise ValueError('Automatic version sample ID must be set within the range 1-12.')
        except:   
            raise TypeError(f'Automatic version sample ID only support single sample id.')


def create_barcode(data):
    sample_type=data['sample']['sample_type'].lower()
    sample_id=str(data['sample']['sample_id'])
    out_path=data['out_dir']
    script_path=data['zUMIs_directory']

    if sample_type=='manual':
        with open(f'{script_path}/manual_barcode_list.yaml', 'r', encoding='utf-8') as y:
            barcode_set = yaml.safe_load(y)
        sample_id=sample_id.split(',')
        with open(f'{out_path}/config/expect_barcode.tsv','w') as pipe_file, open(f'{out_path}/config/expect_bin_barcode.tsv','w') as bin_file, open(f'{out_path}/config/expect_id_barcode.tsv','w') as summary_file:
            print('\t'.join(['wellID','umi_barcodes','internal_barcodes']),file=summary_file)
            for c in sample_id:
                bc_t_i5=barcode_set[c][0];bc_n_i5=barcode_set[c][1];bc_n_i7=barcode_set[c][2]
                for k,v in zip(bc_t_i5,bc_n_i5):
                    for j in bc_n_i7:
                        print('\t'.join([k+j,v+j]),file=pipe_file)
                for kk,vv in zip(bc_t_i5,bc_n_i5):
                    for jj in bc_n_i7:
                        print('\t'.join([kk+jj,vv+jj]),file=bin_file)
                        print('\t'.join(['A'+c,kk+jj,vv+jj]),file=summary_file)
                        break
                    break
    else:
        with open(f'{script_path}/auto_barcode_list.yaml', 'r', encoding='utf-8') as y:
            barcode_set = yaml.safe_load(y)
        with open(f'{out_path}/config/expect_barcode.tsv','w') as pipe_file, open(f'{out_path}/config/expect_id_barcode.tsv','w') as summary_file:
            print('\t'.join(['wellID','umi_barcodes','internal_barcodes']),file=summary_file)
            for k in barcode_set[f'plate{sample_id}']:
                print('\t'.join(barcode_set[f'plate{sample_id}'][k]),file=pipe_file)
                print(k+'\t'+'\t'.join(barcode_set[f'plate{sample_id}'][k]),file=summary_file)    
    
    
def check_file_exists(data):
    files_to_check=[data['sequence_files']['file1']['name'],
                   data['sequence_files']['file2']['name'],
                   data['reference']['STAR_index'],
                   data['reference']['GTF_file']]
    for filepath in files_to_check:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f'{filepath} not exist, please check!')
           

def process_fq(data):
    sample_type = data['sample']['sample_type'].lower()
    script_path = data['zUMIs_directory']
    fq1 = data['sequence_files']['file1']['name']
    fq1_name = os.path.basename(fq1)
    fq2 = data['sequence_files']['file2']['name']
    fq2_name = os.path.basename(fq2) 
    out_path=data['out_dir']

    os.symlink(f'{fq1}', f'{out_path}/data/{fq1_name}')

    if sample_type == 'manual':
        out_fq2_name = f'{fq2_name}.process.fq.gz'
        ham_dist=data['barcodes']['BarcodeBinning']
        print('Process fastq binning......', flush=True)
        process_status = subprocess.run([f'{script_path}/transfer_barcode -i {fq2} -l {out_path}/config/expect_barcode.tsv -o {out_fq2_name} -d {ham_dist} -p {out_path}/data'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        stdout = process_status.stdout
        stderr = process_status.stderr
        returncode = process_status.returncode
        if int(returncode) != 0:
            raise Exception(f'Bin barcode error with {stderr}')
    else:
        os.symlink(f'{fq2}', f'{out_path}/data/{fq2_name}')
        out_fq2_name = fq2_name

    return fq1_name, out_fq2_name


def modify_yaml(data, fq1_names, fq2_names):
    sample_type = data['sample']['sample_type'].lower()
    out_path=data['out_dir']
    
    data['sequence_files']['file1']['name'] = f'{out_path}/data/{fq1_names}'
    data['sequence_files']['file2']['name'] = f'{out_path}/data/{fq2_names}'

    if sample_type == 'manual':
        data['barcodes']['barcode_file'] = f'{out_path}/config/expect_bin_barcode.tsv'
    else:
        data['barcodes']['barcode_file'] = f'{out_path}/config/expect_barcode.tsv'
    

    new_data=copy.deepcopy(data)
    del new_data['sample']
    new_data['out_dir'] = f'{out_path}/analysis'


    class ForceStr:
        def __init__(self, value):
            self.value = value

    def force_str_representer(dumper, data):
        return dumper.represent_scalar('tag:yaml.org,2002:str', data.value, style='\'')

    yaml.add_representer(ForceStr, force_str_representer)

    new_data['counting_opts']['downsampling']  = ForceStr(new_data['counting_opts']['downsampling'])

    with open(f'{out_path}/config/final_config.yaml', 'w') as out_yaml:
        yaml.dump(new_data, out_yaml, default_flow_style=False, sort_keys=False)

    subprocess.run([f'sed -i "s#: false#: no#g" {out_path}/config/final_config.yaml'], shell=True)
    subprocess.run([f'sed -i "s#: true#: yes#g" {out_path}/config/final_config.yaml'], shell=True)
    subprocess.run([f'sed -i "s#: null#: ~#g" {out_path}/config/final_config.yaml'], shell=True)
    

def run_zUMIs(data):
    script_path = data['zUMIs_directory']
    out_path = data['out_dir']
    
    with open(f'{out_path}/analysis/zUMIs_run.log','w') as run_log, open(f'{out_path}/analysis/zUMIs_run_error.log','w') as error_log:
        print(f'run script: {script_path}/zUMIs.sh -c -y {out_path}/config/final_config.yaml', file=run_log, flush=True)
        process_status = subprocess.run([f'{script_path}/zUMIs.sh -c -y {out_path}/config/final_config.yaml'], stdout=run_log, stderr=error_log, text=True, shell=True)

    returncode = process_status.returncode
    if int(returncode) != 0:
        raise Exception(f'Error happen when running zUMIs, see {out_path}/analysis/zUMIs_run_error.log for detail.')


def run_create_reports(data):
    sample_name = data['project']
    sample_species = data['sample']['sample_species'].upper()
    script_path = data['zUMIs_directory']
    result_dir = data['out_dir']

    if sample_species == 'HUMAN':
        sample_species = 'Human'
    else:
        sample_species = 'Mouse'
    
    with open(f'{result_dir}/analysis/report_run.log','w') as run_log, open(f'{result_dir}/analysis/report_run_error.log','w') as error_log:
        process_status = subprocess.run([f'{script_path}/report/create_summary --sample {sample_name} --indir {result_dir}/analysis --species {sample_species} --well {result_dir}/config/expect_id_barcode.tsv'], stdout=run_log, stderr=error_log, text=True, shell=True)

    returncode = process_status.returncode
    if int(returncode) != 0:
        raise Exception(f'Error happen when cerate reports, see {result_dir}/analysis/report_run_error.log for detail.')
    
    shutil.copytree(f'{result_dir}/analysis/summary/div', f'{result_dir}/results/html')
    shutil.copy(f'{result_dir}/analysis/summary/{sample_name}_stat.xls', f'{result_dir}/results/{sample_name}_stat.xls')
    shutil.copy(f'{result_dir}/analysis/zUMIs_output/stats/{sample_name}.features.pdf', f'{result_dir}/results/{sample_name}.features.pdf')


def get_args():
    parser = argparse.ArgumentParser(description='MGIEasy high sensitive full length transcriptome data analysis pipeline.')
    parser.add_argument('-y', '--yaml', required=True, type=str, help='Path of yaml file.')
    args = parser.parse_args()
    return args


def main():

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start analysis.', flush=True)
    args = get_args()
    
    data=load_yaml(args.yaml)    
    check_nested_empty(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} YAML file load and check complete, all required fields are not empty.', flush=True)

    check_species(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Species is support.', flush=True)

    make_dir(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Create output directory finish.', flush=True)

    check_file_exists(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Fastq, reference and annotation files all exist.', flush=True)

    validate_type_id(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Check sample_type and sample_id are matched successfully.', flush=True)

    create_barcode(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Create barcode file finsh.', flush=True)

    fq1_names, fq2_names = process_fq(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Process fastq finish.', flush=True)

    modify_yaml(data, fq1_names, fq2_names)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Create new yaml finish.', flush=True)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start run zUMIs.', flush=True)
    run_zUMIs(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Finish run zUMIs.', flush=True)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Start create results.', flush=True)
    run_create_reports(data)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Finish create results.', flush=True)

    print('All analysis finish, bye.')

    
if __name__ == '__main__':
    main()
