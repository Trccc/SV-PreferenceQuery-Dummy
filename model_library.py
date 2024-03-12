import os
import time
from functools import partial
from hd_SV_mp import calculate_SV_HD_partition_mp
from hd_perm_sample import calculate_SV_HD_permutation_sampling_mp
from hd_enum_coalition_mp import calculate_SV_HD_coalition_enumeration_mp
from hd_SV_cell_mp import calculate_SV_HD_cell_mp
from hyperplane_arrangement import wssv
import numpy as np
import pickle
from utility import partition_parser,cell_parser,permutation_sampling_parser,enumerate_coalitions_parser,enumerate_coalitions_part_parser,parser_raw_data

def find_sv_in_output(fn,baseline_flag = False):
    with open(fn, 'r') as file:
        lines = file.readlines()
    if not baseline_flag:
        last_lines = lines[-3]
    else:
        last_lines = lines[-4]
    return list(map(float,last_lines.strip().split(' ')))


def SV_partition(input_fn,output_fn,partition_exe_fn,topk=1,save_flag=False,method='part_sample_SV',niter = 10, neval = 1000, num_process=2, preprocessed_flag=False, flag_2d = False):
    result = {}
    if not flag_2d:
        if not preprocessed_flag:
            if topk == 1:
                command = f"{partition_exe_fn} -i {input_fn} -o {output_fn} -n -dc -p {num_process} -pp {num_process}"
            else:
                command = f"{partition_exe_fn} -i {input_fn} -o {output_fn} -n -k {topk} -dc -p {num_process}"
            t0 = time.time()
            os.system(command)
            t1 = time.time()
            result['time_for_partition'] = t1-t0
            ret = calculate_SV_HD_partition_mp(output_fn,topk=topk,save_flag=save_flag,method=method,niter = niter, neval = neval, preprocessed_flag=preprocessed_flag)
            result['ret'] = ret
            return result  
        else:
            parsed_data = pickle.load(open(input_fn,'rb'))
            result['time_for_preprocessing'] = parsed_data['time_for_preprocessing']
            ret = calculate_SV_HD_partition_mp(parsed_data['parsed_data'],topk=topk,save_flag=save_flag,method=method,niter = niter, neval = neval, preprocessed_flag=preprocessed_flag)
            result['ret'] = ret
            return result 
    else:
        command = f"{partition_exe_fn} -i {input_fn} -f -k {topk} -p {num_process} -pssv {niter*neval}"
        t0 = time.time()
        os.system(command)
        t1 = time.time()
        result['time_for_computation'] = t1-t0
        output_fn = input_fn + '_output'
        result['ret'] = {'SV':find_sv_in_output(output_fn)}
        if not save_flag:
            os.remove(output_fn)
        return result   



def SV_cell(input_fn,output_fn,arr_exe_fn,topk=1,num_process = 2, niter = 10, neval = 1000, save_flag=False, preprocessed_flag=False,flag_2d=False):
    result = {}
    if not flag_2d:
        if not preprocessed_flag:
            if topk == 1:
                command = f"{arr_exe_fn} -i {input_fn} -o {output_fn} -n -dc -a -ao {output_fn} -p {num_process} -pp {num_process}"
            else:
                command = f"{arr_exe_fn} -i {input_fn} -o {output_fn} -n -k {topk} -dc -a -ao {output_fn} -p {num_process}"
            t0 = time.time()
            os.system(command)
            t1 = time.time()
            result['time_for_partition_arr'] = t1-t0
            output_dir_fn = []
            if os.path.exists(output_fn+f'_{1}_{0}'):
                output_dir_fn.append(output_fn+f'_{1}_{0}')
            else:
                for i in range(num_process):
                    can_fn = output_fn+f'_{num_process}_{i}'
                    if os.path.exists(can_fn):
                        output_dir_fn.append(output_fn+f'_{num_process}_{i}')
                        if not save_flag:
                            os.remove(output_fn+f'_{i}')
            assert (len(output_dir_fn)==num_process) or (len(output_dir_fn)==1)
            tmp_SV = []
            for idx,each_fn in enumerate(output_dir_fn):
                ret = calculate_SV_HD_cell_mp(each_fn,topk=topk,save_flag=save_flag,niter = niter,neval=neval,preprocessed_flag=preprocessed_flag)
                tmp_SV.append(ret['SV'])
            sum_SV = np.sum(tmp_SV,axis=0)
            result['ret'] = {'SV':sum_SV.tolist()}
            return result   
        else:
            parsed_data = pickle.load(open(input_fn,'rb'))
            result['time_for_preprocessing'] = parsed_data['time_for_preprocessing']
            tmp_SV = []
            for idx,each_parsed in enumerate(parsed_data['parsed_data']):
                ret = calculate_SV_HD_cell_mp(each_parsed,topk=topk,save_flag=save_flag,niter = niter, neval = neval, preprocessed_flag=preprocessed_flag)
                tmp_SV.append(ret['SV'])
            sum_SV = np.sum(tmp_SV,axis=0)
            result['ret'] = {'SV':sum_SV.tolist()}
            return result   
    else:
        command = f"{arr_exe_fn} -i {input_fn} -f -k {topk} -p {num_process}"
        t0 = time.time()
        os.system(command)
        t1 = time.time()
        result['time_for_computation'] = t1-t0
        output_fn = input_fn + '_output'
        result['ret'] = {'SV':find_sv_in_output(output_fn)}
        if not save_flag:
            os.remove(output_fn)
        return result   

def SV_enumeration(input_fn,output_fn,enumeration_exe_fn,topk=1,num_process = 2, save_flag=False, part_flag = True, niter = 10, neval = 1000,preprocessed_flag=False,flag_2d=False):
    result = {}
    if not flag_2d:
        if not preprocessed_flag:
            if part_flag == True:
                command = f"{enumeration_exe_fn} -i {input_fn} -o {output_fn} -k {topk} -n -dc -p {num_process} -part"
            elif part_flag == False:
                command = f"{enumeration_exe_fn} -i {input_fn} -o {output_fn} -k {topk} -n -dc -p {num_process}"
            t0 = time.time()
            os.system(command)
            t1 = time.time()
            result['time_for_partition'] = t1-t0
            ret = calculate_SV_HD_coalition_enumeration_mp(output_fn,topk=topk,save_flag=save_flag,part_flag = part_flag,niter = niter,neval=neval,preprocessed_flag=preprocessed_flag)
            result['ret'] = {'SV':ret}
            return result
        else:
            parsed_data = pickle.load(open(input_fn,'rb'))
            result['time_for_preprocessing'] = parsed_data['time_for_preprocessing']
            ret = calculate_SV_HD_coalition_enumeration_mp(parsed_data['parsed_data'],topk=topk,save_flag=save_flag,part_flag = part_flag,niter = niter, neval = neval, preprocessed_flag=preprocessed_flag)
            result['ret'] = {'SV':ret}
            return result   
    else:
        command = f"{enumeration_exe_fn} -i {input_fn} -f -k {topk} -p {num_process}"
        t0 = time.time()
        os.system(command)
        t1 = time.time()
        result['time_for_computation'] = t1-t0
        output_fn = input_fn + '_output'
        print(input_fn)
        print(output_fn)
        result['ret'] = {'SV':find_sv_in_output(output_fn,True)}
        if not save_flag:
            os.remove(output_fn)
        return result   



def SV_sampling(input_fn,output_fn,sampling_exe_fn,topk=1, num_samples = 16, num_process = 2,save_flag=False,niter = 10,neval = 1000,preprocessed_flag=False,flag_2d = False):
    result = {}
    if not flag_2d:
        if not preprocessed_flag:
            if topk == 1:
                command = f"{sampling_exe_fn} -i {input_fn} -o {output_fn} -n -s {num_samples} -dc -p {num_process}"
            else:
                command = f"{sampling_exe_fn} -i {input_fn} -o {output_fn} -n -s {num_samples} -k {topk} -dc -p {num_process}"
            t0 = time.time()
            os.system(command)
            t1 = time.time()
            result['time_for_partition'] = t1-t0
            ret = calculate_SV_HD_permutation_sampling_mp(output_fn,topk=topk,save_flag=save_flag,niter = niter,neval = neval,preprocessed_flag=preprocessed_flag)
            result['ret'] = {'SV':ret}
            result['num_samples'] = num_samples
            return result
        else:
            parsed_data = pickle.load(open(input_fn,'rb'))
            result['time_for_preprocessing'] = parsed_data['time_for_preprocessing']
            ret = calculate_SV_HD_permutation_sampling_mp(parsed_data['parsed_data'],topk=topk,save_flag=save_flag,niter = niter,neval = neval,preprocessed_flag=preprocessed_flag)
            result['ret'] = {'SV':ret}
            return result   
    else:
        command = f"{sampling_exe_fn} -i {input_fn} -f -k {topk} -p {num_process} -s {num_samples}"
        t0 = time.time()
        os.system(command)
        t1 = time.time()
        output_fn = input_fn + '_output'
        result['ret'] = {'SV':find_sv_in_output(output_fn,True)}
        result['num_samples'] = num_samples
        if not save_flag:
            os.remove(output_fn)
        return result   
    

def SV_wssv(input_fn,output_fn,wssv_exe_fn,topk=1, num_process = 2,save_flag=False,niter = 10,neval = 1000,preprocessed_flag=False,flag_2d = False):
    result = {}
    if not flag_2d:
        if not preprocessed_flag:
            command = f"{wssv_exe_fn} -i {input_fn} -o {output_fn} -dc -p {num_process}"
            t0 = time.time()
            os.system(command)
            t1 = time.time()
            result['time_for_partition'] = t1-t0
            parsed_data = parser_raw_data(output_fn,has_index=True,save_flag=save_flag)
            n_owners = parsed_data['num_data_owners']
            dimension = parsed_data['dimension']
            num_points = []
            for i in range(n_owners):
                num_points.append(len(parsed_data['points'][i]['pvals']))
            max_num_points = max(num_points)
            padding_points = [0]*dimension
            combined_pt = np.concatenate(([[parsed_data['points'][i]['pvals']+[padding_points]*(max_num_points-len(parsed_data['points'][i]['pvals']))] for i in range(n_owners)]),axis=0)
            my_wssv = wssv(combined_pt, dimension, n_owners)
            ret = my_wssv.calculate_SV_whole_space(niter,neval)
            result['ret'] = {'SV':ret}
            return result
        else:
            loaded_data = pickle.load(open(input_fn,'rb'))
            result['time_for_preprocessing'] = loaded_data['time_for_preprocessing']
            parsed_data = loaded_data['parsed_data']
            # parsed_data = parser_raw_data(input_fn,save_flag=save_flag)
            n_owners = parsed_data['num_data_owners']
            dimension = parsed_data['dimension']
            num_points = []
            for i in range(n_owners):
                num_points.append(len(parsed_data['points'][i]['pvals']))
            # print(f'num_points is {num_points}')
            max_num_points = max(num_points)
            # print(f'max_num_points is {max_num_points}')
            padding_points = [0]*dimension
            combined_pt = []
            # for i in range(n_owners):
            #     combined_pt.append(parsed_data['points'][i]['pvals']+[padding_points]*(max_num_points-len(parsed_data['points'][i]['pvals'])))
            # print(combined_pt)    
            combined_pt = np.concatenate(([[parsed_data['points'][i]['pvals']+[padding_points]*(max_num_points-len(parsed_data['points'][i]['pvals']))] for i in range(n_owners)]),axis=0)
            my_wssv = wssv(combined_pt, dimension, n_owners)
            ret = my_wssv.calculate_SV_whole_space(niter,neval)
            result['ret'] = {'SV':ret}
            return result
    else:
        command = f"{wssv_exe_fn} -i {input_fn} -f -k {topk} -p {num_process} -wssv {niter*neval}"
        t0 = time.time()
        os.system(command)
        t1 = time.time()
        result['time_for_computation'] = t1-t0
        output_fn = input_fn + '_output'
        result['ret'] = {'SV':find_sv_in_output(output_fn)}
        if not save_flag:
            os.remove(output_fn)
        return result   


def pre_SV(input_fn,exe_fn,method,save_flag=False, num_process=2, num_samples=0, overwrite_flag=False):
    if num_samples:
        assert method == 'perm'
    result = {}
    output_fn = input_fn+'_output'
    if num_samples:
        method = f'permX{num_samples}'
    save_fn = input_fn+f'_{method}_parsed.pk'
    if os.path.exists(save_fn) and not overwrite_flag:
        print(f'parsed data already exists')
        return
    else:
        if method == 'pssv':
            command = f"{exe_fn} -i {input_fn} -o {output_fn} -n -dc -p {num_process} -pp {num_process}"
            parser = partition_parser
        elif method == 'enum':
            command = f"{exe_fn} -i {input_fn} -o {output_fn} -n -dc -p {num_process}"
            parser = enumerate_coalitions_parser
        elif method == 'enpt':
            command = f"{exe_fn} -i {input_fn} -o {output_fn} -n -dc -p {num_process} -part"
            parser = enumerate_coalitions_part_parser
        elif method.startswith('perm'):
            assert num_samples is not None
            command = f"{exe_fn} -i {input_fn} -o {output_fn} -n -s {num_samples} -dc -p {num_process}"
            parser = permutation_sampling_parser
        elif method == 'cell':
            command = f"{exe_fn} -i {input_fn} -o {output_fn} -n -dc -a -ao {output_fn} -p {num_process} -pp {num_process}"
            parser = cell_parser
        elif method == 'wssv':
            command = f"{exe_fn} -i {input_fn} -o {output_fn} -dc -p {num_process}"
            parser = partial(parser_raw_data,has_index=True)
        t0 = time.time()
        os.system(command)
        t1 = time.time()
        if method != 'cell':
            parsed_data = parser(output_fn,save_flag=save_flag)
            result['time_for_preprocessing'] = t1-t0
            result['parsed_data'] = parsed_data
            pickle.dump(result,open(save_fn,'wb'))
            print(f'{method} data saved')
        else:
            output_dir_fn = []
            if os.path.exists(output_fn+f'_{1}_{0}'):
                # print(output_fn+f'_{1}_{0}')
                # print(f"np = {num_process}")
                output_dir_fn.append(output_fn+f'_{1}_{0}')
                if not save_flag:
                    os.remove(output_fn+f'_{0}')

            else:
                # print(f"np2 = {num_process}")
                for i in range(num_process):
                    can_fn = output_fn+f'_{num_process}_{i}'
                    if os.path.exists(can_fn):
                        output_dir_fn.append(output_fn+f'_{num_process}_{i}')
                        if not save_flag:
                            os.remove(output_fn+f'_{i}')
            assert (len(output_dir_fn)==num_process) or (len(output_dir_fn)==1)
            tmp_SV = []
            for idx,each_fn in enumerate(output_dir_fn):
                tmp_SV.append(parser(each_fn,save_flag=save_flag))
            result['time_for_preprocessing'] = t1-t0
            result['parsed_data'] = tmp_SV
            pickle.dump(result,open(save_fn,'wb'))
            print(f'{method} data saved')
 