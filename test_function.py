import unittest
import os.path as osp
from model_library import SV_partition, SV_cell, SV_enumeration, SV_sampling, SV_wssv
import argparse

abspath = osp.abspath(__file__)
directories = abspath.split("/")
for i in range(len(directories) - 2, -1, -1):
    if directories[i] == "PreferenceShapley":
        root_index = i
        break
root_dir = "/".join(directories[:root_index+1])
test3d_fns = ["3d1000n_10p_anti_layer_0","3d1000n_10p_anti_random_0"]
test3d_fps = [osp.join(root_dir,"data/data_for_test/",fn) for fn in test3d_fns]
test2d_fns = ["2d10000n_10p_anti_layer_0","2d10000n_10p_anti_random_0"]
test2d_fps = [osp.join(root_dir,"data/data_for_test/",fn) for fn in test2d_fns]
pssv_exefn = osp.join(root_dir,'main_compute_partition_hd')
cell_exefn = osp.join(root_dir,'main_compute_partition_hd')
wssv_exefn = osp.join(root_dir,'main_compute_skyline')
enum_exefn = osp.join(root_dir,'main_permutation')
perm_exefn = osp.join(root_dir,'main_permutation')

pssv2d_exefn = osp.join(root_dir,'main_compute_sv_2d')
wssv2d_exefn = osp.join(root_dir,'main_compute_sv_2d')
cell2d_exefn = osp.join(root_dir,'main_compute_sv_2d')
enum2d_exefn = osp.join(root_dir,'main_permutation')
perm2d_exefn = osp.join(root_dir,'main_permutation')

class Test_SV_method(unittest.TestCase):

    def test_wssv(self):
        for each_fp in test3d_fps:
            result = SV_wssv(each_fp,each_fp+'output',wssv_exefn,
                                  niter = 5, neval = 1000, num_process=16, preprocessed_flag=False)
            self.assertTrue(len(result['ret']['SV']) == 10)
            self.assertTrue(sum(result['ret']['SV'])>300)
            self.assertTrue(sum([x>20 for x in result['ret']['SV']])==10)
    def test_wssv2d(self):
        for each_fp in test2d_fps:
            result = SV_wssv(each_fp,each_fp+'output',wssv2d_exefn,
                             niter = 5, neval = 1000, num_process=16, preprocessed_flag=False,flag_2d=True)
            self.assertTrue(len(result['ret']['SV']) == 10)
            self.assertTrue(sum(result['ret']['SV'])>3000)
            self.assertTrue(sum([x>300 for x in result['ret']['SV']])==10)
    def test_pssv(self):
        for each_fp in test3d_fps:
            result = SV_partition(each_fp,each_fp+'output',pssv_exefn,method='part_sample_SV',
                                  niter = 5, neval = 10000, num_process=16, preprocessed_flag=False)
            self.assertTrue(len(result['ret']['SV']) == 10)
            self.assertTrue(sum(result['ret']['SV'])>300)
            self.assertTrue(sum([x>20 for x in result['ret']['SV']])==10)
    def test_pssv2d(self):
        for each_fp in test2d_fps:
            result = SV_partition(each_fp,each_fp+'output',pssv2d_exefn,method='part_sample_SV',
                                  niter = 5, neval = 10000, num_process=16, preprocessed_flag=False,flag_2d=True)
            self.assertTrue(len(result['ret']['SV']) == 10)
            self.assertTrue(sum(result['ret']['SV'])>3000)
            self.assertTrue(sum([x>300 for x in result['ret']['SV']])==10)
    def test_cell(self):        
        for each_fp in test3d_fps:
            result = SV_cell(each_fp,each_fp+'output',cell_exefn, niter = 5, 
                             neval = 1000, num_process=16, preprocessed_flag=False)
            self.assertTrue(len(result['ret']['SV']) == 10)
            self.assertTrue(sum(result['ret']['SV'])>300)
            self.assertTrue(sum([x>20 for x in result['ret']['SV']])==10)
    def test_cell2d(self):
        ground_truth = [[420.821, 406.683, 404.811, 407.223, 387.01, 391.067, 392.887, 384.831, 415.438, 375.298],
                   [399.549, 394.307, 407.04, 398.245, 395.202, 398.763, 395.96, 397.542, 402.543, 396.916]]                
        ct = 0
        for each_fp in test2d_fps:
            print(each_fp)
            result = SV_cell(each_fp,each_fp+'output',cell2d_exefn, niter = 5, 
                             neval = 1000, num_process=16, preprocessed_flag=False,flag_2d=True)
            if ct == 0:
                self.assertTrue(sum([i==j for i,j in zip(result['ret']['SV'],ground_truth[0])])==10)
                ct += 1
            elif ct == 1:
                self.assertTrue(sum([i==j for i,j in zip(result['ret']['SV'],ground_truth[1])])==10)
    def test_enum(self):     
        for each_fp in test3d_fps:
            result = SV_enumeration(each_fp,each_fp+'output',enum_exefn, part_flag = False,
                                    niter = 5, neval = 1000, num_process=16, preprocessed_flag=False)
            self.assertTrue(len(result['ret']['SV']) == 10)
            self.assertTrue(sum(result['ret']['SV'])>300)
            self.assertTrue(sum([x>20 for x in result['ret']['SV']])==10)
    def test_enum2d(self):     
        ground_truth = [[420.821, 406.683, 404.811, 407.223, 387.01, 391.067, 392.887, 384.831, 415.438, 375.298],
                   [399.549, 394.307, 407.04, 398.245, 395.202, 398.763, 395.96, 397.542, 402.543, 396.916]]
        ct = 0
        for each_fp in test2d_fps:
            result = SV_enumeration(each_fp,each_fp+'output',enum_exefn, part_flag = False,
                                    niter = 5, neval = 1000, num_process=16, preprocessed_flag=False,flag_2d=True)
            if ct == 0:
                self.assertTrue(sum([i==j for i,j in zip(result['ret']['SV'],ground_truth[0])])==10)
                ct += 1
            elif ct == 1:
                self.assertTrue(sum([i==j for i,j in zip(result['ret']['SV'],ground_truth[1])])==10)
    def test_perm(self):
        for each_fp in test3d_fps:
            result = SV_sampling(each_fp,each_fp+'output',perm_exefn, num_samples=100,
                                    niter = 5, neval = 1000, num_process=16, preprocessed_flag=False)
            self.assertTrue(len(result['ret']['SV']) == 10)
            self.assertTrue(sum(result['ret']['SV'])>300)
            self.assertTrue(sum([x>0 for x in result['ret']['SV']])>5)
    def test_perm2d(self):
        for each_fp in test2d_fps:
            result = SV_sampling(each_fp,each_fp+'output',perm_exefn, num_samples=100,
                                    niter = 5, neval = 1000, num_process=16, preprocessed_flag=False,flag_2d=True)
            self.assertTrue(len(result['ret']['SV']) == 10)
            self.assertTrue(sum(result['ret']['SV'])>3000)
            self.assertTrue(sum([x>0 for x in result['ret']['SV']])>5)

if __name__ == '__main__':
    # unittest.main()

    parser = argparse.ArgumentParser(description='Test SV calculation methods')
    parser.add_argument('-t', '--models_type', choices=['all','pssv','enum','perm','cell','wssv','cell2d','wssv2d','pssv2d','enum2d','perm2d'], default='all', help='Model type selection (all, or pssv, cell, enum, enpt and perm)')
    args = parser.parse_args()

    if args.models_type == 'pssv':
        suite = unittest.TestLoader().loadTestsFromName('test_pssv', Test_SV_method)
    if args.models_type == 'pssv2d':
        suite = unittest.TestLoader().loadTestsFromName('test_pssv2d', Test_SV_method)
    elif args.models_type == 'wssv':
        suite = unittest.TestLoader().loadTestsFromName('test_wssv', Test_SV_method)
    elif args.models_type == 'wssv2d':
        suite = unittest.TestLoader().loadTestsFromName('test_wssv2d', Test_SV_method)
    elif args.models_type == 'cell':
        suite = unittest.TestLoader().loadTestsFromName('test_cell', Test_SV_method)
    elif args.models_type == 'cell2d':
        suite = unittest.TestLoader().loadTestsFromName('test_cell2d', Test_SV_method)
    elif args.models_type == 'enum':
        suite = unittest.TestLoader().loadTestsFromName('test_enum', Test_SV_method)
    elif args.models_type == 'enum2d':
        suite = unittest.TestLoader().loadTestsFromName('test_enum2d', Test_SV_method)
    elif args.models_type == 'perm':
        suite = unittest.TestLoader().loadTestsFromName('test_perm', Test_SV_method)
    elif args.models_type == 'perm2d':
        suite = unittest.TestLoader().loadTestsFromName('test_perm2d', Test_SV_method)
    elif args.models_type == 'all':
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_SV_method)

    # Run the tests
    unittest.TextTestRunner().run(suite)
