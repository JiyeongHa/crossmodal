def load_attributes(attr_file):
    x = os.path.join(attr_file)
    attr = ColumnData(x, header=True)
    # attr = SampleAttributes(x)
    return attr
def load_nii(nii_file, mask_file, attr):
    """load experiment dataset"""

    fds = fmri_dataset(samples=os.path.join(nii_file),
                       targets=attr.xMD, chunks=attr.run,
                       mask=os.path.join(mask_file))
    return fds
def lag_correction(fds, runTRs, lagTRs):
    """correct dataset for hemodynamic lag"""

    # split dataset into runs

    nRuns = len(fds) / float(runTRs)
    if int(nRuns) != nRuns:
        print 'Error! number of TRs per run must be a factor of total TRs'
        raise SystemExit

    nRuns = int(nRuns)

    split_fds = []

    for i in range(nRuns):  # split dataset into separate runs
        split_fds.append(fds[i * runTRs:(i + 1) * runTRs])

    # do the shift for each run

    for i in range(len(split_fds)):
        split_fds[i].sa.targets[lagTRs:] = \
            split_fds[i].sa.targets[:-lagTRs]  # need to shift target labels too

        split_fds[i].sa.censor[lagTRs:] = \
            (split_fds[i].sa.censor[:-lagTRs])  # and censor labels

        split_fds[i].sa.trial[lagTRs:] = \
            split_fds[i].sa.trial[:-lagTRs]  # and trial label

        split_fds[i].sa.chunks[lagTRs:] = \
            split_fds[i].sa.chunks[:-lagTRs]  # and run label

        split_fds[i].sa.Odd[lagTRs:] = \
            split_fds[i].sa.Odd[:-lagTRs]  # and run labe

        split_fds[i].sa.cond[lagTRs:] = \
            split_fds[i].sa.cond[:-lagTRs]  # and run label

        split_fds[i].sa.TR[lagTRs:] = \
            split_fds[i].sa.TR[:-lagTRs]  # and run label

        split_fds[i].sa.Alabel[lagTRs:] = \
            split_fds[i].sa.Alabel[:-lagTRs]  # and run label


        split_fds[i] = (split_fds[i])[lagTRs:]

    ##  merge back datasets
    fds = split_fds[0]
    for i in range(1, len(split_fds)):
        fds.append(split_fds[i])

    return fds
def make_null_dist_plot(dist_samples, empirical):
    pl.hist(dist_samples, bins=20, normed=True, alpha=0.8)
    pl.axvline(empirical, color='red')
    # chance-level for a binary classification with balanced samples
    pl.axvline(0.5, color='black', ls='--')
    # scale x-axis to full range of possible error values
    pl.xlim(0,1)
    pl.xlabel('Average cross-validated classification error')

    _ = pl.figure()

if __debug__:
    from mvpa2.base import debug
# libraries needed by pymvpa
import os
import sys
from mvpa2.suite import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as ss

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# subject info.
#sName = ["01HJH", "02YYH", "03LJA", "04LSY", "05YJG", "06JDH", "07KHY", "08LKY_new", "09JEH", "10NJS", "11KDB", "12SCL", "13NHJ", "14JJY"]
sName = ["01HJH",  "03LJA",  "06JDH", "07KHY", "08LKY_new", "09JEH", "11KDB", "12SCL", "13NHJ", "14JJY"]
sName = ["01HJH", "02YYH", "03LJA", "04LSY", "05YJG", "06JDH", "07KHY", "08LKY_new", "09JEH", "10NJS", "11KDB", "12SCL", "13NHJ", "14JJY"]
sName = ["01HJH", "02YYH", "03LJA", "04LSY", "06JDH", "07KHY", "08LKY_new", "09JEH", "10NJS", "11KDB", "12SCL", "13NHJ", "14JJY"]

nsbj = len(sName)

# experiment info.
cond = ["cong", "incong", "amb"] #, "incong", "amb"
nRun = 10
nTR = 149
# decoding parameters
useTR = 3
lagTR = 2
nTarget = 2  # the number of testing lines
nCond=0

ROI=[ "A1_fmasked", "LBelt_Complex_fmasked", "MBelt_Complex_fmasked", "PBelt_Complex_fmasked",
      "RetroInsular_Cortex_fmasked", "Area_TA2_fmasked", "Auditory_4_Complex_fmasked",
       "Auditory_5_Complex_fmasked"]

roi_names= 'test'
#bootstrapping parameters










nPerm = 10000 #the number of permutation for each subject
n_null_sample = 0

#result array
group_null_dist = np.zeros((nsbj, nPerm))
group_acc = np.zeros((nsbj,))
group_chunk = np.array(np.arange(1,nsbj+1)).reshape(nsbj,)
group_roi_results = np.zeros((nsbj,len(ROI)))
group_pval = np.zeros((nsbj,len(ROI)))
group_shuf_null_dist = np.zeros((nsbj, nPerm))
tmp_mean_acc_null_dist = np.zeros((1,nPerm))
mean_acc_null_dist = np.zeros((1,nPerm))
ROIs_perm_dist = np.array([]).reshape((0, nPerm+1)) #perm. + real acc.

resultDir = '/media/duri/Decoding_CrM/stat'
if not os.path.exists(resultDir):
    os.makedirs(resultDir)

# ===========load stimulus files==============
voxelNum = [1]
xVoxelNum = 0

for xCond in range(1, len(cond)):
    for xROI in range(0, len(ROI)):
        for xSN in range(0, len(sName)):

            # set directories
            basedir = '/home/haji/Desktop/Decoding_CrMExp2/Decoding_' + sName[xSN] + '/'
            os.chdir(basedir)  # current directory: basedir
            attrDir = os.path.abspath('onset')
            #roiDir = os.path.abspath('ROImasks')
            roiDir = os.path.abspath('ROImasks/audROIs_HCPMMP')

            # file Names in locations
            attrFile = attrDir + '/' + sName[xSN] + '_basicOnset_CIN.txt'  # onset file (= attribute file)
            datFile = basedir + sName[xSN] + 'CrM_CIN_FNL.nii'  # fMRI data file
            # result file name
            xSN_result = np.ones((len(cond), len(ROI))) * 111

            # load stimulus files
            attr = load_attributes(attr_file=attrFile)

            fds = load_nii(nii_file=datFile, mask_file=(roiDir + '/' + ROI[xROI] + '.nii'), attr=attr)
            fds.sa['censor'] = attr.censor
            fds.sa['cond'] = attr.xCIN
            fds.sa['TR'] = attr.TR
            fds.sa['trial'] = attr.trial
            fds.sa['Odd'] = attr.Odd
            fds.sa['Alabel'] = attr.xaMD

            ####### for python 2.7.# ########
            fds.samples = asarray(fds.samples)
            fds = lag_correction(fds=fds, runTRs=nTR, lagTRs=lagTR)  # another custom subfunction
            afterlagN = len(fds)
            print "After lag correction: %d " % afterlagN

            ## remove censored points (motion and outlier)
            fds = fds[fds.sa.censor == 1]
            print "Censored points: %d " % (afterlagN - len(fds))

            ## remove oddball trials
            fds = fds[fds.sa.Odd != 1]

            fds_cond = fds[fds.sa.cond == (xCond + 1)]
            print cond[xCond]
            fds_cond = fds_cond[fds_cond.sa.targets != 0]
            nVox = fds_cond.nfeatures
            ## zscore before removing rest TRs
            zscore(fds_cond, chunks_attr='chunks')
            # fds_cond = fds_cond[np.logical_or(fds_cond.sa.TR == 7*(fds_cond.sa.trial-1)+3, fds_cond.sa.TR == 7*(fds_cond.sa.trial-1)+4)]
            print len(fds_cond)
            ## get a dataset with one sample per stimulus category for each run
            averager = mean_group_sample(['targets', 'chunks'])
            fds_cond = fds_cond.get_mapped(averager)

            ###Haxby et al.
            # clf = kNN(k=5)
            clf = LinearCSVMC()  # kNN(k=1, dfx=one_minus_correlation, voting='majority')
            nfeatures = ceil(nVox * voxelNum[xVoxelNum])
            fsel = SensitivityBasedFeatureSelection(OneWayAnova(),
                                                    FixedNElementTailSelector(nfeatures, tail='upper',
                                                                              mode='select', sort='False'))

            fclf = FeatureSelectionClassifier(clf, fsel)
            permutator = AttributePermutator('targets', count=nPerm)
            distr_est = MCNullDist(permutator, tail='right', enable_ca=['dist_samples'])
            partitioner = ChainNode([NFoldPartitioner(cvtype=1),
                                     Balancer(attr='targets',
                                              count=1,  # for real data > 1
                                              limit='partitions',
                                              apply_selection=True
                                              )],
                                    space='partitions')
            cv_mc_corr = CrossValidation(fclf,
                                         partitioner,
                                         errorfx=mean_match_accuracy,
                                         null_dist=distr_est,
                                         postproc=mean_sample(),
                                         enable_ca=['stats'])
            print "%s start %d permutation testing.... " % (sName[xSN], nPerm)
            acc_result = cv_mc_corr(fds_cond)
            # p-values for each participants
            p = cv_mc_corr.ca.null_prob
            group_pval[xSN,xROI] = np.asscalar(p)
            group_null_dist[xSN,0:nPerm] = np.asarray(cv_mc_corr.null_dist.ca.dist_samples)
            group_acc[xSN] = np.mean(acc_result)
            group_roi_results[xSN][xROI] = np.mean(acc_result)
            print "%s: %s %s permutation has finished. p-value: %2.3f, acc: %2.3f " % (sName[xSN], cond[xCond], ROI[xROI], np.asscalar(p), np.asscalar(acc_result))

        #p-values for mean accuracy
        for i in range(0, nsbj):
            group_shuf_null_dist[i,:] = np.random.permutation(group_null_dist[i,:])

        all_dist = np.concatenate((group_shuf_null_dist, group_acc.reshape(-1,1)), axis=1)
        ##wrap up this ROI result..
        # result file name
        fname_mean_acc_null_dist = '%(path)s/%(cond)s_%(roi)s_acc_all_dist.csv' % {"path": resultDir,
                                                                                   "cond": cond[xCond],
                                                                                   "roi": ROI[xROI]}
        npy_mean_acc_null_dist = '%(path)s/%(cond)s_%(roi)s_acc_all_dist.npy' % {"path": resultDir,
                                                                                 "cond": cond[xCond],
                                                                                 "roi": ROI[xROI]}
        np.savetxt(fname_mean_acc_null_dist, all_dist, delimiter=',') #x=roi, y=perm results
        np.save(npy_mean_acc_null_dist, all_dist)  # x=roi, y=perm results
