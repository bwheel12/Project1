#initialize alignment classes
%run ./Project1-main/align/algs.py
#read in negative pairs
with open('./Project1-main/'+ 'scoring_matrices/Negpairs.txt') as f:
    pairs_txt = f.readlines()
    pairs = []
for i in range(0,len(pairs_txt)):
    pair_temp = [pairs_txt[i][0:pairs_txt[i].find(" ")],pairs_txt[i][pairs_txt[i].find(" ")+1:pairs_txt[i].find("\n")]]
    pairs.append(pair_temp)
    
#run local alignment on negative pairs
neg_align_seq = []
neg_align_scr = []
for j in range(len(pairs)):
        temp_align = SmithWaterman(pairs[j][0],pairs[j][1],'./Project1-main/scoring_matrices/BLOSUM50.mat')
        temp_align.read_scoring_mat()
        temp_align.set_gap_penalties(-11,-3)
        temp_align.set_up_align_mats()
        temp_align.step_through()
        temp_align.follow_back()
        neg_align_seq.append([temp_align.final_sequence_1,temp_align.final_sequence_2])
        neg_align_scr.append(temp_align.alignment_score)

        
#read in positive pairs
with open('./Project1-main/'+ 'scoring_matrices/Pospairs.txt') as f:
    pairs_txt = f.readlines()
    pairs = []
for i in range(0,len(pairs_txt)):
    pair_temp = [pairs_txt[i][0:pairs_txt[i].find(" ")],pairs_txt[i][pairs_txt[i].find(" ")+1:pairs_txt[i].find("\n")]]
    pairs.append(pair_temp)
#align positive pairs
pos_align_seq = []
pos_align_scr = []
for j in range(len(pairs)):
        temp_align = SmithWaterman(pairs[j][0],pairs[j][1],'./Project1-main/scoring_matrices/BLOSUM50.mat')
        temp_align.read_scoring_mat()
        temp_align.set_gap_penalties(-11,-3)
        temp_align.set_up_align_mats()
        temp_align.step_through()
        temp_align.follow_back()
        pos_align_seq.append([temp_align.final_sequence_1,temp_align.final_sequence_2])
        pos_align_scr.append(temp_align.alignment_score)
        
numpy.mean(neg_align_scr)
numpy.mean(pos_align_scr)
        
import matplotlib.pyplot as plt

bins = numpy.linspace(0, 400, 100)
plt.hist(neg_align_scr, bins, alpha = 0.5, label = 'negative')
plt.hist(pos_align_scr, bins, alpha = 0.5, label = 'positive')
plt.legend(loc='upper right')
plt.show()

all_scores = neg_align_scr + pos_align_scr
all_score_mean = numpy.mean(all_scores)

confusion = numpy.zeros([2,2])
true_neg = len([i for i in neg_align_scr if i < all_score_mean])
confusion[0,0] = true_neg
false_pos = len([i for i in neg_align_scr if i >= all_score_mean])
confusion[0,1] = false_pos
false_neg = len([i for i in pos_align_scr if i < all_score_mean])
confusion[1,0] = false_neg
true_pos = len([i for i in pos_align_scr if i >= all_score_mean])
confusion[1,1] = true_pos

TPR = true_pos/(false_neg+true_pos)
FPR = false_pos/(false_pos+true_neg)

roc_TPRs = []
roc_FPRs = []

for i in range(400):
    all_score_mean = i
    confusion = numpy.zeros([2,2])
    true_neg = len([i for i in neg_align_scr if i < all_score_mean])
    confusion[0,0] = true_neg
    false_pos = len([i for i in neg_align_scr if i >= all_score_mean])
    confusion[0,1] = false_pos
    false_neg = len([i for i in pos_align_scr if i < all_score_mean])
    confusion[1,0] = false_neg
    true_pos = len([i for i in pos_align_scr if i >= all_score_mean])
    confusion[1,1] = true_pos

    TPR = true_pos/(false_neg+true_pos)
    FPR = false_pos/(false_pos+true_neg)
    roc_TPRs.append(TPR)
    roc_FPRs.append(FPR)
    
plt.plot(roc_FPRs, roc_TPRs, 'o', color = 'black')

#determine ROC_AUC
roc_AUC = 0

for i in range(len(roc_TPRs)-1):
    roc_AUC = roc_AUC + (roc_TPRs[i]+roc_TPRs[i+1])/2*(roc_FPRs[i]-roc_FPRs[i+1])
    
    
#