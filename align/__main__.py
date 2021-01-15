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
plt.savefig('score_distributions.png')

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

import math
all_scr = neg_align_scr + pos_align_scr
roc_TPRs = []
roc_FPRs = []

for i in range(math.floor(min(all_scr))-1,math.floor(max(all_scr))+2):
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
plt.savefig('ROC.png')

#determine ROC_AUC
roc_AUC = 0

for i in range(len(roc_TPRs)-1):
    roc_AUC = roc_AUC + (roc_TPRs[i]+roc_TPRs[i+1])/2*(roc_FPRs[i]-roc_FPRs[i+1])
    
    
###optimize gap

gap_roc_AUCs = numpy.zeros([20,5])
with open('./Project1-main/'+ 'scoring_matrices/Negpairs.txt') as f:
    pairs_txt_neg = f.readlines()
    pairs_neg = []
for i in range(0,len(pairs_txt_neg)):
    pair_temp_neg = [pairs_txt_neg[i][0:pairs_txt_neg[i].find(" ")],pairs_txt_neg[i][pairs_txt_neg[i].find(" ")+1:pairs_txt_neg[i].find("\n")]]
    pairs_neg.append(pair_temp_neg)

#read in positive pairs
with open('./Project1-main/'+ 'scoring_matrices/Pospairs.txt') as f:
    pairs_txt_pos = f.readlines()
    pairs_pos = []
for i in range(0,len(pairs_txt_pos)):
    pair_temp_pos = [pairs_txt_pos[i][0:pairs_txt_pos[i].find(" ")],pairs_txt_pos[i][pairs_txt_pos[i].find(" ")+1:pairs_txt_pos[i].find("\n")]]
    pairs_pos.append(pair_temp_pos)
            
for z in range(1,21):
    for zed in range(1,6):
        #run local alignment on negative pairs
        neg_align_seq = []
        neg_align_scr = []
        for j in range(len(pairs_neg)):
                temp_align = SmithWaterman(pairs_neg[j][0],pairs_neg[j][1],'./Project1-main/scoring_matrices/BLOSUM50.mat')
                temp_align.read_scoring_mat()
                temp_align.set_gap_penalties(-z,-zed)
                temp_align.set_up_align_mats()
                temp_align.step_through()
                temp_align.follow_back()
                neg_align_seq.append([temp_align.final_sequence_1,temp_align.final_sequence_2])
                neg_align_scr.append(temp_align.alignment_score)


        
        #align positive pairs
        pos_align_seq = []
        pos_align_scr = []
        for j in range(len(pairs_pos)):
                temp_align = SmithWaterman(pairs_pos[j][0],pairs_pos[j][1],'./Project1-main/scoring_matrices/BLOSUM50.mat')
                temp_align.read_scoring_mat()
                temp_align.set_gap_penalties(-z,-zed)
                temp_align.set_up_align_mats()
                temp_align.step_through()
                temp_align.follow_back()
                pos_align_seq.append([temp_align.final_sequence_1,temp_align.final_sequence_2])
                pos_align_scr.append(temp_align.alignment_score)
                
        import math
        all_scr = neg_align_scr + pos_align_scr
        roc_TPRs = []
        roc_FPRs = []

        for i in range(math.floor(min(all_scr))-1,math.floor(max(all_scr))+2):
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
        
        print(z,zed,end='\r')

        #determine ROC_AUC
        roc_AUC = 0

        for i in range(len(roc_TPRs)-1):
            roc_AUC = roc_AUC + (roc_TPRs[i]+roc_TPRs[i+1])/2*(roc_FPRs[i]-roc_FPRs[i+1])
        
        gap_roc_AUCs[z-1,zed-1] = roc_AUC
        
        
        
##Evaluate the different alignment scoring matrices
score_roc_AUCs = numpy.zeros([4])
with open('./Project1-main/'+ 'scoring_matrices/Negpairs.txt') as f:
    pairs_txt_neg = f.readlines()
    pairs_neg = []
for i in range(0,len(pairs_txt_neg)):
    pair_temp_neg = [pairs_txt_neg[i][0:pairs_txt_neg[i].find(" ")],pairs_txt_neg[i][pairs_txt_neg[i].find(" ")+1:pairs_txt_neg[i].find("\n")]]
    pairs_neg.append(pair_temp_neg)

#read in positive pairs
with open('./Project1-main/'+ 'scoring_matrices/Pospairs.txt') as f:
    pairs_txt_pos = f.readlines()
    pairs_pos = []
for i in range(0,len(pairs_txt_pos)):
    pair_temp_pos = [pairs_txt_pos[i][0:pairs_txt_pos[i].find(" ")],pairs_txt_pos[i][pairs_txt_pos[i].find(" ")+1:pairs_txt_pos[i].find("\n")]]
    pairs_pos.append(pair_temp_pos)
    
test_matrices = ['BLOSUM50.mat','BLOSUM62.mat','PAM100.mat','PAM250.mat']


for zz in range(4):
    matrix_location = './Project1-main/scoring_matrices/' + test_matrices[zz]
    z = 6 #the gap opening penalty
    zed = 5 #the gap extension penalty
    #run local alignment on negative pairs
    neg_align_seq = []
    neg_align_scr = []
    for j in range(len(pairs_neg)):
            temp_align = NeedlemanWunsch(pairs_neg[j][0],pairs_neg[j][1],matrix_location)
            temp_align.read_scoring_mat()
            temp_align.set_gap_penalties(-z,-zed)
            temp_align.set_up_align_mats()
            temp_align.step_through()
            temp_align.follow_back()
            neg_align_seq.append([temp_align.final_sequence_1,temp_align.final_sequence_2])
            neg_align_scr.append(temp_align.alignment_score)



    #align positive pairs
    pos_align_seq = []
    pos_align_scr = []
    for j in range(len(pairs_pos)):
            temp_align = NeedlemanWunsch(pairs_pos[j][0],pairs_pos[j][1],matrix_location)
            temp_align.read_scoring_mat()
            temp_align.set_gap_penalties(-z,-zed)
            temp_align.set_up_align_mats()
            temp_align.step_through()
            temp_align.follow_back()
            pos_align_seq.append([temp_align.final_sequence_1,temp_align.final_sequence_2])
            pos_align_scr.append(temp_align.alignment_score)
    
    import math
    all_scr = neg_align_scr + pos_align_scr
    roc_TPRs = []
    roc_FPRs = []

    for i in range(math.floor(min(all_scr))-1,math.floor(max(all_scr))+2):
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

    print(z,zed,end='\r')

    #determine ROC_AUC
    roc_AUC = 0

    for i in range(len(roc_TPRs)-1):
        roc_AUC = roc_AUC + (roc_TPRs[i]+roc_TPRs[i+1])/2*(roc_FPRs[i]-roc_FPRs[i+1])

    score_roc_AUCs[zz] = roc_AUC