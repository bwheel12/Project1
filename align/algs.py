
class PairwiseAligner():
    
    def __init__(self,sequence_1, sequence_2, scoring_matrix_name):
        with open('./Project1-main/'+ sequence_1) as f:
            lines = f.readlines()
        sequence_temp_1 = ""
        for i in range(1,len(lines)):
            sequence_temp_1 = sequence_temp_1 + lines[i][0:(lines[i].find('\n'))]
        
        with open('./Project1-main/'+ sequence_2) as g:
            lines_1 = g.readlines()
        sequence_temp_2 = ""
        for i in range(1,len(lines_1)):
            sequence_temp_2 = sequence_temp_2 + lines_1[i][0:(lines_1[i].find('\n'))]
        
        self.sequence_1 = " " + sequence_temp_1
        self.sequence_2 = " " + sequence_temp_2
        self.scoring_matrix_name = scoring_matrix_name
        
    #define function for reading the score matrices
    def read_scoring_mat(self):
        import numpy as np
        scoring_matrix_temp = np.loadtxt(self.scoring_matrix_name,dtype ='str')
        self.scoring_matrix = np.array([[int(scoring_matrix_temp[j,i]) for i in range(0,24)] for j in range(1,25)])
        self.scoring_lookup = scoring_matrix_temp[0,]
    
    def set_gap_penalties(self,opened,extended):
        self.opened = opened
        self.extended = extended
        
    
    #define function to initialize gap matrices
    #define function to initialize result matrices. Resulting score and path 
    def set_up_align_mats(self):
        import numpy
        self.mat_Ia = numpy.zeros([len(self.sequence_1),len(self.sequence_2)])
        self.mat_Ib = numpy.zeros([len(self.sequence_1),len(self.sequence_2)])
        self.mat_Ms = numpy.zeros([len(self.sequence_1),len(self.sequence_2)])
        self.mat_Md = numpy.zeros([len(self.sequence_1),len(self.sequence_2)])
        
    def check_score_mat(self,AA_1,AA_2):
        import numpy
        AA_1_loc = int(numpy.where(self.scoring_lookup == AA_1.upper())[0])
        AA_2_loc = int(numpy.where(self.scoring_lookup == AA_2.upper())[0])
        match_score = self.scoring_matrix[AA_1_loc,AA_2_loc]
        return match_score
        
class SmithWaterman(PairwiseAligner):
    #local
    #define function for stepping through and assigning scores and directions
    ##currently global alignement, needs updating
    def step_through(self):
        import numpy as np
        m = len(self.sequence_1)
        n = len(self.sequence_2)
        for i in range(m):
            self.mat_Ms[i,0] = 0 
            self.mat_Md[i,0] = 1
            self.mat_Ia[i,0] = 0
            self.mat_Ib[i,0] = 0
        for j in range(n):
            self.mat_Ms[0,j] = 0
            self.mat_Md[0,j] = 2
            self.mat_Ia[0,j] = 0
            self.mat_Ib[0,j] = 0
        for k in range(1,m):
            for l in range(1,n):
                #idea here is to d all calculations for decisions in 2x2 matrices for Ia and Ib and 3 x 2 for Ms. Max value can then be decided  with ia and ib, and passed to decision ms, upon determining the ideal way to proceed pass appropriate score and decision to final Ms and Md matrices
                if (k == 0 ) and (l == 0):
                    self.mat_Ms[k,l] = 0
                    self.mat_Md[k,l] = 0
                    self.mat_Ia[k,l] = 0
                    self.mat_Ib[k,l] = 0
                 
                if (k > 0) or (l > 0):
                    decision_Ia = np.array([self.mat_Ms[k-1,l]+self.opened,self.mat_Ia[k-1,l]+self.extended])
                    self.mat_Ia[k,l]     = max(decision_Ia)
                    decision_Ib = np.array([self.mat_Ms[k,l-1]+self.opened,self.mat_Ib[k,l-1]+self.extended])
                    self.mat_Ib[k,l]     = max(decision_Ib)
                    decision_Ms = np.array([self.mat_Ms[k-1,l-1]+self.check_score_mat(self.sequence_1[k],self.sequence_2[l]),self.mat_Ia[k,l],self.mat_Ib[k,l],0])
                    #print(decision_Ms)
                    self.mat_Ms[k,l]     = max(decision_Ms)
                    #direction matrix code 0 = match, 1 = up, 2 = left, 3 = terminated
                    self.mat_Md[k,l]     =int(min(np.where(decision_Ms == self.mat_Ms[k,l]))[0])
                    
                    
                    
                    
    #define function for tracinng back and giving the correctly aligned sequence and score
    def follow_back(self):
        import numpy
        i_possible = numpy.where(self.mat_Ms == numpy.amax(self.mat_Ms))[0]
        j_possible = numpy.where(self.mat_Ms == numpy.amax(self.mat_Ms))[1]
        i = int(i_possible[len(i_possible)-1])
        j = int(j_possible[len(j_possible)-1])
        self.alignment_score = self.mat_Ms[i,j]
        end = False
        return_sequence_1 = ""
        return_sequence_2 = ""
        while end == False:
            if self.mat_Md[i,j] == 0:
                i_temp = i-1
                j_temp = j-1
                return_sequence_1 = self.sequence_1[i] + return_sequence_1
                return_sequence_2 = self.sequence_2[j] + return_sequence_2
            if self.mat_Md[i,j] == 1:
                i_temp = i-1
                j_temp = j
                return_sequence_1 = self.sequence_1[i] + return_sequence_1
                return_sequence_2 = "-" + return_sequence_2
            if self.mat_Md[i,j] == 2:
                i_temp = i
                j_temp = j-1
                return_sequence_1 = "-" + return_sequence_1
                return_sequence_2 = self.sequence_2[j] + return_sequence_2    
            i = i_temp
            j = j_temp
            if ((i == 0) & (j == 0)) or self.mat_Ms[i,j] == 0:
                end = True
                self.final_sequence_1 = return_sequence_1
                self.final_sequence_2 = return_sequence_2
    
    
class NeedlemanWunsch(PairwiseAligner):
    #global
    
    #define function for stepping through and assigning scores and directions
    def step_through(self):
        import numpy as np
        m = len(self.sequence_1)
        n = len(self.sequence_2)
        for i in range(m):
            self.mat_Ms[i,0] = 0 + i*self.extended
            self.mat_Md[i,0] = 1
            self.mat_Ia[i,0] = -100
            self.mat_Ib[i,0] = -100
        for j in range(n):
            self.mat_Ms[0,j] = 0 + j*self.extended
            self.mat_Md[0,j] = 2
            self.mat_Ia[0,j] = -100
            self.mat_Ib[0,j] = -100
        for k in range(1,m):
            for l in range(1,n):
                #idea here is to d all calculations for decisions in 2x2 matrices for Ia and Ib and 3 x 2 for Ms. Max value can then be decided  with ia and ib, and passed to decision ms, upon determining the ideal way to proceed pass appropriate score and decision to final Ms and Md matrices
                if (k == 0 ) and (l == 0):
                    self.mat_Ms[k,l] = 0
                    self.mat_Md[k,l] = 0
                    self.mat_Ia[k,l] = 0
                    self.mat_Ib[k,l] = 0
                 
                if (k > 0) or (l > 0):
                    decision_Ia = np.array([self.mat_Ms[k-1,l]+self.opened,self.mat_Ia[k-1,l]+self.extended])
                    self.mat_Ia[k,l]     = max(decision_Ia)
                    decision_Ib = np.array([self.mat_Ms[k,l-1]+self.opened,self.mat_Ib[k,l-1]+self.extended])
                    self.mat_Ib[k,l]     = max(decision_Ib)
                    decision_Ms = np.array([self.mat_Ms[k-1,l-1]+self.check_score_mat(self.sequence_1[k],self.sequence_2[l]),self.mat_Ia[k,l],self.mat_Ib[k,l]])
                    #print(decision_Ms)
                    self.mat_Ms[k,l]     = max(decision_Ms)
                    #direction matrix code 0 = match, 1 = up, 2 = left
                    self.mat_Md[k,l]     =int(min(np.where(decision_Ms == self.mat_Ms[k,l]))[0])

        #define function for tracinng back and giving the correctly aligned sequence and score
        #idea put follow_Back into parent. pass start i,j. each child has unique function to update end
    def follow_back(self):
        import numpy
        i = len(self.sequence_1) - 1
        j = len(self.sequence_2) - 1
        end = False
        return_sequence_1 = ""
        return_sequence_2 = ""
        self.alignment_score = self.mat_Ms[i,j]
        while end == False:
            if self.mat_Md[i,j] == 0:
                i_temp = i-1
                j_temp = j-1
                return_sequence_1 = self.sequence_1[i] + return_sequence_1
                return_sequence_2 = self.sequence_2[j] + return_sequence_2
            if self.mat_Md[i,j] == 1:
                i_temp = i-1
                j_temp = j
                return_sequence_1 = self.sequence_1[i] + return_sequence_1
                return_sequence_2 = "-" + return_sequence_2
            if self.mat_Md[i,j] == 2:
                i_temp = i
                j_temp = j-1
                return_sequence_1 = "-" + return_sequence_1
                return_sequence_2 = self.sequence_2[j] + return_sequence_2    
            i = i_temp
            j = j_temp
            if (i == 0) & (j == 0):
                end = True
                self.final_sequence_1 = return_sequence_1
                self.final_sequence_2 = return_sequence_2
                    