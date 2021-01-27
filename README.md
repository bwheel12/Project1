# Project 1 - Sequence Alignment Benjamin Wheeler
## Due 01/27/2021

![BuildStatus](https://github.com/bwheel12/Project1/workflows/HW1/badge.svg?event=workflow_dispatch)

In this assignment, you will implement two classical alignment algorithms and then evaluate each algorithm’s performance with a range of parameters. There are two parts to this assignment and Part 2 requires completion of Part 1. We recommend reading through both Part 1 and Part 2 before beginning this assignment. 

* Part 1 - API and implementation
* Part 2 - Evaluating alignments

### main
Runs all code in align/\_\_main\_\_.py, useful for part 2
```
python -m align
```

### testing
Testing is as simple as running
```
python -m pytest
```
from the root directory of this project.

## To align two given sequences the workflow is as follows
```
from align import algs
temp_align = algs.SmithWaterman('sequence_1_location.txt','sequence_2_location.txt','scoring_matrix_location.mat')
temp_align.read_scoring_mat()
temp_align.set_gap_penalties(-11,-3)
temp_align.set_up_align_mats()
temp_align.step_through()
temp_align.follow_back()
```


## Main Contents

The main file contains all the lines of code to run the exercises in part two of the assignment. When single values were requested the are generally printed to the terminal. Matrices and graphs generated are then saved to the main directory of the project. 

## Functions Inherited from Parent PairwiseAligner Class
The class PairwiseAligner is not intended to be called. Rather its two childclasses are to be interacted with directly. The following functions are inherited in each class.

### Initialization of an alignment class
To align two sequences, either a smith-waterman or NeedlemanWunsch object must be initialized. Upon initialization three variables must be passed: a string pointing to the location of a fasta containing sequence 1, a string pointing to the location of a fasta containing sequence 2, and a string pointing to the location of the scoring matrix to be used.

### read_scoring_mat
This function does not take any arguments and simply instructs the software to import the actual values from the scoring matrix file. No argument is given.

### set_gap_penalties
This function is used to set the gap opening and extension penalties to be used. These values are not provided at initilization, because the intent is that the same pair of sequences would possibly be run with multiple gap penalties. This allows gap penalties to be reset quickly and the alignement rerun with out the need to load the sequences again. This function takes two integer arguments. The first is presumed to be the gap opening penalty and the second is the gap extension penalty. These values should be negative.

### set_up_alignment_mats
This function creates the initial empty matrices for the alignment process. Should be called by the user before and inbetween alignment runs. No argument is given.

### check_score_mat
This function looks up the score for a given amino acid pair from the provided score matrix and returns that score. The amino acid pair to be score should each be provided as individual character arguments. While possibe, this is not intended for the user to call. This will be used internally by the step_through function

## SmithWaterman Class and its private functions
This class is a child of the parent PairwiseAligner class. It inherits all the functions above and implements two functions of its own. These fucntions deal with the actual alignment implementation and save the results to atributes within the SmithWaterman object. These alignments are designed for local alignment using the SmithWaterman algorithm.

### step_through
This function interates through the alignment matrices determining the appropriate score and path according to the SmithWaterman algorithm. The function should be called by the user with no arguments passed. The function will update the alignment matrices as atributes of the alignment object. 

### follow_back
This function finds the optimal alignment within the score matrix. It then follows it back to yield the appropriate, locally aligned sequences these are assigned to sequence_1 and sequence_2 atrributes of the SmithWaterman Object. It further assigns the score to the alignment_score attribute of the SmithWaterman object. The user should call this function with no arguments.

## NeedlemanWunsch Class and its private functions
This class is a child of the parent PairwiseAligner class. It inherits all the functions above and implements its own versions of step_through and follow_back. These functions deal with the actual alignment implementation and save the results to attributes within the NeedlemanWunsch object. These alignments are desigend for global alignment using the NeedlemanWunsch algorithm.

### step_through
This function interates through the alignment matrices determining the appropriate score and path according to the NeedlemanWunsch algorithm. The function should be called by the user with no arguments passed. The function will update the alignment matrices as atributes of the alignment object. 

### follow_back
This function finds the optimal alignment within the score matrix. It then follows it back to yield the appropriate, globally aligned sequences these are assigned to sequence_1 and sequence_2 atrributes of the NeedlemanWunsch Object. It further assigns the score to the alignment_score attribute of the NeedlemanWunsch object. The user should call this function with no arguments.

## Dependencies
These sets of code depend on the following packages: numpy, matplotlib, pandas, math