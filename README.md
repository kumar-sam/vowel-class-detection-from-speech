----------------------------------------------------------------------------------------
    CONSOLE APPLICATION : 214101048_vowelRecognition ( Assignment 3 : Vowel Recognition)
----------------------------------------------------------------------------------------


Input files :-
-----------
All the data for training are testing are stored in the resources/recordings folder. And all the intermediate files are stored inside processed folder.
Following is the brief of each folder.
(i)   resoures/recordings - contains all 100 recording files
(ii)  resources/processed/normalization - Intermidate generated file: contains the normalized data for all 100 files.
(iii) resources/processed/cc_calc - Intermidate generated file: contains the capstral coefficients of all 100 files. the length of each files is 5*12 = 600 entry.
(iv)  resources/processed/training - Intermidate generated file: contains the capstral coefficients of 5 represntative vowel files. dimension 5*12 entry.   

----------------------------------------------------------------------------------------

How to run the Program :-
----------------------
1) Default Input file is already available and program works well without changing anything just compile and run the program.

2) To Change the input file apply changes as belows:

   o  upload the 100 new files inside the resoures/recordings/ folder. name (214101048_a_1.txt) should follow the following convention.
      [ Name should be 2141010xx_v_1.txt to 214101010xx_v_20.txt and V = a, e, i, o, u  and xx is the roll number.]

   o Inside the program change #define statement as follows :
      
      #define PREFIX_FILE "214101048_" to   #define PREFIX_FILE "2141010XX_"

3) Now simply run the program on the Visual Code.
---------------------------------------------------------------------------------------------------
