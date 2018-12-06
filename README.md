# Notes
* This project aims to model the reaction time distributions of true cognitive processing and error guessing. 
* This project is collaborated with Prof. Yang Zhang at SooChow University
* The model included here is mainly based on the ref.

	> Glickman, M. E., Gray, J. R., & Morales, C. J. (2005). Combining speed and accuracy to assess error-free cognitive processes. psychometrika, 70(3), 405-425. 


# Files
### model files
* runfitting.m
	* The main program to organize the data and supply to fitting code

* fitRTmodel_optimize
	* The main program that implements model fitting for a single subject 
	* The optimization here utilizes matlab internal optimization process

* fitRTmodel_EM
	* (to be continued..) 	

### data files

* taskDifficultyOnIORExp1.mat
	* mat data file, which contains the data for 54 subjects 	
