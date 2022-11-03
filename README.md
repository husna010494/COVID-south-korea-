# COVID-south-korea-

1. The dataset
patient.csv: https://www.kaggle.com/datasets/kimjihoo/coronavirusdataset?select=PatientInfo.csv. (patientinfo.csv)

 a. patient.id: patient id from south korea earlier 2020.
 
 b. sex: 1: male, 2: female
 
 c. age: the random number from original dataset (the original dataset is categorical data (10s, 20, ...), so we make random number where 10s is random number from 10-19)
 
 d. province: the province of patient
 
 e. state: released (uncensored dataset), isolated (censored dataset)
 
 f. time: the difference date from confirm date and released date (the maximum date is 5th August 2022)
 
 g. status: 1: uncensored , 0: censored
 
2. The library: MASS, survival, randomForestSRC, alabama, tidyverse, survivalROC.
3. The code for stacked method is in COVID.R. Please, change the root for specific folder.
4. To know the weight, just write stacked.est$alphas
5. Result 

a. The qqplot shows that the lognormal distribution is better distribution for the dataset

b. The ROC of stacked method and rsf has similar result
 
#Installation instructions / Run Instructions

Step 1: Install all libraries in your R.

Step 2a: Open "COVID.R"

Step 2b: Run all code in "Covid model.R" after change the spesific root.

Step 3: Run all code "ROC.R" and then change the marker by type model.

Step 4: write summary(mod.ln) and summary(mod.cox) for result of each model.

Step 5: write stacked.est$alphas to weight from each model for stacked model.

#Requirements

This was created and tested on an R X64. 4.0.3.
