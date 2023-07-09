# Dataset Outline 
1. **data_dump_density_preprocessed_V2.pk** -----> All data containing 4050 data samples

2. **data_dump_density_preprocessed_train.pk** ----->  up to 3550 samples selected randomly (via random splitting method) for training.   
    **data_dump_density_preprocessed_test.pk** -----> the remaining 500 samples are used for testing      

3. **data_dump_density_preprocessed_train_excludeX.pk** -----> Exclude N(X) samples involving [X] as input parameters (via deterministic splitting method), thus generating a training sample set of size = 4050 - N(X)    
    **data_dump_density_preprocessed_test_excludeX.pk** -----> Corresponding excluded testing dataset has N(X) samples   

4. **scaler_new.pkl** -----> Scaler modules fitting with the maximum and the minimum of the original input parameter spaces  
