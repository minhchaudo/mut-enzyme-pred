# Mutated enzyme efficiency prediction (MEEP)

### 1. Introduction:

This is a package to predict the change in enzyme activity upon point mutation. Our classification model was trained on over 25000 mutations retrieved from [D3DistalMutation](https://www.d3pharma.com/D3DistalMutation/index.php), with are labeled as either enhancing (class 1) or debilitating (class 0). Our regression model was trained on over 2600 data points retrieved from [IntEnzyDB](https://intenzydb.accre.vanderbilt.edu/) to predict the change in enzyme efficiency, which is calculated as $log_{10}\frac{kcat_{mut}/Km_{mut}}{kcat_{wt}/Km_{wt}}$.

### 2. Installation:

To install the project, run the following code in the command line:

```
git clone https://github.com/minhchaudo/meep.git
cd meep
pip install -r requirements.txt
cd ..
```

### 3. Usage:

MEEP can be used for either prediction or training with your own dataset.
For prediction using our trained model, provide at least the path to your data and the path to write the output predictions.

```
python -m meep -data your/data/path.csv -o your/prediction/path.csv
```

Your data file is expected to be a .csv file with separator "," and contains the following columns: "id", "wt_seqs", "mut_seqs".

The default task is classification. You can include the -reg flag to switch to regression.

```
python -m meep -reg -subs -data your/data/path.csv -o your/prediction/path.csv
```

For further customization, you can specify the path to your own model, the path to your feature selector, and the path to the file that details the substrate descriptors you want to include.

```
python -m meep -data your/data/path.csv -o your/prediction/path.csv -m your/path/to/model.pkl -s your/path/to/selector.pkl -c your/path/to/descriptors.txt
```

Our default classification model doesn't require substrate smiles for prediction, but our regression model does. Therefore, if you want to use our models for prediction, please include the -subs tag for the regression model and exclude it for the classification model. If you are using your own model for prediction, add the flag when appropriate. When the -subs flag is added, you must also include the "smiles" column in your data.

For training a classification model, please specify at least the path to the training data and the path to save the output model.

```
python -m meep -train -traindata your/data/path.csv -ms your/model/save/path.pkl
```

For training a regression model, please include the -reg flag.

```
python -m meep -train -reg -traindata your/data/path.csv -ms your/model/save/path.pkl
```

Your training data file is expected to be a .csv file with separator "," and contains the following columns: "id", "wt_seqs", "mut_seqs", and "val". The "val" column should contain float values in case of regression and binary values (0 and 1) in the case of classification. To illustrate, please take a look at the files data_li/df_li_class.csv and data_li/df_li_reg.csv.

If you want to incorporate substrate information into your model, please add the -subs tag. Your training data is then expected to contain, in addition to the above, the "smiles" column.

For model validation and testing, you can also specify the proportion of the test set.

```
python -m meep -train -reg -traindata your/data/path.csv -ms your/model/save/path.pkl -ts 0.1
```

For a list of complete arguments, please run

```
python -m meep -h
```
