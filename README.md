# Mutated enzyme efficiency prediction (MEEP)

### 1. Introduction:

This is a package to predict the change in enzyme efficiency upon point mutation. Change in enzyme efficiency was caculated as $log_{10}\frac{kcat_{mut}/Km_{mut}}{kcat_{wt}/Km_{wt}}$

Our model was trained on over 2500 data points retrieved from [IntEnzyDB](https://intenzydb.accre.vanderbilt.edu) on 2 tasks, regression and binary classification (whether the mutation is enhancing, which is positive class 1, or debilitating, which is negative class 0). Our regression model achieved R2 = 0.49 and MSE = 0.96. Our classification model achieved accuracy 0.87 and MCC 0.49.

### 2. Installation:

To install the project, run the following code in the command line:

```
git clone https://github.com/minhchaudo/meep.git
cd meep
pip install -r requirements.txt
cd ..
```

### 3. Usage:

MEEP can be used for either prediction or training.
For prediction using our trained model, provide at least the path to your data and the path to write the output predictions.

```
python -m meep -data your/data/path.csv -o your/prediction/path.csv
```

The default task is classification. You can include the -reg flag to switch to regression.

```
python -m meep -reg -data your/data/path.csv -o your/prediction/path.csv
```

For further customization, you can specify the path to your own model, the path to your feature selector, and the path to the file that details the substrate descriptors you want to include.

```
python -m meep -data your/data/path.csv -o your/prediction/path.csv -m your/path/to/model.pkl -s your/path/to/selector.pkl -c your/path/to/descriptors.txt
```

For training a classification model, please specify at least the path to the training data and the path to save the output model.

```
python -m meep -traindata your/data/path.csv -ms your/model/save/path.pkl
```

For training a regression model, please include the -reg flag.

```
python -m meep -reg -traindata your/data/path.csv -ms your/model/save/path.pkl
```

For model validation and testing, can also specify the proportion of the test set.

```
python -m meep -reg -traindata your/data/path.csv -ms your/model/save/path.pkl -ts 0.1
```

For a list of complete arguments, please run

```
python -m meep -h
```
