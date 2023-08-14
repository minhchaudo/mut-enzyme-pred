import sys
import os
PATH = os.path.dirname(__file__)
import argparse
import torch
sys.path.append(PATH)

arg_list = sys.argv[1:]
parser = argparse.ArgumentParser()
parser.add_argument("-train", action="store_true", help="If flag is present, switch to training mode (prediction is the default)")
parser.add_argument("-reg", action="store_true", help="If flag is present, switch to regression task (classification is the default)")
parser.add_argument("-subs", action="store_true", help="If flag is present, substrate information (SMILES) is incorporated.")
parser.add_argument("-data", help="Path to the csv file for predictions")
parser.add_argument("-o", help="Path for the output prediction file")
parser.add_argument("-m", help="Path to model, defaults to /models/xgb_reg.pkl for regression and /models/xgb_class.pkl for classification")
parser.add_argument("-s", help="Path to feature selector, defaults to /models/xgb_selector.pkl")
parser.add_argument("-c", help="Path to descriptors columns, defaults to /models/descriptors_columns.txt")
parser.add_argument("-traindata", help="Path to training data")
parser.add_argument("-ms", help="Path to save model, defaults to /models/my_model.pkl")
parser.add_argument("-ss", help="Path to save feature selector, defaults to /models/my_selector.pkl")
parser.add_argument("-cs", help="Path to save descriptors columns, defaults to /models/my_desc_cols.txt")
parser.add_argument("-ts", help="Size of the test dataset (defaults to 0.1)")

args = parser.parse_args()

if __name__ == "__main__":
    if torch.cuda.is_available():
        print("Using GPU")
        device = torch.device('cuda:0')
    else:
        print("No GPU detected. Using CPU, which will take significantly longer.")
        device = torch.device('cpu')
    if args.train:
        try:
            assert args.traindata and args.ms
            ms = args.ms if args.ms is not None else PATH + "/models/my_model.pkl"
            ss = args.ss if args.ss is not None else PATH + "/models/my_selector.pkl"
            cs = args.cs if args.cs is not None else PATH + "/models/my_desc_cols.txt"
            ts = args.ts if args.ts is not None else .1
            from train import train
            train(args.traindata, device, ms, ss, cs, args.reg, args.subs, ts)
        except AssertionError:
            print("Please provide at least path to training data (-traindata) and path to save the trained model (-ms)!")

    else:
        try:
            assert args.data and args.o
            from predict import predict
            m = args.m if args.m is not None else PATH + "/models/xgb_reg.pkl" if args.reg is True else PATH + "/models/xgb_class.pkl"
            s = args.s if args.s is not None else PATH + "/models/xgb_selector.pkl"
            c = None if args.subs is False else args.c if args.c is not None else PATH + "/models/descriptors_columns.txt"
            predict(args.data, args.o, device, m, s, c)
        except AssertionError:
            print("Please provide at least the path to the data (-data) and the path to save the output predictions (-o)!")