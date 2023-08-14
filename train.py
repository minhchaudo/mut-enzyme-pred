from transformers import T5EncoderModel, T5Tokenizer
import torch
import pickle
import pandas as pd
from tqdm import tqdm
from mordred import Calculator, descriptors
from rdkit import Chem
from utils import *
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import train_test_split, RepeatedKFold, RandomizedSearchCV, cross_val_score
from xgboost import XGBRegressor, XGBClassifier
import numpy as np


def get_model(device):
        print("\nDownloading model...")
        model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc").to(device)
        model = model.eval()
        tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', do_lower_case=False)
        return model, tokenizer

def get_embeddings( model, tokenizer, seqs, device,
                   max_residues=4000, max_seq_len=1000, max_batch=100):
    print("\nGenerating embeddings...\n")
    protein_embs = dict()
    seq_dict   = sorted( seqs.items(), key=lambda kv: len( seqs[kv[0]] ), reverse=True )
    batch = list()
    for seq_idx, (pdb_id, seq) in enumerate(tqdm(seq_dict),1):
        seq = seq
        seq_len = len(seq)
        seq = ' '.join(list(seq))
        batch.append((pdb_id,seq,seq_len))
        n_res_batch = sum([ s_len for  _, _, s_len in batch ]) + seq_len
        if len(batch) >= max_batch or n_res_batch>=max_residues or seq_idx==len(seq_dict) or seq_len>max_seq_len:
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = list()
            token_encoding = tokenizer.batch_encode_plus(seqs, add_special_tokens=True, padding="longest")
            input_ids      = torch.tensor(token_encoding['input_ids']).to(device)
            attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)
            try:
                with torch.no_grad():
                    embedding_repr = model(input_ids, attention_mask=attention_mask)
            except RuntimeError:
                print("RuntimeError during embedding for {} (L={})".format(pdb_id, seq_len))
                continue
            for batch_idx, identifier in enumerate(pdb_ids):
                s_len = seq_lens[batch_idx]
                emb = embedding_repr.last_hidden_state[batch_idx,:s_len]
                protein_emb = emb.mean(dim=0)
                protein_embs[identifier] = protein_emb.detach().cpu().numpy().squeeze()
    return pd.DataFrame.from_dict(protein_embs, 'index')

def preprocess_descriptors(df, drop_pct):
    df = df.copy()
    to_drop = []
    for col in df.columns:
        try:
            df[col] = df[col].apply(lambda x: float(x))
        except:
            to_drop.append(col)
    df.drop(columns=to_drop, inplace=True)
    to_drop = []
    for col in df.columns:
        if df[col].isna().sum() > len(df[col])*drop_pct:
            to_drop.append(col)
    df.drop(columns=to_drop, inplace=True)
    to_drop = []
    for col in df.columns:
        if min(df[col]) == max(df[col]):
            to_drop.append(col)
    df.drop(columns=to_drop, inplace=True)
    return df    


def get_descriptors(smiles, id, desc_columns_save_path):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    calc = Calculator(descriptors)
    df = calc.pandas(mols)
    df.index = id
    df = preprocess_descriptors(df, .1)
    pd.Series(df.columns).to_csv(desc_columns_save_path, index=False, header=False)
    return df

def feature_select(X_train, X_test, Y_train, algo, selector_save_path):
    model = algo()
    model.fit(X_train, Y_train)
    selector = SelectFromModel(model, prefit=True)
    X_train_new = selector.transform(X_train)
    X_test_new = selector.transform(X_test)
    pickle.dump(selector, open(selector_save_path, "wb"))
    return X_train_new, X_test_new

def optimize(X_train, Y_train, algo, n_fold=5, n_repeats=3):
    model = algo()
    max_depth = range(4,8)
    eta = [0.01, 0.1, 0.2, 0.3]
    min_child_weight = range(1,6,2)
    grid = dict(max_depth=max_depth, eta=eta, min_child_weight=min_child_weight)
    cv = RepeatedKFold(n_splits=n_fold, n_repeats=n_repeats)
    randomSearch = RandomizedSearchCV(estimator=model,
        cv=cv, param_distributions=grid,
        scoring="r2")
    searchResults = randomSearch.fit(X_train, Y_train)
    opt_model = algo(**searchResults.best_params_)
    opt_model.fit(X_train, Y_train)
    return opt_model, searchResults.best_params_

def cross_val_and_test(model, X_train, Y_train, X_test, Y_test, scoring='r2', n_fold=5, n_repeats=3):
    cv = RepeatedKFold(n_splits=n_fold, n_repeats=n_repeats)
    scores = cross_val_score(model, X_train, Y_train, scoring=scoring, cv=cv)
    scores = abs(scores)
    print(f"{scoring} of optimized model over {n_fold}-fold cross validation with {n_repeats} repetitions: mean {scores.mean()}, std {scores.std()}")
    print(f"Score of optimized model on test set: {model.score(X_test, Y_test)}")

def train_full_model(X, Y, algo, params, model_save_path):
    model = algo(**params)
    algo.fit(X, Y)
    pickle.dump(model, open(model_save_path, "wb"))

def train_and_save_model(X, Y, algo, test_size, model_save_path, selector_save_path):
    X = X.to_numpy()
    Y = Y.to_numpy()
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size)
    X_train_new, X_test_new = feature_select(X_train, X_test, Y_train, algo, selector_save_path)
    opt_model, opt_params = optimize(X_train_new, Y_train, algo)
    cross_val_and_test(opt_model, X_train_new, Y_train, X_test_new, Y_test)
    train_full_model(np.concatenate([X_train_new, X_test_new]), np.concatenate([Y_train, Y_test]), algo, opt_params, model_save_path)
    print(f"Full model saved to {model_save_path}")


def train(datapath, device, model_save_path, selector_save_path, desc_columns_save_path, reg, subs, test_size):
    try:
        df = pd.read_csv(datapath).set_index("id")
        wt_seqs = df["wt_seqs"].to_dict()
        mut_seqs = df["mut_seqs"].to_dict()
        model, tokenizer = get_model(device)
        wt_embs_df = get_embeddings(model, tokenizer, wt_seqs, device)
        mut_embs_df = get_embeddings(model, tokenizer, mut_seqs, device)
        if subs is True:
            if "smiles" in df.columns:
                smiles, id = df["smiles"], df.index
                descriptors = get_descriptors(smiles, id, desc_columns_save_path)
                full_df = pd.concat([wt_embs_df, mut_embs_df, descriptors, df["val"]], axis=1)
            else:
                print("Warning: Substrate information is expected to be incorporated, but your data doesn't contain the 'smiles' column. Ignoring the -subs flag.")
        else:
            full_df = pd.concat([wt_embs_df, mut_embs_df, df["val"]], axis=1)
        if reg:
            algo = XGBRegressor
        else:
            algo = XGBClassifier
        train_and_save_model(full_df.iloc[:,:-1], full_df["val"], algo, test_size, model_save_path, selector_save_path)
    except Exception as error:
        print("An error occurred!!!")
        print(error)
