from transformers import T5EncoderModel, T5Tokenizer
import torch
import pickle
import pandas as pd
from tqdm import tqdm
from mordred import Calculator, descriptors
from rdkit import Chem
from utils import *

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

def get_descriptors(smiles, id, desc_cols_path):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    calc = Calculator(descriptors)
    df = calc.pandas(mols)
    descriptors_cols = list(pd.read_csv(desc_cols_path, header=None)[0])
    df = df[descriptors_cols]
    df.index = id
    return df

def get_preds_df(model_path, selector_path, wt_embs, mut_embs, descriptors=None):
    if descriptors is not None:
        X = pd.concat([wt_embs, mut_embs, descriptors], axis=1)
    else:
        X = pd.concat([wt_embs, mut_embs], axis=1)
    index = pd.Series(X.index)
    X = X.to_numpy()
    selector = pickle.load(open(selector_path, "rb"))
    X_new = selector.transform(X)
    model = pickle.load(open(model_path, "rb"))
    preds = pd.Series(model.predict(X_new))
    preds_df = pd.concat([index, preds], axis=1)
    preds_df.columns = ["ID", "predictions"]
    return preds_df

def predict(data_path, out_path, device, model_path, selector_path, desc_cols_path=None):
    try:
        df = pd.read_csv(data_path).set_index("id")
        wt_seqs = df["wt_seqs"].to_dict()
        mut_seqs = df["mut_seqs"].to_dict()
        model, tokenizer = get_model(device)
        wt_embs_df = get_embeddings(model, tokenizer, wt_seqs, device)
        mut_embs_df = get_embeddings(model, tokenizer, mut_seqs, device)
        if desc_cols_path is not None:
            if "smiles" in df.columns:
                smiles, id = df["smiles"], df.index
                descriptors = get_descriptors(smiles, id, desc_cols_path)
                preds_df = get_preds_df(model_path, selector_path, wt_embs_df, mut_embs_df, descriptors)
            else:
                print("Warning: Substrate information is expected, but the data doesn't contain the 'smiles' column. Ignoring the '-subs' flag.")
        else:
            preds_df = get_preds_df(model_path, selector_path, wt_embs_df, mut_embs_df)
        preds_df.to_csv(out_path, index=False)
        print(f"Predictions saved to {out_path}")
    except Exception as error:
        print("Oops, an error occurred!")
        print(error)

