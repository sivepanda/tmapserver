from tmapserver.sanitize import sanitize 
import argparse
import anndata as ad

def main():
    parser = argparse.ArgumentParser(description='Sanitize STHD files and spin up a TissUUmaps viewing server locally.')
    parser.add_argument('adata_path', help='Path to read adata object.')
    parser.add_argument('output_path', help='Path to save sanitized adata object.')
    args = parser.parse_args()

    adata = ad.read_h5ad(args.adata_path)

    adata.obs.columns = [sanitize(c) for c in adata.obs.columns]
    
    adata.write_h5ad(args.output_path)

if __name__ == "__main__":
    main()
