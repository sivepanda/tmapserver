from tmapserver.sanitize import sanitize 
import argparse
import anndata as ad

def main():
    parser = argparse.ArgumentParser(description='Sanitize STHD files and spin up a TissUUmaps viewing server locally.')
    parser.add_argument('adata_path', help='Path to save the enrichment analysis results CSV file.')
    args = parser.parse_args()

    adata = ad.read_h5ad(args.adata_path)

    adata.obs.columns = [sanitize(c) for c in adata.obs.columns]

if __name__ == "__main__":
    main()
