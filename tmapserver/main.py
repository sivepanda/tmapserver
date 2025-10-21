from tmapserver.sanitize import sanitize
import argparse
import anndata as ad
import sys

def main():
    parser = argparse.ArgumentParser(description='Sanitize AnnData column names by removing/replacing special characters.')
    parser.add_argument('adata_path', help='Path to read adata object.')
    parser.add_argument('output_path', help='Path to save sanitized adata object.')
    args = parser.parse_args()

    try:
        print(f"Reading AnnData from: {args.adata_path}")
        adata = ad.read_h5ad(args.adata_path)

        print("Sanitizing column names...")
        original_columns = list(adata.obs.columns)
        adata.obs.columns = [sanitize(c) for c in adata.obs.columns]

        # Show what was changed
        changes = []
        for orig, new in zip(original_columns, adata.obs.columns):
            if orig != new:
                changes.append(f"  '{orig}' -> '{new}'")

        if changes:
            print("Column name changes:")
            for change in changes:
                print(change)
        else:
            print("No column names needed sanitization.")

        print(f"Writing sanitized AnnData to: {args.output_path}")
        adata.write_h5ad(args.output_path)
        print("Sanitization complete!")

        return 0

    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 