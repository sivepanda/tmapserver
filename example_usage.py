#!/usr/bin/env python3
"""
Example usage of adata-sanitize package for AnnData column name sanitization.
"""

import anndata as ad
from tmapserver import sanitize

def example_programmatic_usage():
    """Example of using adata-sanitize programmatically."""

    # Load your AnnData file
    adata = ad.read_h5ad("your_data.h5ad")

    # Sanitize column names
    print("Original columns:", list(adata.obs.columns))
    adata.obs.columns = [sanitize(c) for c in adata.obs.columns]
    print("Sanitized columns:", list(adata.obs.columns))

    # Save sanitized data
    adata.write_h5ad("sanitized_data.h5ad")

    # You can also sanitize individual strings
    test_string = "Cell Type (CD45+/CD3-)"
    sanitized = sanitize(test_string)
    print(f"'{test_string}' -> '{sanitized}'")

def example_cli_usage():
    """Example CLI usage (run these commands in terminal)."""

    examples = [
        # Basic sanitization
        "adata-sanitize input.h5ad output.h5ad",
    ]

    print("CLI Usage Examples:")
    for example in examples:
        print(f"  {example}")

if __name__ == "__main__":
    print("adata-sanitize Usage Examples")
    print("=" * 40)

    print("\nCLI Usage:")
    example_cli_usage()

    print("\n\nProgrammatic Usage:")
    print("See the example_programmatic_usage() function in this file.")
    print("Note: Uncomment and modify the function to run with your data.")

    # Uncomment the line below to run the programmatic example
    # example_programmatic_usage() 