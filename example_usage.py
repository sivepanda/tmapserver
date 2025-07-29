#!/usr/bin/env python3
"""
Example usage of tmapserver package for AnnData sanitization and TissUUmaps server integration.
"""

import anndata as ad
from tmapserver import sanitize, run_tissuumaps_server, TissUUmapsServer

def example_programmatic_usage():
    """Example of using tmapserver programmatically."""
    
    # Load your AnnData file
    adata = ad.read_h5ad("your_data.h5ad")
    
    # Method 1: Simple sanitization only
    print("Original columns:", list(adata.obs.columns))
    adata.obs.columns = [sanitize(c) for c in adata.obs.columns]
    print("Sanitized columns:", list(adata.obs.columns))
    
    # Save sanitized data
    adata.write_h5ad("sanitized_data.h5ad")
    
    # Method 2: Sanitize and start server (simple)
    server = run_tissuumaps_server(adata, port=5000, open_browser=True)
    print(f"Server running at: {server.server_url}")
    
    # Method 3: More control with TissUUmapsServer class
    with TissUUmapsServer() as server:
        # Ensure TissUUmaps is installed
        server.ensure_tissuumaps_available()
        
        # Convert data and start server
        project_path, csv_path = server.convert_anndata_for_tissuumaps(adata, "/tmp/tissuumaps_data")
        url = server.start_server(project_path, port=5001)
        
        print(f"Server started at: {url}")
        
        # Do something while server is running
        input("Press Enter to stop server...")
        
    # Server automatically stops when exiting the 'with' block

def example_cli_usage():
    """Example CLI usage (run these commands in terminal)."""
    
    examples = [
        # Basic sanitization only
        "tmapserver input.h5ad output.h5ad",
        
        # Sanitize and start server (auto-opens browser)
        "tmapserver input.h5ad output.h5ad --server",
        
        # Start server on specific port, don't open browser
        "tmapserver input.h5ad output.h5ad --server --port 8080 --no-browser",
        
        # Keep server running (use Ctrl+C to stop)
        "tmapserver input.h5ad output.h5ad --server --keep-server",
        
        # Full example with all options
        "tmapserver input.h5ad output.h5ad --server --port 3000 --no-browser --keep-server"
    ]
    
    print("CLI Usage Examples:")
    for example in examples:
        print(f"  {example}")

if __name__ == "__main__":
    print("tmapserver Usage Examples")
    print("=" * 40)
    
    print("\nCLI Usage:")
    example_cli_usage()
    
    print("\n\nProgrammatic Usage:")
    print("See the example_programmatic_usage() function in this file.")
    print("Note: Uncomment and modify the function to run with your data.")
    
    # Uncomment the line below to run the programmatic example
    # example_programmatic_usage() 