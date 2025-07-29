from tmapserver.sanitize import sanitize 
from tmapserver.tissuumaps_server import run_tissuumaps_server, TissUUmapsServer
import argparse
import anndata as ad
import sys
import signal
import atexit

def signal_handler(signum, frame):
    """Handle Ctrl+C gracefully."""
    print("\nShutting down server...")
    sys.exit(0)

def main():
    parser = argparse.ArgumentParser(description='Sanitize STHD files and spin up a TissUUmaps viewing server locally.')
    parser.add_argument('adata_path', help='Path to read adata object.')
    parser.add_argument('output_path', help='Path to save sanitized adata object.')
    parser.add_argument('--server', action='store_true', 
                       help='Start TissUUmaps server after sanitization')
    parser.add_argument('--port', type=int, default=5000,
                       help='Port for TissUUmaps server (default: 5000)')
    parser.add_argument('--no-browser', action='store_true',
                       help='Don\'t automatically open browser')
    parser.add_argument('--keep-server', action='store_true',
                       help='Keep server running (press Ctrl+C to stop)')
    args = parser.parse_args()

    # Register signal handler for graceful shutdown
    signal.signal(signal.SIGINT, signal_handler)
    
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

        # Start TissUUmaps server if requested
        if args.server:
            print("\nStarting TissUUmaps server...")
            try:
                server = run_tissuumaps_server(
                    adata, 
                    port=args.port, 
                    open_browser=not args.no_browser
                )
                
                # Register cleanup function
                atexit.register(lambda: server.stop_server())
                
                if args.keep_server:
                    print("\nServer is running. Press Ctrl+C to stop.")
                    print(f"Access your data at: {server.server_url}")
                    try:
                        # Keep the main thread alive
                        signal.pause()
                    except KeyboardInterrupt:
                        print("\nStopping server...")
                else:
                    print(f"\nServer started at: {server.server_url}")
                    print("Note: Server will stop when this program exits.")
                    print("Use --keep-server flag to keep it running.")
                    
            except Exception as e:
                print(f"Error starting TissUUmaps server: {e}")
                print("The data has been sanitized and saved successfully.")
                return 1
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
