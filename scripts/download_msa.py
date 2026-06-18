#!/usr/bin/env python3
"""
download_msa.py
---------------
Bootstrap script to download and extract the compressed MSA datasets (20,295 files, ~1.3 GB)
from Google Drive to a local target folder (default: msa/).

Usage:
  python3 toga_protocols/scripts/download_msa.py --output msa
"""

import os
import sys
import argparse
import zipfile
import shutil

DEFAULT_FILE_ID = "1aVLhdJCAefvuRxOasuP5QHutTDl9p9NA"

def main():
    parser = argparse.ArgumentParser(description="Bootstrap download of TOGA MSA datasets from Google Drive.")
    parser.add_argument("--file_id", type=str, default=DEFAULT_FILE_ID,
                        help="Google Drive File ID for the MSA zip archive.")
    parser.add_argument("--output", type=str, default="msa",
                        help="Local target directory to extract the alignments to (default: msa).")
    args = parser.parse_args()

    # 1. Automatically check and install gdown if not present
    try:
        import gdown
    except ImportError:
        print("[*] 'gdown' package is required to handle Google Drive downloads. Installing 'gdown'...")
        import subprocess
        try:
            subprocess.run([sys.executable, "-m", "pip", "install", "gdown"], check=True)
            import gdown
        except Exception as e:
            print(f"[!] Error: Failed to automatically install 'gdown': {e}")
            print("[!] Please run 'pip install gdown' manually and run this script again.")
            sys.exit(1)

    # 2. Setup paths
    args.output = os.path.abspath(args.output)
    output_dir = os.path.dirname(args.output)
    os.makedirs(output_dir, exist_ok=True)
    
    zip_path = os.path.join(output_dir, "toga_msa.zip")
    url = f"https://drive.google.com/uc?id={args.file_id}"

    # 3. Download from Google Drive
    print(f"[*] Downloading MSA archive from Google Drive (File ID: {args.file_id})...")
    try:
        gdown.download(url, zip_path, quiet=False)
    except Exception as e:
        print(f"[!] Download failed: {e}")
        sys.exit(1)

    # 4. Extract Zip File and strip common path prefixes robustly
    print(f"[*] Extracting archive to '{args.output}'...")
    if not os.path.exists(zip_path):
        print(f"[!] Error: Downloaded zip file not found at {zip_path}")
        sys.exit(1)

    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            namelist = zip_ref.namelist()
            if not namelist:
                print("[!] Error: The zip archive is empty.")
                sys.exit(1)

            # Detect common path prefixes in the zip (e.g. Users/sergei/Dropbox/TOGA2026/msa/)
            split_paths = [name.split('/') for name in namelist if name.strip()]
            min_len = min(len(p) for p in split_paths)
            common_parts = []
            for i in range(min_len - 1):
                part = split_paths[0][i]
                if all(p[i] == part for p in split_paths):
                    common_parts.append(part)
                else:
                    break
            
            common_prefix = ""
            if common_parts:
                common_prefix = "/".join(common_parts) + "/"
                print(f"[*] Detected directory prefix in zip: '{common_prefix}' (stripping during extraction)")
            else:
                print("[*] No common directory prefix detected in zip. Extracting directly.")

            os.makedirs(args.output, exist_ok=True)

            # Extract members
            count = 0
            for member in zip_ref.infolist():
                if member.is_dir():
                    continue
                
                # Strip prefix
                rel_path = member.filename
                if common_prefix and rel_path.startswith(common_prefix):
                    rel_path = rel_path[len(common_prefix):]
                
                target_file_path = os.path.join(args.output, rel_path)
                os.makedirs(os.path.dirname(target_file_path), exist_ok=True)
                
                with zip_ref.open(member) as source, open(target_file_path, "wb") as target:
                    shutil.copyfileobj(source, target)
                count += 1

            print(f"🎉 Successfully extracted {count} alignment files to '{args.output}'")
            
    except Exception as e:
        print(f"[!] Extraction failed: {e}")
        sys.exit(1)
    finally:
        # Clean up temporary zip file
        if os.path.exists(zip_path):
            print("[*] Cleaning up temporary zip file...")
            os.remove(zip_path)

if __name__ == "__main__":
    main()
