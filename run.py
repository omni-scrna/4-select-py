#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

# Add src directory to Python path
src_dir = Path(__file__).parent / "src"
sys.path.insert(0, str(src_dir))

from main import process_data

def parse_args():
    parser = argparse.ArgumentParser(description='OmniBenchmark module')

    # Required by OmniBenchmark
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory for results')
    parser.add_argument('--name', type=str, required=True,
                       help='Module name/identifier')
    # Stage-specific inputs
    parser.add_argument('--normalized.h5', nargs='+', dest='normalized_h5', required=True,
                       help='Input: normalized.h5')
    parser.add_argument('--selection_type', type=str, required=True,
                       choices = ["scanpy_hvg"], #, "pearson_residuals", "giniclust3"],
                       help='Selection method')
    parser.add_argument('--number_selected', type=int, required=True,
                       help='Input: number_selected')
    return parser.parse_args()

def main():
    args = parse_args()

    print(f"Output directory: {args.output_dir}")
    print(f"Module name: {args.name}")
    print(f"normalized.h5: {args.normalized_h5}")
    print(f"selection_type: {args.selection_type}")
    print(f"number_selected: {args.number_selected}")

    # TODO: Implement your module logic
    # Process the data using main function
    process_data(args)

if __name__ == "__main__":
    main()
