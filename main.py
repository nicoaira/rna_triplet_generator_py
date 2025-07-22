#!/usr/bin/env python3
"""
RNA Triplet Generator - Python Version

This is a modular Python implementation of the RNA triplet generator,
originally written in Rust. It generates anchor/positive/negative RNA
triplets with structural modifications for machine learning datasets.
"""

import logging
import sys
from pathlib import Path
from typing import Optional

from config.cli_parser import parse_arguments
from core.dataset_generator import DatasetGenerator
from utils.logger import setup_logging
from utils.output_handler import OutputHandler


def main() -> int:
    """Main entry point for the RNA triplet generator."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Setup logging
        setup_logging(debug=args.debug)
        logger = logging.getLogger(__name__)
        
        logger.info("Starting RNA Triplet Generator (Python version)")
        logger.info(f"Generating {args.num_structures} triplets")
        logger.info(f"Output directory: {args.output_dir}")
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(exist_ok=True)
        
        # Initialize dataset generator
        generator = DatasetGenerator(args)
        
        # Generate dataset
        logger.info("Generating triplet dataset...")
        triplets = generator.generate_dataset()
        
        # Handle output
        output_handler = OutputHandler(args)
        output_handler.save_dataset(triplets, output_dir)
        
        logger.info("Dataset generation completed successfully!")
        return 0
        
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())