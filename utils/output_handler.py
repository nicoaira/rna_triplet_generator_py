"""
Output handling utilities for saving datasets and metadata.
"""

import csv
import json
import logging
import random
from pathlib import Path
from typing import List, Tuple

logger = logging.getLogger(__name__)


class OutputHandler:
    """Handles saving datasets and metadata to various formats."""
    
    def __init__(self, args):
        """Initialize the output handler with configuration."""
        self.args = args
    
    def save_dataset(self, triplets: List, output_dir: Path) -> None:
        """
        Save the complete dataset with metadata.
        
        Args:
            triplets: List of RNA triplets to save
            output_dir: Directory to save files
        """
        logger.info(f"Saving dataset with {len(triplets)} triplets to {output_dir}")
        
        # Create metadata
        from core.models import DatasetMetadata
        metadata = DatasetMetadata.create(self.args, len(triplets))
        
        # Save main dataset
        if self.args.split:
            self._save_split_dataset(triplets, metadata, output_dir)
        else:
            self._save_single_dataset(triplets, metadata, output_dir)
        
        # Save metadata
        self._save_metadata(metadata, output_dir / "metadata.json")
        
        logger.info("Dataset saved successfully!")
    
    def _save_single_dataset(self, triplets: List, 
                           metadata, output_dir: Path) -> None:
        """Save dataset as a single file."""
        csv_path = output_dir / "rna_triplets.csv"
        self._save_triplets_csv(triplets, csv_path, metadata)
    
    def _save_split_dataset(self, triplets: List, 
                          metadata, output_dir: Path) -> None:
        """Save dataset split into train and validation sets."""
        # Shuffle triplets
        shuffled_triplets = triplets.copy()
        random.shuffle(shuffled_triplets)
        
        # Calculate split sizes
        total_size = len(shuffled_triplets)
        train_size = int(total_size * self.args.train_fraction)
        
        # Split the data
        train_triplets = shuffled_triplets[:train_size]
        val_triplets = shuffled_triplets[train_size:]
        
        logger.info(f"Split dataset: {len(train_triplets)} train, {len(val_triplets)} validation")
        
        # Save splits
        train_path = output_dir / "rna_triplets_train.csv"
        val_path = output_dir / "rna_triplets_val.csv"
        full_path = output_dir / "rna_triplets.csv"
        
        self._save_triplets_csv(train_triplets, train_path, metadata)
        self._save_triplets_csv(val_triplets, val_path, metadata)
        self._save_triplets_csv(triplets, full_path, metadata)  # Also save full dataset
    
    def _save_triplets_csv(self, triplets: List, 
                          file_path: Path, metadata) -> None:
        """
        Save triplets to CSV file with metadata header.
        
        Args:
            triplets: List of triplets to save
            file_path: Output file path
            metadata: Dataset metadata
        """
        logger.info(f"Saving {len(triplets)} triplets to {file_path}")
        
        with open(file_path, 'w', newline='') as csvfile:
            # Write metadata as comment
            metadata_json = json.dumps(metadata.__dict__, indent=None)
            csvfile.write(f"# Metadata: {metadata_json}\n")
            
            # Write CSV data
            fieldnames = [
                'triplet_id', 'anchor_seq', 'anchor_structure',
                'positive_seq', 'positive_structure', 'negative_seq', 'negative_structure',
                'total_modifications', 'stem_modifications', 'hloop_modifications',
                'iloop_modifications', 'bulge_modifications', 'mloop_modifications',
                'total_len_stems', 'total_len_hloops', 'total_len_iloops',
                'total_len_bulges', 'total_len_mloops',
                'stem_insertions', 'stem_deletions',
                'hloop_insertions', 'hloop_deletions',
                'iloop_insertions', 'iloop_deletions',
                'bulge_insertions', 'bulge_deletions',
                'mloop_insertions', 'mloop_deletions',
                'total_insertions', 'total_deletions'
            ]
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for triplet in triplets:
                writer.writerow({
                    'triplet_id': triplet.triplet_id,
                    'anchor_seq': triplet.anchor_seq,
                    'anchor_structure': triplet.anchor_structure,
                    'positive_seq': triplet.positive_seq,
                    'positive_structure': triplet.positive_structure,
                    'negative_seq': triplet.negative_seq,
                    'negative_structure': triplet.negative_structure,
                    'total_modifications': triplet.total_modifications,
                    'stem_modifications': triplet.stem_modifications,
                    'hloop_modifications': triplet.hloop_modifications,
                    'iloop_modifications': triplet.iloop_modifications,
                    'bulge_modifications': triplet.bulge_modifications,
                    'mloop_modifications': triplet.mloop_modifications,
                    'total_len_stems': triplet.total_len_stems,
                    'total_len_hloops': triplet.total_len_hloops,
                    'total_len_iloops': triplet.total_len_iloops,
                    'total_len_bulges': triplet.total_len_bulges,
                    'total_len_mloops': triplet.total_len_mloops,
                    'stem_insertions': triplet.stem_insertions,
                    'stem_deletions': triplet.stem_deletions,
                    'hloop_insertions': triplet.hloop_insertions,
                    'hloop_deletions': triplet.hloop_deletions,
                    'iloop_insertions': triplet.iloop_insertions,
                    'iloop_deletions': triplet.iloop_deletions,
                    'bulge_insertions': triplet.bulge_insertions,
                    'bulge_deletions': triplet.bulge_deletions,
                    'mloop_insertions': triplet.mloop_insertions,
                    'mloop_deletions': triplet.mloop_deletions,
                    'total_insertions': triplet.total_insertions,
                    'total_deletions': triplet.total_deletions
                })
    
    def _save_metadata(self, metadata, file_path: Path) -> None:
        """Save metadata to JSON file."""
        with open(file_path, 'w') as f:
            json.dump(metadata.__dict__, f, indent=2)
        
        logger.info(f"Metadata saved to {file_path}")


class DatasetAnalyzer:
    """Analyze and generate statistics for RNA triplet datasets."""
    
    @staticmethod
    def analyze_dataset(triplets: List) -> dict:
        """
        Analyze a dataset and return statistics.
        
        Args:
            triplets: List of RNA triplets
            
        Returns:
            Dictionary containing dataset statistics
        """
        if not triplets:
            return {}
        
        # Basic statistics
        num_triplets = len(triplets)
        
        # Sequence length statistics
        anchor_lengths = [len(t.anchor_seq) for t in triplets]
        positive_lengths = [len(t.positive_seq) for t in triplets]
        negative_lengths = [len(t.negative_seq) for t in triplets]
        
        # Modification statistics
        total_mods = [t.total_modifications for t in triplets]
        stem_mods = [t.stem_modifications for t in triplets]
        hloop_mods = [t.hloop_modifications for t in triplets]
        iloop_mods = [t.iloop_modifications for t in triplets]
        bulge_mods = [t.bulge_modifications for t in triplets]
        mloop_mods = [t.mloop_modifications for t in triplets]
        
        stats = {
            'num_triplets': num_triplets,
            'sequence_lengths': {
                'anchor': {
                    'mean': sum(anchor_lengths) / len(anchor_lengths),
                    'min': min(anchor_lengths),
                    'max': max(anchor_lengths)
                },
                'positive': {
                    'mean': sum(positive_lengths) / len(positive_lengths),
                    'min': min(positive_lengths),
                    'max': max(positive_lengths)
                },
                'negative': {
                    'mean': sum(negative_lengths) / len(negative_lengths),
                    'min': min(negative_lengths),
                    'max': max(negative_lengths)
                }
            },
            'modifications': {
                'total': {
                    'mean': sum(total_mods) / len(total_mods),
                    'min': min(total_mods),
                    'max': max(total_mods)
                },
                'stem': {
                    'mean': sum(stem_mods) / len(stem_mods),
                    'min': min(stem_mods),
                    'max': max(stem_mods)
                },
                'hairpin': {
                    'mean': sum(hloop_mods) / len(hloop_mods),
                    'min': min(hloop_mods),
                    'max': max(hloop_mods)
                },
                'internal': {
                    'mean': sum(iloop_mods) / len(iloop_mods),
                    'min': min(iloop_mods),
                    'max': max(iloop_mods)
                },
                'bulge': {
                    'mean': sum(bulge_mods) / len(bulge_mods),
                    'min': min(bulge_mods),
                    'max': max(bulge_mods)
                },
                'multi': {
                    'mean': sum(mloop_mods) / len(mloop_mods),
                    'min': min(mloop_mods),
                    'max': max(mloop_mods)
                }
            }
        }
        
        return stats
    
    @staticmethod
    def print_statistics(stats: dict) -> None:
        """Print dataset statistics in a readable format."""
        print("\n=== Dataset Statistics ===")
        print(f"Number of triplets: {stats['num_triplets']}")
        
        print("\nSequence Lengths:")
        for seq_type, length_stats in stats['sequence_lengths'].items():
            print(f"  {seq_type.title()}: mean={length_stats['mean']:.1f}, "
                  f"min={length_stats['min']}, max={length_stats['max']}")
        
        print("\nModifications:")
        for mod_type, mod_stats in stats['modifications'].items():
            print(f"  {mod_type.title()}: mean={mod_stats['mean']:.2f}, "
                  f"min={mod_stats['min']}, max={mod_stats['max']}")
        print("=" * 30)