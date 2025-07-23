"""
Main dataset generator that orchestrates the creation of RNA triplets.
"""

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List
from tqdm import tqdm

from .models import RnaTriplet, DatasetMetadata
from .rna_generator import RnaGenerator, BulgeGraphParser
from .modification_engine import ModificationEngine
from .negative_generator import NegativeSampleGenerator, AppendingEngine
from .models import NodeType, classify_node

logger = logging.getLogger(__name__)


class DatasetGenerator:
    """Main class for generating RNA triplet datasets."""
    
    def __init__(self, args):
        """Initialize the dataset generator with configuration."""
        self.args = args
        self.rna_generator = RnaGenerator()
        self.bulge_parser = BulgeGraphParser()
        self.modification_engine = ModificationEngine(args)
        self.negative_generator = NegativeSampleGenerator(args)
        self.appending_engine = AppendingEngine(args)
    
    def generate_dataset(self) -> List[RnaTriplet]:
        """
        Generate the complete dataset of RNA triplets.
        
        Returns:
            List of RNA triplets
        """
        logger.info(f"Generating {self.args.num_structures} triplets using {self.args.num_workers} workers")
        
        if self.args.num_workers == 1:
            # Single-threaded generation
            return self._generate_sequential()
        else:
            # Multi-threaded generation
            return self._generate_parallel()
    
    def _generate_sequential(self) -> List[RnaTriplet]:
        """Generate triplets sequentially."""
        triplets = []
        
        with tqdm(total=self.args.num_structures, desc="Generating triplets") as pbar:
            for i in range(self.args.num_structures):
                triplet = self._generate_single_triplet(i)
                triplets.append(triplet)
                pbar.update(1)
        
        return triplets
    
    def _generate_parallel(self) -> List[RnaTriplet]:
        """Generate triplets in parallel using multiprocessing."""
        triplets = []
        
        # Calculate work distribution
        total_work = self.args.num_structures
        batch_size = self.args.batch_size
        num_batches = (total_work + batch_size - 1) // batch_size
        
        with tqdm(total=total_work, desc="Generating triplets") as pbar:
            with ProcessPoolExecutor(max_workers=self.args.num_workers) as executor:
                # Submit batch jobs
                future_to_batch = {}
                
                for batch_idx in range(num_batches):
                    start_idx = batch_idx * batch_size
                    end_idx = min(start_idx + batch_size, total_work)
                    batch_indices = list(range(start_idx, end_idx))
                    
                    future = executor.submit(self._generate_batch, batch_indices)
                    future_to_batch[future] = batch_idx
                
                # Collect results
                for future in as_completed(future_to_batch):
                    batch_triplets = future.result()
                    triplets.extend(batch_triplets)
                    pbar.update(len(batch_triplets))
        
        # Sort by triplet_id to maintain order
        triplets.sort(key=lambda x: x.triplet_id)
        return triplets
    
    def _generate_batch(self, indices: List[int]) -> List[RnaTriplet]:
        """Generate a batch of triplets (used in multiprocessing)."""
        batch_triplets = []
        
        for idx in indices:
            triplet = self._generate_single_triplet(idx)
            batch_triplets.append(triplet)

        return batch_triplets

    def _calculate_structure_lengths(self, bulge_graph) -> tuple[int, int, int, int, int]:
        """Calculate total nucleotide lengths for each structure type on the anchor."""
        totals = {
            NodeType.STEM: 0,
            NodeType.HAIRPIN: 0,
            NodeType.INTERNAL: 0,
            NodeType.BULGE: 0,
            NodeType.MULTI: 0,
        }

        for node_name, coords in bulge_graph.node_mapping().items():
            ntype = classify_node(node_name, coords)
            if not coords:
                size = 0
            elif ntype == NodeType.STEM:
                size = coords[1] - coords[0] + 1 if len(coords) >= 2 else 0
            elif ntype == NodeType.HAIRPIN:
                size = coords[1] - coords[0] + 1 if len(coords) >= 2 else 0
            elif ntype == NodeType.BULGE:
                size = coords[1] - coords[0] + 1 if len(coords) >= 2 else 0
            elif ntype == NodeType.MULTI:
                size = coords[1] - coords[0] + 1 if len(coords) == 2 else 0
            elif ntype == NodeType.INTERNAL:
                size = (
                    coords[1] - coords[0] + coords[3] - coords[2] + 2
                    if len(coords) >= 4
                    else 0
                )
            else:
                size = 0
            if ntype in totals:
                totals[ntype] += size

        return (
            totals[NodeType.STEM],
            totals[NodeType.HAIRPIN],
            totals[NodeType.INTERNAL],
            totals[NodeType.BULGE],
            totals[NodeType.MULTI],
        )
    
    def _generate_single_triplet(self, triplet_id: int) -> RnaTriplet:
        """
        Generate a single RNA triplet.
        
        Args:
            triplet_id: Unique identifier for the triplet
            
        Returns:
            Complete RNA triplet
        """
        # Step 1: Generate anchor sequence
        anchor_seq, anchor_struct = self._generate_anchor()
        
        # Step 2: Parse anchor structure into bulge graph
        anchor_graph = self.bulge_parser.parse_structure(anchor_struct)

        stem_len, hloop_len, iloop_len, bulge_len, mloop_len = self._calculate_structure_lengths(anchor_graph)
        
        # Step 3: Sample modifications and generate positive
        sampled_mods = self.modification_engine.sample_modifications(anchor_graph)
        pos_seq, pos_struct, mod_counts, action_counts = self.modification_engine.apply_modifications(
            anchor_seq, anchor_struct, anchor_graph, sampled_mods
        )
        
        # Step 4: Generate negative sample
        neg_seq, neg_struct = self.negative_generator.generate_negative_sample(anchor_seq)
        
        # Step 5: Maybe apply appending events
        pos_seq, pos_struct, neg_seq, neg_struct = self.appending_engine.maybe_append_sequences(
            pos_seq, pos_struct, neg_seq, neg_struct, len(anchor_seq)
        )
        
        # Step 6: Create and return the triplet
        return RnaTriplet(
            triplet_id=triplet_id,
            anchor_seq=anchor_seq,
            anchor_structure=anchor_struct,
            positive_seq=pos_seq,
            positive_structure=pos_struct,
            negative_seq=neg_seq,
            negative_structure=neg_struct,
            total_modifications=mod_counts.total,
            stem_modifications=mod_counts.stem,
            hloop_modifications=mod_counts.hloop,
            iloop_modifications=mod_counts.iloop,
            bulge_modifications=mod_counts.bulge,
            mloop_modifications=mod_counts.mloop,
            total_len_stems=stem_len,
            total_len_hloops=hloop_len,
            total_len_iloops=iloop_len,
            total_len_bulges=bulge_len,
            total_len_mloops=mloop_len,
            stem_insertions=action_counts.stem_insertions,
            stem_deletions=action_counts.stem_deletions,
            hloop_insertions=action_counts.hloop_insertions,
            hloop_deletions=action_counts.hloop_deletions,
            iloop_insertions=action_counts.iloop_insertions,
            iloop_deletions=action_counts.iloop_deletions,
            bulge_insertions=action_counts.bulge_insertions,
            bulge_deletions=action_counts.bulge_deletions,
            mloop_insertions=action_counts.mloop_insertions,
            mloop_deletions=action_counts.mloop_deletions,
            total_insertions=action_counts.total_insertions,
            total_deletions=action_counts.total_deletions
        )
    
    def _generate_anchor(self) -> tuple[str, str]:
        """
        Generate an anchor sequence and its structure.
        
        Returns:
            Tuple of (sequence, structure)
        """
        # Choose sequence length
        length = self.rna_generator.choose_sequence_length(
            self.args.seq_len_distribution,
            self.args.seq_min_len,
            self.args.seq_max_len,
            self.args.seq_len_mean,
            self.args.seq_len_sd
        )
        
        # Generate random sequence
        sequence = self.rna_generator.generate_random_sequence(length)
        
        # Fold sequence to get structure
        structure = self.rna_generator.fold_rna(sequence)
        
        return sequence, structure