"""
RNA modification engine for applying structural changes.
"""

import logging
import random
from typing import Dict, List, Tuple, Optional
from scipy.stats import truncnorm

from .models import (
    BulgeGraph, NodeType, ModificationType, ModificationCounts, 
    SampledModifications, classify_node
)

logger = logging.getLogger(__name__)


class ModificationEngine:
    """Handles applying modifications to RNA sequences and structures."""
    
    def __init__(self, args):
        """Initialize the modification engine with configuration."""
        self.args = args
        self.modification_counts = {}
        self.node_actions = {}
    
    def _format_bulge_graph_summary(self, bulge_graph: BulgeGraph, sequence: str, structure: str) -> str:
        """Format a detailed summary of the bulge graph state with nodes sorted by first coordinate."""
        node_mapping = bulge_graph.node_mapping()
        
        # Sort nodes by their first coordinate for consistent output
        sorted_nodes = []
        for node_name, coords in node_mapping.items():
            first_coord = coords[0] if coords else float('inf')
            sorted_nodes.append((first_coord, node_name, coords))
        
        sorted_nodes.sort(key=lambda x: x[0])
        
        summary_lines = [
            f"Sequence:  {sequence}",
            f"Structure: {structure}",
            f"BulgeGraph Summary ({len(node_mapping)} nodes):"
        ]
        
        for _, node_name, coords in sorted_nodes:
            node_type = classify_node(node_name, coords)
            mod_count = self.modification_counts.get(node_name, 0)
            summary_lines.append(f"  {node_name:<4}: {coords} → {node_type.value} (modified {mod_count} times)")
        
        return "\n".join(summary_lines)

    def sample_modifications(self) -> SampledModifications:
        """Sample modification counts from truncated normal distributions."""
        return SampledModifications(
            n_stem_indels=self._sample_truncated_normal(
                self.args.n_stem_indels_mean, self.args.n_stem_indels_sd,
                self.args.n_stem_indels_min, self.args.n_stem_indels_max
            ),
            n_hloop_indels=self._sample_truncated_normal(
                self.args.n_hloop_indels_mean, self.args.n_hloop_indels_sd,
                self.args.n_hloop_indels_min, self.args.n_hloop_indels_max
            ),
            n_iloop_indels=self._sample_truncated_normal(
                self.args.n_iloop_indels_mean, self.args.n_iloop_indels_sd,
                self.args.n_iloop_indels_min, self.args.n_iloop_indels_max
            ),
            n_bulge_indels=self._sample_truncated_normal(
                self.args.n_bulge_indels_mean, self.args.n_bulge_indels_sd,
                self.args.n_bulge_indels_min, self.args.n_bulge_indels_max
            ),
            n_mloop_indels=self._sample_truncated_normal(
                self.args.n_mloop_indels_mean, self.args.n_mloop_indels_sd,
                self.args.n_mloop_indels_min, self.args.n_mloop_indels_max
            )
        )
    
    def _sample_truncated_normal(self, mean: float, sd: float, min_val: int, max_val: int) -> int:
        """Sample from a truncated normal distribution."""
        if min_val >= max_val:
            return min_val
        
        a, b = (min_val - mean) / sd, (max_val - mean) / sd
        sample = truncnorm.rvs(a, b, loc=mean, scale=sd)
        return int(round(sample))
    
    def apply_modifications(self, sequence: str, structure: str, 
                          bulge_graph: BulgeGraph, 
                          sampled_mods: SampledModifications) -> Tuple[str, str, ModificationCounts]:
        """
        Apply modifications to create a positive sample.
        
        Args:
            sequence: Original RNA sequence
            structure: Original dot-bracket structure
            bulge_graph: Parsed bulge graph
            sampled_mods: Sampled modification counts
            
        Returns:
            Tuple of (modified_sequence, modified_structure, modification_counts)
        """
        pos_seq = sequence
        pos_struct = structure
        self.modification_counts = {}
        self.node_actions = {}
        
        # Normalize modification counts if enabled
        if self.args.mod_normalization:
            factor = max(1.0, len(sequence) / self.args.normalization_len)
            sampled_mods = SampledModifications(
                n_stem_indels=int(sampled_mods.n_stem_indels * factor),
                n_hloop_indels=int(sampled_mods.n_hloop_indels * factor),
                n_iloop_indels=int(sampled_mods.n_iloop_indels * factor),
                n_bulge_indels=int(sampled_mods.n_bulge_indels * factor),
                n_mloop_indels=int(sampled_mods.n_mloop_indels * factor)
            )
        
        logger.debug(f"Applying modifications: {sampled_mods}")
        
        # Show initial bulge graph state
        logger.debug("=== Initial bulge graph state ===")
        logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
        
        # Apply stem modifications
        logger.debug(f"=== Starting stem modifications ({sampled_mods.n_stem_indels} planned) ===")
        for i in range(sampled_mods.n_stem_indels):
            logger.debug(f"--- Stem modification {i + 1}/{sampled_mods.n_stem_indels} ---")
            logger.debug(f"Before stem modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
            
            new_seq, new_struct = self._modify_stems(pos_seq, pos_struct, bulge_graph)

            # Bulge graph updated in place
            pos_seq, pos_struct = new_seq, new_struct
            
            logger.debug(f"After stem modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
        
        # Apply loop modifications
        logger.debug(f"=== Starting hairpin loop modifications ({sampled_mods.n_hloop_indels} planned) ===")
        for i in range(sampled_mods.n_hloop_indels):
            logger.debug(f"--- Hairpin modification {i + 1}/{sampled_mods.n_hloop_indels} ---")
            logger.debug(f"Before hairpin modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
            
            new_seq, new_struct = self._modify_loops(pos_seq, pos_struct, bulge_graph, NodeType.HAIRPIN)

            # Bulge graph updated in place
            pos_seq, pos_struct = new_seq, new_struct
            
            logger.debug(f"After hairpin modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
        
        logger.debug(f"=== Starting internal loop modifications ({sampled_mods.n_iloop_indels} planned) ===")
        for i in range(sampled_mods.n_iloop_indels):
            logger.debug(f"--- Internal loop modification {i + 1}/{sampled_mods.n_iloop_indels} ---")
            logger.debug(f"Before internal loop modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
            
            new_seq, new_struct = self._modify_loops(pos_seq, pos_struct, bulge_graph, NodeType.INTERNAL)

            # Bulge graph updated in place
            pos_seq, pos_struct = new_seq, new_struct
            
            logger.debug(f"After internal loop modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
        
        logger.debug(f"=== Starting bulge modifications ({sampled_mods.n_bulge_indels} planned) ===")
        for i in range(sampled_mods.n_bulge_indels):
            logger.debug(f"--- Bulge modification {i + 1}/{sampled_mods.n_bulge_indels} ---")
            logger.debug(f"Before bulge modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
            
            new_seq, new_struct = self._modify_loops(pos_seq, pos_struct, bulge_graph, NodeType.BULGE)

            # Bulge graph updated in place
            pos_seq, pos_struct = new_seq, new_struct
            
            logger.debug(f"After bulge modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
        
        logger.debug(f"=== Starting multi loop modifications ({sampled_mods.n_mloop_indels} planned) ===")
        for i in range(sampled_mods.n_mloop_indels):
            logger.debug(f"--- Multi loop modification {i + 1}/{sampled_mods.n_mloop_indels} ---")
            logger.debug(f"Before multi loop modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
            
            new_seq, new_struct = self._modify_loops(pos_seq, pos_struct, bulge_graph, NodeType.MULTI)

            # Bulge graph updated in place
            pos_seq, pos_struct = new_seq, new_struct
            
            logger.debug(f"After multi loop modification {i + 1}:")
            logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
        
        # Show final bulge graph state
        logger.debug("=== Final bulge graph state ===")
        logger.debug(self._format_bulge_graph_summary(bulge_graph, pos_seq, pos_struct))
        logger.debug(f"Length changed from {len(sequence)} to {len(pos_seq)}")
        
        # Calculate final modification counts
        mod_counts = self._calculate_modification_counts(bulge_graph)
        
        logger.debug("=== Modification summary ===")
        logger.debug(f"Total modifications: {mod_counts.total}")
        logger.debug(f"Stem modifications: {mod_counts.stem}")
        logger.debug(f"Hairpin modifications: {mod_counts.hloop}")
        logger.debug(f"Internal loop modifications: {mod_counts.iloop}")
        logger.debug(f"Bulge modifications: {mod_counts.bulge}")
        logger.debug(f"Multi loop modifications: {mod_counts.mloop}")
        
        return pos_seq, pos_struct, mod_counts
    
    def _modify_stems(self, sequence: str, structure: str, bulge_graph: BulgeGraph) -> Tuple[str, str]:
        """Apply a single stem modification."""
        eligible_nodes = self._get_eligible_nodes(bulge_graph, NodeType.STEM, 
                                                 self.args.stem_min_size, self.args.stem_max_size,
                                                 self.args.same_stem_max_n_mod)
        
        if not eligible_nodes:
            logger.debug("No eligible stem nodes found for modification")
            return sequence, structure
        
        node_name, coords = random.choice(eligible_nodes)
        action = self._choose_action(node_name, len(coords), 
                                   self.args.stem_min_size, self.args.stem_max_size)
        
        logger.debug(f"Modifying stem node '{node_name}' with {len(coords)} coordinates: {coords}")
        logger.debug(f"Action: {action.value}")
        
        if action == ModificationType.INSERT:
            sequence, structure = self._insert_stem_pair(sequence, structure, node_name, coords, bulge_graph)
            logger.debug("Inserted complementary base pair in stem")
        elif action == ModificationType.DELETE:
            sequence, structure = self._delete_stem_pair(sequence, structure, node_name, coords, bulge_graph)
            logger.debug("Deleted complementary base pair from stem")
        
        self.modification_counts[node_name] = self.modification_counts.get(node_name, 0) + 1
        return sequence, structure
    
    def _modify_loops(self, sequence: str, structure: str, bulge_graph: BulgeGraph, 
                     node_type: NodeType) -> Tuple[str, str]:
        """Apply a single loop modification."""
        # Get parameters based on loop type
        if node_type == NodeType.HAIRPIN:
            min_size, max_size, max_mods = (self.args.hloop_min_size, 
                                          self.args.hloop_max_size, 
                                          self.args.same_hloop_max_n_mod)
        elif node_type == NodeType.INTERNAL:
            min_size, max_size, max_mods = (self.args.iloop_min_size,
                                          self.args.iloop_max_size,
                                          self.args.same_iloop_max_n_mod)
        elif node_type == NodeType.BULGE:
            min_size, max_size, max_mods = (self.args.bulge_min_size,
                                          self.args.bulge_max_size,
                                          self.args.same_bulge_max_n_mod)
        elif node_type == NodeType.MULTI:
            min_size, max_size, max_mods = (self.args.mloop_min_size,
                                          self.args.mloop_max_size,
                                          self.args.same_mloop_max_n_mod)
        else:
            return sequence, structure
        
        eligible_nodes = self._get_eligible_nodes(bulge_graph, node_type, min_size, max_size, max_mods)
        
        if not eligible_nodes:
            logger.debug(f"No eligible {node_type.value} nodes found for modification")
            return sequence, structure
        
        node_name, coords = random.choice(eligible_nodes)
        action = self._choose_action(node_name, len(coords), min_size, max_size)
        
        logger.debug(f"Modifying {node_type.value} node '{node_name}' with {len(coords)} coordinates: {coords}")
        logger.debug(f"Action: {action.value}")
        
        if action == ModificationType.INSERT:
            sequence, structure = self._insert_loop_base(sequence, structure, node_name, coords, bulge_graph)
            logger.debug(f"Inserted base in {node_type.value}")
        elif action == ModificationType.DELETE:
            sequence, structure = self._delete_loop_base(sequence, structure, node_name, coords, min_size, bulge_graph)
            logger.debug(f"Deleted base from {node_type.value}")
        
        self.modification_counts[node_name] = self.modification_counts.get(node_name, 0) + 1
        return sequence, structure
    
    def _get_eligible_nodes(self, bulge_graph: BulgeGraph, node_type: NodeType,
                           min_size: int, max_size: int, max_mods: int) -> List[Tuple[str, List[int]]]:
        """Get nodes eligible for modification."""
        eligible = []
        nodes_by_type = bulge_graph.get_nodes_by_type(node_type)
        
        for node_name, coords in nodes_by_type:
            # Skip nodes with empty coordinates - they cannot be modified
            if not coords or len(coords) == 0:
                logger.debug(f"Skipping node '{node_name}' with empty coordinates")
                continue
                
            # Check modification count limit
            if self.modification_counts.get(node_name, 0) >= max_mods:
                continue
            
            current_size = len(coords)
            
            # Check if node can be modified
            if self._can_modify_node(node_name, current_size, min_size, max_size):
                eligible.append((node_name, coords))
        
        return eligible
    
    def _can_modify_node(self, node_name: str, current_size: int, min_size: int, max_size: int) -> bool:
        """Check if a node can be modified."""
        if node_name in self.node_actions:
            # Action already committed
            action = self.node_actions[node_name]
            if action == ModificationType.INSERT:
                return current_size < max_size
            elif action == ModificationType.DELETE:
                return current_size > min_size
        else:
            # Can modify if insertion or deletion is possible
            return current_size < max_size or current_size > min_size
        
        return False
    
    def _choose_action(self, node_name: str, current_size: int, min_size: int, max_size: int) -> ModificationType:
        """Choose insertion or deletion for a node."""
        if node_name in self.node_actions:
            return self.node_actions[node_name]
        
        can_insert = current_size < max_size
        can_delete = current_size > min_size
        
        if can_insert and can_delete:
            action = random.choice([ModificationType.INSERT, ModificationType.DELETE])
        elif can_insert:
            action = ModificationType.INSERT
        elif can_delete:
            action = ModificationType.DELETE
        else:
            action = ModificationType.INSERT  # Fallback
        
        self.node_actions[node_name] = action
        return action
    
    def _insert_stem_pair(self, sequence: str, structure: str, node_name: str,
                          coords: List[int], bulge_graph: BulgeGraph) -> Tuple[str, str]:
        """Insert a complementary base pair in a stem at random positions."""
        if len(coords) < 2:
            return sequence, structure
        
        sorted_coords = sorted(coords)
        n = len(sorted_coords)
        
        # Split into left and right halves
        left_half = sorted_coords[:n//2]
        right_half = sorted_coords[n//2:]
        
        # Choose a random position in the left half to insert after
        # We can insert between any adjacent positions in the left half
        left_insert_options = []
        for i in range(len(left_half)):
            if i == 0:
                # Can insert before the first position
                left_insert_options.append(left_half[i] - 1)
            # Can insert after each position
            left_insert_options.append(left_half[i])
        
        # Choose random insertion position in left half
        left_pos = random.choice(left_insert_options)
        
        # Calculate corresponding position in right half
        # If we insert at position left_pos, we need to find the corresponding right position
        # The pairing logic: if left_half[i] pairs with right_half[-(i+1)]
        # We need to find which left position corresponds to our insertion
        
        # Find which "gap" in the left half we're inserting into
        left_gap_index = 0
        for i, pos in enumerate(left_half):
            if left_pos <= pos:
                left_gap_index = i
                break
        else:
            left_gap_index = len(left_half)
        
        # The corresponding right position is the mirror of the left gap
        if left_gap_index == 0:
            # Inserting before first left base, so insert after last right base
            right_pos = right_half[-1]
        elif left_gap_index == len(left_half):
            # Inserting after last left base, so insert before first right base
            right_pos = right_half[0] - 1
        else:
            # Inserting between left positions, insert between corresponding right positions
            # If inserting between left_half[i-1] and left_half[i], 
            # insert between right_half[-(i)] and right_half[-(i+1)]
            right_mirror_index = len(right_half) - left_gap_index
            if right_mirror_index > 0:
                right_pos = right_half[right_mirror_index - 1]
            else:
                right_pos = right_half[0] - 1
        
        # Choose complementary bases
        complement_pairs = [('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')]
        base_left, base_right = random.choice(complement_pairs)
        
        logger.debug(f"Inserting stem pair at positions {left_pos} ('{base_left}') and {right_pos} ('{base_right}')")
        
        # Insert bases (higher index first to avoid position shifts)
        if left_pos < right_pos:
            sequence = sequence[:right_pos] + base_right + sequence[right_pos:]
            structure = structure[:right_pos] + ')' + structure[right_pos:]
            sequence = sequence[:left_pos] + base_left + sequence[left_pos:]
            structure = structure[:left_pos] + '(' + structure[left_pos:]
        else:
            sequence = sequence[:left_pos] + base_left + sequence[left_pos:]
            structure = structure[:left_pos] + '(' + structure[left_pos:]
            sequence = sequence[:right_pos] + base_right + sequence[right_pos:]
            structure = structure[:right_pos] + ')' + structure[right_pos:]

        # Update bulge graph indices
        from .bulge_graph_updater import BulgeGraphUpdater
        BulgeGraphUpdater.insert_stem_pair(bulge_graph, node_name, left_pos, right_pos)

        return sequence, structure
    
    def _delete_stem_pair(self, sequence: str, structure: str, node_name: str,
                          coords: List[int], bulge_graph: BulgeGraph) -> Tuple[str, str]:
        """Delete a complementary base pair from a stem at random positions."""
        if len(coords) < 4:  # Need at least 4 bases to safely delete a pair
            return sequence, structure
        
        sorted_coords = sorted(coords)
        n = len(sorted_coords)
        
        # Split into left and right halves
        left_half = sorted_coords[:n//2]
        right_half = sorted_coords[n//2:]
        
        # Choose a random position from the left half to delete
        left_delete_idx = random.randint(0, len(left_half) - 1)
        left_pos = left_half[left_delete_idx] - 1  # Convert to 0-based
        
        # Calculate the corresponding position in the right half
        # The pairing logic: left_half[i] pairs with right_half[-(i+1)]
        # So left_half[left_delete_idx] pairs with right_half[-(left_delete_idx+1)]
        right_delete_idx = len(right_half) - 1 - left_delete_idx
        right_pos = right_half[right_delete_idx] - 1  # Convert to 0-based
        
        logger.debug(f"Deleting stem pair at positions {left_pos} and {right_pos}")
        logger.debug(f"Left half: {left_half}, deleting index {left_delete_idx}")
        logger.debug(f"Right half: {right_half}, deleting index {right_delete_idx}")
        
        # Delete bases (higher index first to avoid position shifts)
        if right_pos < len(sequence) and left_pos < len(sequence) and right_pos > left_pos:
            sequence = sequence[:right_pos] + sequence[right_pos+1:]
            structure = structure[:right_pos] + structure[right_pos+1:]
            sequence = sequence[:left_pos] + sequence[left_pos+1:]
            structure = structure[:left_pos] + structure[left_pos+1:]
        elif left_pos < len(sequence) and right_pos < len(sequence) and left_pos > right_pos:
            # Handle case where positions might be reversed
            sequence = sequence[:left_pos] + sequence[left_pos+1:]
            structure = structure[:left_pos] + structure[left_pos+1:]
            sequence = sequence[:right_pos] + sequence[right_pos+1:]
            structure = structure[:right_pos] + structure[right_pos+1:]

        # Update bulge graph indices
        from .bulge_graph_updater import BulgeGraphUpdater
        BulgeGraphUpdater.delete_stem_pair(bulge_graph, node_name, left_pos, right_pos)

        return sequence, structure
    
    def _insert_loop_base(self, sequence: str, structure: str, node_name: str,
                          coords: List[int], bulge_graph: BulgeGraph) -> Tuple[str, str]:
        """Insert a base in a loop region."""
        if not coords:
            return sequence, structure
        
        # Choose random position within the loop (convert to 0-based)
        pos = random.choice(coords) - 1
        
        # Insert random base
        base = random.choice(['A', 'C', 'G', 'U'])
        sequence = sequence[:pos] + base + sequence[pos:]
        structure = structure[:pos] + '.' + structure[pos:]

        from .bulge_graph_updater import BulgeGraphUpdater
        BulgeGraphUpdater.insert_loop_base(bulge_graph, node_name, pos)

        return sequence, structure
    
    def _delete_loop_base(self, sequence: str, structure: str, node_name: str,
                          coords: List[int], min_size: int, bulge_graph: BulgeGraph) -> Tuple[str, str]:
        """Delete a base from a loop region."""
        if len(coords) <= min_size or not coords:
            return sequence, structure
        
        # Choose random position to delete (convert to 0-based)
        pos = random.choice(coords) - 1
        
        if pos < len(sequence):
            sequence = sequence[:pos] + sequence[pos+1:]
            structure = structure[:pos] + structure[pos+1:]
            from .bulge_graph_updater import BulgeGraphUpdater
            BulgeGraphUpdater.delete_loop_base(bulge_graph, node_name, pos)

        return sequence, structure
    
    def _calculate_modification_counts(self, bulge_graph: BulgeGraph) -> ModificationCounts:
        """Calculate modification counts by structural element type."""
        counts = ModificationCounts()
        node_mapping = bulge_graph.node_mapping()
        
        for node_name, count in self.modification_counts.items():
            if node_name in node_mapping:
                coords = node_mapping[node_name]
                node_type = classify_node(node_name, coords)
                
                if node_type == NodeType.STEM:
                    counts.stem += count
                elif node_type == NodeType.HAIRPIN:
                    counts.hloop += count
                elif node_type == NodeType.INTERNAL:
                    counts.iloop += count
                elif node_type == NodeType.BULGE:
                    counts.bulge += count
                elif node_type == NodeType.MULTI:
                    counts.mloop += count
        
        counts.total = sum(self.modification_counts.values())
        return counts