#!/usr/bin/env python3
"""
Structure plotting utilities for RNA sequences and structures.

This module provides functionality to visualize RNA secondary structures
using matplotlib, similar to the rnartistcore functionality from the Rust version.
"""

import logging
import math
import random
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
import numpy as np

logger = logging.getLogger(__name__)


class RnaStructurePlotter:
    """Plots RNA secondary structures using matplotlib."""
    
    def __init__(self, figsize: Tuple[int, int] = (12, 8)):
        """Initialize the plotter with figure settings."""
        self.figsize = figsize
        self.base_colors = {
            'A': '#FF6B6B',  # Red
            'U': '#4ECDC4',  # Teal
            'G': '#45B7D1',  # Blue
            'C': '#FFA07A',  # Orange
            'N': '#808080'   # Gray for unknown
        }
        self.bond_color = '#2C3E50'  # Dark blue-gray
        self.backbone_color = '#34495E'  # Darker gray
        
    def plot_triplet(self, triplet, output_path: Path, title: str = None) -> None:
        """
        Plot an RNA triplet (anchor, positive, negative) in a single figure.
        
        Args:
            triplet: RnaTriplet object containing sequences and structures
            output_path: Path to save the plot
            title: Optional title for the plot
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        # Plot each structure
        self._plot_structure(axes[0], triplet.anchor_seq, triplet.anchor_structure, 
                           f"Anchor (ID: {triplet.triplet_id})")
        self._plot_structure(axes[1], triplet.positive_seq, triplet.positive_structure, 
                           f"Positive ({triplet.total_modifications} mods)")
        self._plot_structure(axes[2], triplet.negative_seq, triplet.negative_structure, 
                           "Negative")
        
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.debug(f"Saved triplet plot to {output_path}")
    
    def plot_structure(self, sequence: str, structure: str, output_path: Path, 
                      title: str = None) -> None:
        """
        Plot a single RNA structure.
        
        Args:
            sequence: RNA sequence string
            structure: Dot-bracket structure string
            output_path: Path to save the plot
            title: Optional title for the plot
        """
        fig, ax = plt.subplots(1, 1, figsize=self.figsize)
        self._plot_structure(ax, sequence, structure, title)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.debug(f"Saved structure plot to {output_path}")
    
    def _plot_structure(self, ax, sequence: str, structure: str, title: str = None) -> None:
        """
        Plot RNA structure on given axes using circular layout.
        
        Args:
            ax: Matplotlib axes object
            sequence: RNA sequence string
            structure: Dot-bracket structure string
            title: Optional title for the subplot
        """
        if len(sequence) != len(structure):
            logger.error(f"Sequence and structure length mismatch: {len(sequence)} vs {len(structure)}")
            return
        
        n = len(sequence)
        if n == 0:
            ax.text(0.5, 0.5, "Empty sequence", ha='center', va='center', transform=ax.transAxes)
            if title:
                ax.set_title(title)
            return
        
        # Calculate circular positions
        positions = self._calculate_circular_positions(n)
        
        # Find base pairs
        pairs = self._find_base_pairs(structure)
        
        # Draw backbone
        self._draw_backbone(ax, positions)
        
        # Draw base pairs
        self._draw_base_pairs(ax, positions, pairs)
        
        # Draw bases
        self._draw_bases(ax, positions, sequence)
        
        # Set up axes
        ax.set_aspect('equal')
        ax.axis('off')
        
        if title:
            ax.set_title(title, fontsize=12, fontweight='bold')
        
        # Add sequence info as text
        info_text = f"Length: {n} nt"
        if pairs:
            info_text += f", {len(pairs)} pairs"
        ax.text(0.02, 0.02, info_text, transform=ax.transAxes, 
               fontsize=8, alpha=0.7)
    
    def _calculate_circular_positions(self, n: int, radius: float = 1.0) -> List[Tuple[float, float]]:
        """Calculate positions for bases in a circular layout."""
        positions = []
        for i in range(n):
            angle = 2 * math.pi * i / n - math.pi / 2  # Start from top
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            positions.append((x, y))
        return positions
    
    def _find_base_pairs(self, structure: str) -> List[Tuple[int, int]]:
        """Find base pairs from dot-bracket notation."""
        pairs = []
        stack = []
        
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    j = stack.pop()
                    pairs.append((j, i))
        
        return pairs
    
    def _draw_backbone(self, ax, positions: List[Tuple[float, float]]) -> None:
        """Draw the RNA backbone connecting consecutive bases."""
        if len(positions) < 2:
            return
        
        # Create line segments for backbone
        segments = []
        for i in range(len(positions) - 1):
            segments.append([positions[i], positions[i + 1]])
        
        lc = LineCollection(segments, colors=self.backbone_color, linewidths=2, alpha=0.6)
        ax.add_collection(lc)
    
    def _draw_base_pairs(self, ax, positions: List[Tuple[float, float]], 
                        pairs: List[Tuple[int, int]]) -> None:
        """Draw lines representing base pairs."""
        for i, j in pairs:
            if i < len(positions) and j < len(positions):
                x1, y1 = positions[i]
                x2, y2 = positions[j]
                ax.plot([x1, x2], [y1, y2], color=self.bond_color, 
                       linewidth=1.5, alpha=0.8, zorder=1)
    
    def _draw_bases(self, ax, positions: List[Tuple[float, float]], sequence: str) -> None:
        """Draw individual bases as colored circles with letters."""
        for i, (x, y) in enumerate(positions):
            if i < len(sequence):
                base = sequence[i].upper()
                color = self.base_colors.get(base, self.base_colors['N'])
                
                # Draw circle
                circle = patches.Circle((x, y), 0.05, facecolor=color, 
                                      edgecolor='white', linewidth=1, zorder=2)
                ax.add_patch(circle)
                
                # Draw letter
                ax.text(x, y, base, ha='center', va='center', 
                       fontsize=8, fontweight='bold', color='white', zorder=3)


class StructurePlotManager:
    """Manages plotting of RNA structures for the dataset generator."""
    
    def __init__(self, args):
        """Initialize with command line arguments."""
        self.args = args
        self.plotter = RnaStructurePlotter()
        self.plot_dir = None
        
    def setup_plotting(self, output_dir: Path) -> None:
        """Setup the plotting directory."""
        if self.args.plot:
            self.plot_dir = output_dir / "plots"
            self.plot_dir.mkdir(exist_ok=True)
            logger.info(f"Plot directory created: {self.plot_dir}")
    
    def plot_triplets(self, triplets: List, num_plots: Optional[int] = None) -> None:
        """
        Plot a selection of triplets from the dataset.
        
        Args:
            triplets: List of RnaTriplet objects
            num_plots: Number of triplets to plot (defaults to args.num_plots)
        """
        if not self.args.plot or not self.plot_dir:
            return
        
        if not triplets:
            logger.warning("No triplets to plot")
            return
        
        num_plots = num_plots or self.args.num_plots
        num_plots = min(num_plots, len(triplets))
        
        logger.info(f"Plotting {num_plots} triplets...")
        
        # Select triplets to plot (random sampling if more than requested)
        if len(triplets) <= num_plots:
            selected_triplets = triplets
        else:
            selected_triplets = random.sample(triplets, num_plots)
        
        for i, triplet in enumerate(selected_triplets):
            output_path = self.plot_dir / f"triplet_{triplet.triplet_id:04d}.png"
            title = f"RNA Triplet {triplet.triplet_id}"
            
            try:
                self.plotter.plot_triplet(triplet, output_path, title)
                
                if (i + 1) % 10 == 0 or i == len(selected_triplets) - 1:
                    logger.info(f"Plotted {i + 1}/{len(selected_triplets)} triplets")
                    
            except Exception as e:
                logger.error(f"Failed to plot triplet {triplet.triplet_id}: {e}")
                continue
        
        logger.info(f"Plotting completed. Plots saved to: {self.plot_dir}")
    
    def plot_sample_structures(self, triplets: List, sample_size: int = 5) -> None:
        """
        Plot individual structures from a sample of triplets.
        
        Args:
            triplets: List of RnaTriplet objects
            sample_size: Number of sample structures to plot
        """
        if not self.args.plot or not self.plot_dir:
            return
        
        if not triplets:
            return
        
        sample_size = min(sample_size, len(triplets))
        sampled_triplets = random.sample(triplets, sample_size)
        
        structures_dir = self.plot_dir / "individual_structures"
        structures_dir.mkdir(exist_ok=True)
        
        for triplet in sampled_triplets:
            base_name = f"triplet_{triplet.triplet_id:04d}"
            
            # Plot anchor
            self.plotter.plot_structure(
                triplet.anchor_seq, triplet.anchor_structure,
                structures_dir / f"{base_name}_anchor.png",
                f"Anchor - Triplet {triplet.triplet_id}"
            )
            
            # Plot positive
            self.plotter.plot_structure(
                triplet.positive_seq, triplet.positive_structure,
                structures_dir / f"{base_name}_positive.png",
                f"Positive - Triplet {triplet.triplet_id} ({triplet.total_modifications} mods)"
            )
            
            # Plot negative
            self.plotter.plot_structure(
                triplet.negative_seq, triplet.negative_structure,
                structures_dir / f"{base_name}_negative.png",
                f"Negative - Triplet {triplet.triplet_id}"
            )
        
        logger.info(f"Individual structure plots saved to: {structures_dir}")