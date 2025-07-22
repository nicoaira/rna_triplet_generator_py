#!/usr/bin/env python3
"""
Structure plotting utilities for RNA sequences and structures using rnartistcore.

This module provides functionality to visualize RNA secondary structures
using rnartistcore, similar to the functionality from the Rust version.
"""

import logging
import random
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional
import shutil

logger = logging.getLogger(__name__)


class RnartistCorePlotter:
    """Plots RNA secondary structures using rnartistcore."""
    
    def __init__(self):
        """Initialize the plotter and check for rnartistcore availability."""
        self.rnartistcore_available = self._check_rnartistcore()
        self.montage_available = self._check_montage()
        
    def _check_rnartistcore(self) -> bool:
        """Check if rnartistcore is available in the system."""
        try:
            result = subprocess.run(['rnartistcore', '--version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info("rnartistcore found and available")
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        logger.warning("rnartistcore not found. Structure plotting will be disabled.")
        return False
    
    def _check_montage(self) -> bool:
        """Check if ImageMagick montage is available for combining images."""
        try:
            result = subprocess.run(['montage', '-version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info("ImageMagick montage found and available")
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        logger.warning("ImageMagick montage not found. Will create individual plots only.")
        return False
    
    def plot_triplet(self, triplet, output_path: Path, title: str = None) -> None:
        """
        Plot an RNA triplet (anchor, positive, negative) using rnartistcore.
        
        Args:
            triplet: RnaTriplet object containing sequences and structures
            output_path: Path to save the plot
            title: Optional title for the plot (not used in rnartistcore)
        """
        if not self.rnartistcore_available:
            logger.error("Cannot plot: rnartistcore not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Generate individual structure plots
            anchor_png = temp_path / f"triplet_{triplet.triplet_id}_anchor.png"
            positive_png = temp_path / f"triplet_{triplet.triplet_id}_positive.png"
            negative_png = temp_path / f"triplet_{triplet.triplet_id}_negative.png"
            
            # Create rnartistcore script for all three structures
            batch_script = self._create_triplet_batch_script(
                triplet, temp_path.resolve()
            )
            
            # Write batch script
            script_path = temp_path / "batch_plot.kts"
            with open(script_path, 'w') as f:
                f.write(batch_script)
            
            # Run rnartistcore
            try:
                result = subprocess.run(['rnartistcore', str(script_path)], 
                                      capture_output=True, text=True, timeout=60)
                if result.returncode != 0:
                    logger.error(f"rnartistcore failed: {result.stderr}")
                    return
            except subprocess.TimeoutExpired:
                logger.error("rnartistcore timed out")
                return
            
            # Check if all individual images were created
            if not (anchor_png.exists() and positive_png.exists() and negative_png.exists()):
                logger.error(f"One or more images for triplet {triplet.triplet_id} not found")
                return
            
            # Combine images using montage if available
            if self.montage_available:
                self._combine_images_with_montage([anchor_png, positive_png, negative_png], 
                                                output_path)
            else:
                # If montage not available, just copy the anchor image as fallback
                shutil.copy2(anchor_png, output_path)
                logger.warning(f"Montage not available. Saved only anchor for triplet {triplet.triplet_id}")
        
        logger.debug(f"Saved triplet plot to {output_path}")
    
    def plot_structure(self, sequence: str, structure: str, output_path: Path, 
                      title: str = None, structure_name: str = "structure") -> None:
        """
        Plot a single RNA structure using rnartistcore.
        
        Args:
            sequence: RNA sequence string
            structure: Dot-bracket structure string
            output_path: Path to save the plot
            title: Optional title for the plot (not used in rnartistcore)
            structure_name: Name for the structure in rnartistcore
        """
        if not self.rnartistcore_available:
            logger.error("Cannot plot: rnartistcore not available")
            return
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create rnartistcore script
            script_content = self._create_single_structure_script(
                sequence, structure, temp_path.resolve(), structure_name
            )
            
            # Write script
            script_path = temp_path / "plot.kts"
            with open(script_path, 'w') as f:
                f.write(script_content)
            
            # Run rnartistcore
            try:
                result = subprocess.run(['rnartistcore', str(script_path)], 
                                      capture_output=True, text=True, timeout=30)
                if result.returncode != 0:
                    logger.error(f"rnartistcore failed: {result.stderr}")
                    return
            except subprocess.TimeoutExpired:
                logger.error("rnartistcore timed out")
                return
            
            # Copy result to output path
            temp_png = temp_path / f"{structure_name}.png"
            if temp_png.exists():
                shutil.copy2(temp_png, output_path)
            else:
                logger.error(f"rnartistcore did not create expected output: {temp_png}")
        
        logger.debug(f"Saved structure plot to {output_path}")
    
    def _create_triplet_batch_script(self, triplet, temp_dir: Path) -> str:
        """Create rnartistcore batch script for a triplet."""
        anchor_block = f'''rnartist {{
  png {{
    path = "{temp_dir}"
    width = 500.0
    height = 500.0
  }}
  ss {{
    bn {{
      value = "{triplet.anchor_structure}"
      seq = "{triplet.anchor_seq}"
      name = "triplet_{triplet.triplet_id}_anchor"
    }}
  }}
  theme {{
    details {{
      value = 5
    }}
    scheme {{
      value = "Pumpkin Vegas"
    }}
  }}
}}'''
        
        positive_block = f'''rnartist {{
  png {{
    path = "{temp_dir}"
    width = 500.0
    height = 500.0
  }}
  ss {{
    bn {{
      value = "{triplet.positive_structure}"
      seq = "{triplet.positive_seq}"
      name = "triplet_{triplet.triplet_id}_positive"
    }}
  }}
  theme {{
    details {{
      value = 5
    }}
    scheme {{
      value = "Pumpkin Vegas"
    }}
  }}
}}'''
        
        negative_block = f'''rnartist {{
  png {{
    path = "{temp_dir}"
    width = 500.0
    height = 500.0
  }}
  ss {{
    bn {{
      value = "{triplet.negative_structure}"
      seq = "{triplet.negative_seq}"
      name = "triplet_{triplet.triplet_id}_negative"
    }}
  }}
  theme {{
    details {{
      value = 5
    }}
    scheme {{
      value = "Pumpkin Vegas"
    }}
  }}
}}'''
        
        return f"{anchor_block}\n{positive_block}\n{negative_block}\n"
    
    def _create_single_structure_script(self, sequence: str, structure: str, 
                                      temp_dir: Path, name: str) -> str:
        """Create rnartistcore script for a single structure."""
        return f'''rnartist {{
  png {{
    path = "{temp_dir}"
    width = 500.0
    height = 500.0
  }}
  ss {{
    bn {{
      value = "{structure}"
      seq = "{sequence}"
      name = "{name}"
    }}
  }}
  theme {{
    details {{
      value = 5
    }}
    scheme {{
      value = "Pumpkin Vegas"
    }}
  }}
}}'''
    
    def _combine_images_with_montage(self, image_paths: List[Path], output_path: Path) -> None:
        """Combine multiple images into one using ImageMagick montage."""
        try:
            cmd = ['montage'] + [str(p) for p in image_paths] + [
                '-tile', '3x1',
                '-geometry', '+10+10',
                str(output_path)
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.returncode != 0:
                logger.error(f"montage failed: {result.stderr}")
        except subprocess.TimeoutExpired:
            logger.error("montage timed out")


class StructurePlotManager:
    """Manages plotting of RNA structures for the dataset generator."""
    
    def __init__(self, args):
        """Initialize with command line arguments."""
        self.args = args
        self.plotter = RnartistCorePlotter()
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
        
        if not self.plotter.rnartistcore_available:
            logger.error("Cannot plot: rnartistcore not available")
            return
        
        if not triplets:
            logger.warning("No triplets to plot")
            return
        
        num_plots = num_plots or self.args.num_plots
        num_plots = min(num_plots, len(triplets))
        
        logger.info(f"Plotting {num_plots} triplets using rnartistcore...")
        
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
                
                if (i + 1) % 5 == 0 or i == len(selected_triplets) - 1:
                    logger.info(f"Plotted {i + 1}/{len(selected_triplets)} triplets")
                    
            except Exception as e:
                logger.error(f"Failed to plot triplet {triplet.triplet_id}: {e}")
                continue
        
        logger.info(f"Plotting completed. Plots saved to: {self.plot_dir}")
    
    def plot_sample_structures(self, triplets: List, sample_size: int = 5) -> None:
        """
        Plot individual structures from a sample of triplets using rnartistcore.
        
        Args:
            triplets: List of RnaTriplet objects
            sample_size: Number of sample structures to plot
        """
        if not self.args.plot or not self.plot_dir:
            return
        
        if not self.plotter.rnartistcore_available:
            logger.error("Cannot plot individual structures: rnartistcore not available")
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
                f"Anchor - Triplet {triplet.triplet_id}",
                f"{base_name}_anchor"
            )
            
            # Plot positive
            self.plotter.plot_structure(
                triplet.positive_seq, triplet.positive_structure,
                structures_dir / f"{base_name}_positive.png",
                f"Positive - Triplet {triplet.triplet_id} ({triplet.total_modifications} mods)",
                f"{base_name}_positive"
            )
            
            # Plot negative
            self.plotter.plot_structure(
                triplet.negative_seq, triplet.negative_structure,
                structures_dir / f"{base_name}_negative.png",
                f"Negative - Triplet {triplet.triplet_id}",
                f"{base_name}_negative"
            )
        
        logger.info(f"Individual structure plots saved to: {structures_dir}")