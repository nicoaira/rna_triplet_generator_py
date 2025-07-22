# RNA Triplet Generator - Python Version

A modular Python implementation of the RNA triplet generator for creating machine learning datasets with anchor/positive/negative RNA sequence triplets.

## Features

- **Modular Design**: Clean separation of concerns with dedicated modules for different functionalities
- **RNA Sequence Generation**: Random RNA sequence generation with configurable length distributions
- **Structure Folding**: Integration with ViennaRNA's RNAfold for secondary structure prediction
- **Structural Modifications**: Apply insertions/deletions to stems, hairpin loops, internal loops, bulges, and multi-loops
- **Negative Sample Generation**: Dinucleotide shuffling to preserve sequence composition
- **Appending Events**: Optional sequence appending with linkers
- **Parallel Processing**: Multi-threaded generation for improved performance
- **Dataset Splitting**: Automatic train/validation split functionality
- **Comprehensive Logging**: Detailed logging with configurable levels

## Project Structure

```
python_version/
├── main.py                    # Main entry point
├── requirements.txt           # Python dependencies
├── README.md                 # This file
├── config/                   # Configuration modules
│   ├── __init__.py
│   └── cli_parser.py         # Command-line argument parsing
├── core/                     # Core functionality modules
│   ├── __init__.py
│   ├── models.py             # Data models and classes
│   ├── rna_generator.py      # RNA sequence generation and folding
│   ├── modification_engine.py # RNA modification logic
│   ├── negative_generator.py  # Negative sample generation
│   └── dataset_generator.py  # Main dataset orchestration
└── utils/                    # Utility modules
    ├── __init__.py
    ├── logger.py             # Logging configuration
    └── output_handler.py     # Dataset saving and analysis
```

## Installation

1. **Prerequisites**: Install ViennaRNA package for RNA folding
   ```bash
   # On Ubuntu/Debian
   sudo apt-get install vienna-rna
   
   # On macOS with Homebrew
   brew install viennarna
   
   # Or conda
   conda install -c bioconda viennarna
   ```

2. **Python Dependencies**:
   ```bash
   cd python_version
   pip install -r requirements.txt
   ```

## Usage

### Basic Usage

Generate 100 RNA triplets with default settings:

```bash
python main.py --num_structures 100
```

### Advanced Usage

Generate a larger dataset with custom parameters:

```bash
python main.py \
    --num_structures 10000 \
    --seq_min_len 50 \
    --seq_max_len 200 \
    --seq_len_distribution norm \
    --seq_len_mean 100 \
    --seq_len_sd 20 \
    --n_stem_indels_min 1 \
    --n_stem_indels_max 3 \
    --n_hloop_indels_min 1 \
    --n_hloop_indels_max 2 \
    --mod_normalization \
    --normalization_len 100 \
    --num_workers 8 \
    --batch_size 64 \
    --output_dir output_python \
    --split \
    --train_fraction 0.8 \
    --val_fraction 0.2 \
    --debug
```

### Key Parameters

- `--num_structures`: Number of triplets to generate
- `--seq_min_len/--seq_max_len`: Sequence length range
- `--seq_len_distribution`: Length distribution (unif/norm)
- `--n_*_indels_min/max`: Modification count ranges for different structural elements
- `--mod_normalization`: Scale modifications by sequence length
- `--num_workers`: Number of parallel workers
- `--split`: Enable train/validation splitting
- `--debug`: Enable debug logging

## Module Overview

### Core Modules

1. **models.py**: Data structures and enums
   - `RnaTriplet`: Complete triplet data structure
   - `BulgeGraph`: RNA secondary structure representation
   - `NodeType`: Structural element classification
   - `ModificationCounts`: Tracking applied modifications

2. **rna_generator.py**: RNA sequence operations
   - `RnaGenerator`: Random sequence generation and folding
   - `BulgeGraphParser`: Structure parsing integration

3. **modification_engine.py**: Structural modifications
   - `ModificationEngine`: Applies insertions/deletions to structural elements
   - Handles stems, hairpin loops, internal loops, bulges, and multi-loops

4. **negative_generator.py**: Negative sample creation
   - `NegativeSampleGenerator`: Dinucleotide shuffling
   - `AppendingEngine`: Optional sequence appending

5. **dataset_generator.py**: Main orchestration
   - `DatasetGenerator`: Coordinates all components for triplet generation
   - Supports both sequential and parallel processing

### Utility Modules

1. **logger.py**: Logging configuration
2. **output_handler.py**: Dataset saving and analysis
   - CSV export with metadata
   - Train/validation splitting
   - Dataset statistics

## Output Format

The generator produces CSV files with the following columns:

- `triplet_id`: Unique identifier
- `anchor_seq/anchor_structure`: Original sequence and structure
- `positive_seq/positive_structure`: Modified sequence and structure
- `negative_seq/negative_structure`: Shuffled sequence and structure
- `*_modifications`: Counts of modifications by structural element type

Metadata is saved as JSON and included as a comment in CSV files.

## Comparison with Rust Version

This Python version offers several improvements over the original Rust implementation:

1. **Modularity**: Better separation of concerns with dedicated modules
2. **Readability**: More accessible codebase for researchers
3. **Extensibility**: Easier to add new features and modifications
4. **Integration**: Simpler integration with Python-based ML pipelines
5. **Debugging**: More straightforward debugging and profiling

## Performance Notes

- Use `--num_workers > 1` for parallel processing on multi-core systems
- Adjust `--batch_size` based on available memory
- Enable `--debug` for detailed execution logging
- The ViennaRNA dependency is required for RNA folding functionality

## Dependencies

- **ViennaRNA**: RNA secondary structure prediction
- **NumPy/SciPy**: Numerical computations and distributions
- **tqdm**: Progress bars
- **Python 3.8+**: Modern Python features

## License

This project follows the same license as the original Rust implementation.