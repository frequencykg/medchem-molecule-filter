#!/usr/bin/env python3
"""
Generate images for all SMARTS patterns used in the filters.

This script creates visual representations of all SMARTS patterns
for inclusion in the documentation.
"""

import os
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw

# Add parent directory to path to import medchem_filter
sys.path.insert(0, str(Path(__file__).parent.parent))

from medchem_filter.data.pains_patterns import PAINS_PATTERNS
from medchem_filter.data.reactive_patterns import REACTIVE_PATTERNS
from medchem_filter.data.heterocycle_patterns import HETEROCYCLE_PATTERNS


def generate_pattern_image(name, smarts, output_dir, img_size=(300, 300)):
    """
    Generate an image for a SMARTS pattern.
    
    Args:
        name: Pattern name
        smarts: SMARTS string
        output_dir: Directory to save the image
        img_size: Size of the output image (width, height)
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Create molecule from SMARTS
        mol = Chem.MolFromSmarts(smarts)
        if mol is None:
            print(f"Warning: Could not parse SMARTS for {name}: {smarts}")
            return False
        
        # Generate 2D coordinates
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)
        
        # Draw molecule
        img = Draw.MolToImage(mol, size=img_size)
        
        # Save image
        output_path = output_dir / f"{name}.png"
        img.save(output_path)
        print(f"Generated image: {output_path}")
        return True
        
    except Exception as e:
        print(f"Error generating image for {name}: {e}")
        return False


def generate_all_images(base_output_dir="_static/pattern_images"):
    """
    Generate images for all filter patterns.
    
    Args:
        base_output_dir: Base directory for output images
    """
    # Create output directories
    docs_dir = Path(__file__).parent
    output_dir = docs_dir / base_output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    pains_dir = output_dir / "pains"
    reactive_dir = output_dir / "reactive"
    heterocycle_dir = output_dir / "heterocycle"
    
    for d in [pains_dir, reactive_dir, heterocycle_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Generate PAINS pattern images
    print("\nGenerating PAINS pattern images...")
    print("=" * 60)
    pains_success = 0
    for name, smarts in PAINS_PATTERNS.items():
        if generate_pattern_image(name, smarts, pains_dir):
            pains_success += 1
    print(f"PAINS: Generated {pains_success}/{len(PAINS_PATTERNS)} images")
    
    # Generate Reactive pattern images
    print("\nGenerating Reactive pattern images...")
    print("=" * 60)
    reactive_success = 0
    for name, smarts in REACTIVE_PATTERNS.items():
        if generate_pattern_image(name, smarts, reactive_dir):
            reactive_success += 1
    print(f"Reactive: Generated {reactive_success}/{len(REACTIVE_PATTERNS)} images")
    
    # Generate Heterocycle pattern images
    print("\nGenerating Heterocycle pattern images...")
    print("=" * 60)
    heterocycle_success = 0
    for name, smarts in HETEROCYCLE_PATTERNS.items():
        if generate_pattern_image(name, smarts, heterocycle_dir):
            heterocycle_success += 1
    print(f"Heterocycle: Generated {heterocycle_success}/{len(HETEROCYCLE_PATTERNS)} images")
    
    print("\n" + "=" * 60)
    print(f"Total: Generated {pains_success + reactive_success + heterocycle_success} images")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    generate_all_images()
