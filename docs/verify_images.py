#!/usr/bin/env python3
"""
Test script to verify all generated pattern images are valid.
"""

import sys
from pathlib import Path
from PIL import Image

def verify_images():
    """Verify all pattern images are valid PNG files."""
    docs_dir = Path(__file__).parent
    pattern_dir = docs_dir / "_static" / "pattern_images"
    
    if not pattern_dir.exists():
        print(f"Error: Pattern directory not found: {pattern_dir}")
        return False
    
    errors = []
    success_count = 0
    
    # Check all subdirectories
    for filter_type in ["pains", "reactive", "heterocycle"]:
        subdir = pattern_dir / filter_type
        if not subdir.exists():
            errors.append(f"Missing directory: {subdir}")
            continue
        
        # Check all PNG files in the directory
        for img_file in subdir.glob("*.png"):
            try:
                # Try to open and verify the image
                img = Image.open(img_file)
                img.verify()
                # Check dimensions
                img = Image.open(img_file)  # Re-open after verify
                width, height = img.size
                if width != 300 or height != 300:
                    errors.append(f"{img_file.name}: Wrong size {width}x{height} (expected 300x300)")
                else:
                    success_count += 1
            except Exception as e:
                errors.append(f"{img_file.name}: {str(e)}")
    
    # Print results
    print(f"\nImage Verification Results")
    print("=" * 60)
    print(f"Valid images: {success_count}")
    
    if errors:
        print(f"Errors found: {len(errors)}")
        for error in errors:
            print(f"  - {error}")
        return False
    else:
        print("All images are valid!")
        return True

if __name__ == "__main__":
    success = verify_images()
    sys.exit(0 if success else 1)
