"""snapgene_edit - Read, edit, and write SnapGene .dna files.

This script demonstrates the complete workflow.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from snapgene_edit import SnapGeneEditor

# ========== Example 1: Load and inspect a .dna file ==========
demo_file = r"d:\1lab_document\plasmid\za_sequence.dna"
if os.path.exists(demo_file):
    print("=== Example 1: Read a .dna file ===")
    editor = SnapGeneEditor(demo_file)
    print(f"  File:     {os.path.basename(demo_file)}")
    print(f"  Sequence: {editor.sequence[:50]}... ({editor.length} bp)")
    print(f"  Circular: {editor.is_circular}")
    print(f"  Features: {len(editor.features)}")

    for f in editor.list_features():
        print(f"    - {f['name']} ({f['type']}) [{f['start']}-{f['end']}] "
              f"strand={f['strand']}")
    print()

# ========== Example 2: Edit a sequence and features ==========
print("=== Example 2: Edit sequence and features ===")
editor2 = SnapGeneEditor()
editor2.name = "my_plasmid"
editor2.set_sequence("ATGCGATCGTAGCTAGCTAGCTAGCGCGATCGATCGTACGTAGC")
editor2.is_circular = True

# Add features / annotations
editor2.add_feature("GFP", 0, 15, feature_type="CDS", strand=1, color="#00FF00")
editor2.add_feature("Promoter", 20, 35, feature_type="promoter", strand=1, color="#FF0000")
editor2.add_feature("Terminator", 40, 50, feature_type="terminator", strand=-1, color="#0000FF")

print(f"  Created plasmid: {editor2.name}, {editor2.length} bp")
print(f"  Features: {len(editor2.features)}")

# Edit features
editor2.rename_feature("GFP", "mGFP")
editor2.update_feature_coords("mGFP", 5, 20)
print(f"  After rename: {editor2.list_features()[0]['name']} ({editor2.list_features()[0]['start']}-{editor2.list_features()[0]['end']})")

# Edit sequence
editor2.insert(15, "AAAA")
editor2.replace_region(10, 14, "GGGG")
print(f"  After edits: {editor2.sequence}")

# Save the modified file
output_path = os.path.join(os.path.dirname(demo_file) if os.path.exists(demo_file) else os.getcwd(),
                           "output_edited.dna")
editor2.save(output_path)
print(f"\n  Saved to: {output_path}")
print()

# ========== Example 3: Verify round-trip preserves data ==========
if os.path.exists(demo_file):
    print("=== Example 3: Round-trip (read -> write -> read) ===")
    import tempfile
    tmp = os.path.join(tempfile.gettempdir(), "rt_test.dna")
    editor3 = SnapGeneEditor(demo_file)
    editor3.save(tmp)
    editor4 = SnapGeneEditor(tmp)
    seq_match = editor3.sequence == editor4.sequence
    feat_match = len(editor3.features) == len(editor4.features)
    circ_match = editor3.is_circular == editor4.is_circular
    print(f"  Sequence preserved: {seq_match}")
    print(f"  Features preserved: {feat_match}")
    print(f"  Circular flag preserved: {circ_match}")
    os.remove(tmp)

print("\nDone!")
