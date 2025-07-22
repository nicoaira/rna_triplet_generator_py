import json
from pathlib import Path
from difflib import SequenceMatcher
import sys

import pytest

# Ensure the repository root is on the Python path so that `core` can be imported
ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))

from core.models import BulgeGraph, GraphNode
from core.bulge_graph_updater import BulgeGraphUpdater


def _build_graph(mapping):
    elements = {}
    for name, coords in mapping.items():
        elements[name] = GraphNode(
            positions=coords.copy(),
            start=min(coords) if coords else None,
            end=max(coords) if coords else None,
        )
    return BulgeGraph(elements=elements)


def _diff_indices(pre: str, post: str):
    sm = SequenceMatcher(None, pre, post)
    inserted, deleted = [], []
    for tag, a1, a2, b1, b2 in sm.get_opcodes():
        if tag == "insert":
            inserted.extend(range(a1, a1 + (b2 - b1)))
        elif tag == "delete":
            deleted.extend(range(a1, a2))
    return inserted, deleted


def load_cases():
    path = Path(__file__).parent / "test_cases.json"
    with open(path) as f:
        data = json.load(f)
    return data["test_cases"]


def _format_modification_info(case, inserted, deleted):
    """Format detailed information about the RNA modification."""
    info_lines = []
    
    # Basic modification info
    info_lines.append(f"Modification: {case['mod_action']} on {case['mod_node_type']} '{case['mod_node']}'")
    info_lines.append(f"Original coordinates: {case.get('original_coords', 'N/A')}")
    
    # Sequence change info
    if inserted:
        info_lines.append(f"Inserted at positions: {inserted}")
    if deleted:
        info_lines.append(f"Deleted from positions: {deleted}")
    
    info_lines.append(f"Sequence length change: {case.get('sequence_length_change', 'N/A')}")
    
    # Structure info
    info_lines.append(f"Pre-modification sequence:  {case.get('pre_mod_seq', 'N/A')}")
    info_lines.append(f"Pre-modification structure: {case.get('pre_mod_dotbracket', 'N/A')}")
    info_lines.append(f"Post-modification sequence: {case.get('post_mod_seq', 'N/A')}")
    info_lines.append(f"Post-modification structure:{case.get('post_mod_dotbracket', 'N/A')}")
    
    return '\n'.join(info_lines)


def _format_graph_differences(result, expected):
    """Format the differences between actual and expected graph results."""
    diff_lines = []
    
    all_nodes = set(result.keys()) | set(expected.keys())
    
    for node in sorted(all_nodes):
        result_pos = result.get(node, [])
        expected_pos = expected.get(node, [])
        
        if result_pos != expected_pos:
            diff_lines.append(f"  {node}: got {result_pos}, expected {expected_pos}")
        elif len(diff_lines) < 10:  # Show some matching ones too, but limit output
            diff_lines.append(f"  {node}: {result_pos} ✓")
    
    return '\n'.join(diff_lines)


@pytest.mark.parametrize("case", load_cases())
def test_bulge_graph_update(case):
    bg = _build_graph(case["pre_mod_bulgegraph"])
    inserted, deleted = _diff_indices(case["pre_mod_seq"], case["post_mod_seq"])

    node = case["mod_node"]
    action = case["mod_action"]

    if action == "insert_pair":
        left, right = sorted(inserted)
        BulgeGraphUpdater.insert_stem_pair(bg, node, left, right)
    elif action == "delete_pair":
        left, right = sorted(deleted)
        BulgeGraphUpdater.delete_stem_pair(bg, node, left, right)
    elif action == "insert":
        BulgeGraphUpdater.insert_loop_base(bg, node, inserted[0])
    elif action == "delete":
        BulgeGraphUpdater.delete_loop_base(bg, node, deleted[0])

    result = {n: g.positions for n, g in bg.elements.items()}
    
    # Create detailed error message if assertion fails
    if result != case["post_mod_bulgegraph"]:
        modification_info = _format_modification_info(case, inserted, deleted)
        graph_differences = _format_graph_differences(result, case["post_mod_bulgegraph"])
        
        error_msg = f"""
RNA Structure Modification Test Failed
=====================================

{modification_info}

Graph Node Differences:
{graph_differences}

Case ID: {case.get('case_id', 'N/A')}
Attempt: {case.get('attempt', 'N/A')}
"""
        pytest.fail(error_msg)
    
    assert result == case["post_mod_bulgegraph"]