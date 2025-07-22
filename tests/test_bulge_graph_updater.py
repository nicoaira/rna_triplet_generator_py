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
    assert result == case["post_mod_bulgegraph"]