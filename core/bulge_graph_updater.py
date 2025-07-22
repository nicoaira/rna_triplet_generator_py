"""Utilities for updating bulge graph indices after sequence modifications."""

from .models import BulgeGraph, GraphNode


class BulgeGraphUpdater:
    """Update BulgeGraph positions without reparsing the structure."""

    @staticmethod
    def _shift_indices(bg: BulgeGraph, index: int, delta: int) -> None:
        """Shift all coordinates after a given index by delta.

        Args:
            bg: BulgeGraph to modify.
            index: 0-based position where the modification occurred.
            delta: Amount to shift indices (>0 for insert, <0 for delete).
        """
        threshold = index + 1
        for node in bg.elements.values():
            node.positions = [
                pos + delta if pos >= threshold else pos
                for pos in node.positions
            ]
            if node.start is not None and node.start >= threshold:
                node.start += delta
            if node.end is not None and node.end >= threshold:
                node.end += delta

    @staticmethod
    def _recompute_bounds(node: GraphNode) -> None:
        if node.positions:
            node.start = min(node.positions)
            node.end = max(node.positions)
        else:
            node.start = None
            node.end = None

    @classmethod
    def insert_loop_base(
        cls, bg: BulgeGraph, node_name: str, index: int
    ) -> None:
        """Update graph for a loop base insertion."""
        cls._shift_indices(bg, index, 1)
        node = bg.elements.get(node_name)
        if node is not None:
            node.positions.append(index + 1)
            node.positions.sort()
            cls._recompute_bounds(node)

    @classmethod
    def delete_loop_base(
        cls, bg: BulgeGraph, node_name: str, index: int
    ) -> None:
        """Update graph for a loop base deletion."""
        node = bg.elements.get(node_name)
        if node is not None and (index + 1) in node.positions:
            node.positions.remove(index + 1)
            cls._recompute_bounds(node)
        cls._shift_indices(bg, index, -1)

    @classmethod
    def insert_stem_pair(
        cls, bg: BulgeGraph, node_name: str, left_index: int, right_index: int
    ) -> None:
        """Update graph for a stem pair insertion."""
        node = bg.elements.get(node_name)
        orig = None
        if node is not None and len(node.positions) == 4:
            orig = list(node.positions)

        if left_index < right_index:
            cls._shift_indices(bg, right_index, 1)
            cls._shift_indices(bg, left_index, 1)
        else:
            cls._shift_indices(bg, left_index, 1)
            cls._shift_indices(bg, right_index, 1)

        if orig is not None:
            ls, le, rs, re = orig
            node.positions = [ls, le + 1, rs + 1, re + 2]
            node.start = node.positions[0]
            node.end = node.positions[-1]

    @classmethod
    def delete_stem_pair(
        cls, bg: BulgeGraph, node_name: str, left_index: int, right_index: int
    ) -> None:
        """Update graph for a stem pair deletion."""
        node = bg.elements.get(node_name)
        orig = None
        if node is not None and len(node.positions) == 4:
            orig = list(node.positions)

        if left_index < right_index:
            cls._shift_indices(bg, left_index, -1)
            cls._shift_indices(bg, right_index - 1, -1)
        else:
            cls._shift_indices(bg, right_index, -1)
            cls._shift_indices(bg, left_index - 1, -1)

        if node is None or orig is None:
            return

        ls, le, rs, re = orig
        is_single = ls == le and rs == re

        if is_single:
            # Remove the stem and renumber subsequent stems
            index = int(node_name[1:])
            del bg.elements[node_name]
            i = index + 1
            while f"s{i}" in bg.elements:
                bg.elements[f"s{i-1}"] = bg.elements.pop(f"s{i}")
                i += 1
        else:
            node.positions = [ls, le - 1, rs - 1, re - 2]
            node.start = node.positions[0]
            node.end = node.positions[-1]
