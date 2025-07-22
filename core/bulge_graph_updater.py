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
        original = node.positions.copy() if node is not None else None

        if left_index < right_index:
            cls._shift_indices(bg, right_index, 1)
            cls._shift_indices(bg, left_index, 1)
        else:
            cls._shift_indices(bg, left_index, 1)
            cls._shift_indices(bg, right_index, 1)

        if node is not None and original and len(original) == 4:
            l_start, l_end, r_start, r_end = original
            left_pos = left_index + 1
            right_pos = right_index + 1

            l_start_new, l_end_new = node.positions[0], node.positions[1]
            r_start_new, r_end_new = node.positions[2], node.positions[3]

            if left_pos <= l_start:
                l_start_new = left_pos
            elif left_pos > l_end:
                l_end_new = left_pos

            if right_pos >= r_end:
                r_end_new = right_pos
            elif right_pos < r_start:
                r_start_new = right_pos

            node.positions = [l_start_new, l_end_new, r_start_new, r_end_new]
            cls._recompute_bounds(node)
        elif node is not None:
            node.positions.extend([left_index + 1, right_index + 1])
            node.positions.sort()
            cls._recompute_bounds(node)

    @classmethod
    def delete_stem_pair(
        cls, bg: BulgeGraph, node_name: str, left_index: int, right_index: int
    ) -> None:
        """Update graph for a stem pair deletion."""
        node = bg.elements.get(node_name)
        original = node.positions.copy() if node is not None else None

        if node is not None and original and len(original) == 4:
            l_start, l_end, r_start, r_end = original
            left_pos = left_index + 1
            right_pos = right_index + 1

            if left_pos == l_start:
                l_start += 1
            if left_pos == l_end:
                l_end -= 1
            if right_pos == r_start:
                r_start += 1
            if right_pos == r_end:
                r_end -= 1

            node.positions = [l_start, l_end, r_start, r_end]
            cls._recompute_bounds(node)
        elif node is not None:
            for pos in (left_index + 1, right_index + 1):
                if pos in node.positions:
                    node.positions.remove(pos)
            cls._recompute_bounds(node)

        if left_index < right_index:
            cls._shift_indices(bg, right_index, -1)
            cls._shift_indices(bg, left_index, -1)
        else:
            cls._shift_indices(bg, left_index, -1)
            cls._shift_indices(bg, right_index, -1)
