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
        """Update graph for a loop base deletion.

        ``BulgeGraph`` nodes often store only the start and end coordinates of
        a loop (e.g. hairpins, bulges or multiloops).  When a base inside such a
        loop is removed we need to shrink the loop while keeping the start/end
        representation consistent.  This helper therefore expands two/"four"
        coordinate representations to explicit ranges, removes the deleted base,
        applies the global coordinate shift and then converts the positions back
        to the condensed form used by ``forgi``.
        """

        node = bg.elements.get(node_name)
        pos = index + 1

        if node is not None:
            original = list(node.positions)

            # Expand start/end style representations so we can operate on
            # explicit coordinates.
            if len(original) == 2:
                expanded = list(range(original[0], original[1] + 1))
            elif len(original) == 4:
                ls, le, rs, re = original
                expanded = list(range(ls, le + 1)) + list(range(rs, re + 1))
                left_size_orig = le - ls + 1
                right_size_orig = re - rs + 1
            else:
                expanded = original.copy()

            if pos in expanded:
                expanded.remove(pos)

            # Temporarily assign the expanded list so that the subsequent shift
            # operation updates this node along with all others.
            node.positions = expanded

        # Shift coordinates for all nodes, including the modified one.
        cls._shift_indices(bg, index, -1)

        if node is not None:
            updated = sorted(node.positions)

            if len(original) == 2:
                node.positions = [updated[0], updated[-1]] if updated else []
            elif len(original) == 4:
                # Determine the length of each side from the original coords
                ls, le, rs, re = original
                left_len = left_size_orig
                right_len = right_size_orig

                # Adjust expected lengths depending on which side lost a base
                if ls <= pos <= le:
                    left_len -= 1
                elif rs <= pos <= re:
                    right_len -= 1

                left = updated[:left_len]
                right = updated[left_len:left_len + right_len]

                if left and right:
                    node.positions = [left[0], left[-1], right[0], right[-1]]
                elif left:
                    node.positions = [left[0], left[-1]]
                elif right:
                    node.positions = [right[0], right[-1]]
                else:
                    node.positions = []
            else:
                node.positions = updated

            cls._recompute_bounds(node)

    @classmethod
    def insert_stem_pair(
        cls, bg: BulgeGraph, node_name: str, left_index: int, right_index: int
    ) -> None:
        """Update graph for a stem pair insertion."""
        node = bg.elements.get(node_name)
        orig = None
        if node is not None and len(node.positions) == 4:
            orig = list(node.positions)
            ls, le, rs, re = orig
            # Clamp diff-based indices to expected ranges
            left_index = min(max(left_index, ls - 1), le)
            right_index = min(max(right_index, rs - 1), re)

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
            ls, le, rs, re = orig
            left_index = min(max(left_index, ls - 1), le - 1)
            right_index = min(max(right_index, rs - 1), re - 1)

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

            # Merge adjacent internal loops when possible
            left_name = f"i{index-2}"
            right_name = f"i{index-1}"
            left_loop = bg.elements.get(left_name)
            right_loop = bg.elements.get(right_name)
            if left_loop and right_loop:
                merged = sorted(set(left_loop.positions + right_loop.positions))
                if len(merged) > 2:
                    merged = [merged[0], merged[-1]]
                left_loop.positions = merged
                cls._recompute_bounds(left_loop)
                del bg.elements[right_name]
                j = index
                while f"i{j}" in bg.elements:
                    bg.elements[f"i{j-1}"] = bg.elements.pop(f"i{j}")
                    j += 1
        else:
            node.positions = [ls, le - 1, rs - 1, re - 2]
            node.start = node.positions[0]
            node.end = node.positions[-1]
