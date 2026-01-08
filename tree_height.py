#!/usr/bin/env python3
"""
Calculate tree height (maximum root-to-tip distance) from a Newick tree.
Optionally rescale branch lengths and save to a new file.

Usage:
    python3 tree_height.py tree.newick
    python3 tree_height.py tree.newick --scale 0.0000000004
    python3 tree_height.py tree.newick --scale 0.0000000004 --output rescaled.newick
    cat tree.newick | python3 tree_height.py --scale 2.0
"""
import re
import sys
import argparse

def parse_newick_detailed(newick, scale=1.0, make_ultrametric=False):
    """Parse Newick and find the deepest leaf with its path"""

    class Node:
        def __init__(self, name="", length=0):
            self.name = name
            self.length = length
            self.children = []

        def to_newick(self):
            """Convert node back to Newick format"""
            if not self.children:
                # Leaf node
                if self.length > 0:
                    return f"{self.name}:{self.length}"
                else:
                    return self.name
            else:
                # Internal node
                children_str = ",".join([child.to_newick() for child in self.children])
                if self.name:
                    if self.length > 0:
                        return f"({children_str}){self.name}:{self.length}"
                    else:
                        return f"({children_str}){self.name}"
                else:
                    if self.length > 0:
                        return f"({children_str}):{self.length}"
                    else:
                        return f"({children_str})"

    def tokenize(s):
        """Tokenize Newick string"""
        tokens = []
        current = ""
        for char in s:
            if char in "(),;:":
                if current:
                    tokens.append(current)
                    current = ""
                if char != ";":
                    tokens.append(char)
            else:
                current += char
        if current:
            tokens.append(current)
        return tokens

    def build_tree(tokens, idx=0):
        """Build tree from tokens"""
        node = Node()

        if tokens[idx] == "(":
            idx += 1
            # Parse children
            while True:
                child, idx = build_tree(tokens, idx)
                node.children.append(child)
                if tokens[idx] == ",":
                    idx += 1
                elif tokens[idx] == ")":
                    idx += 1
                    break

        # Parse label
        if idx < len(tokens) and tokens[idx] not in "(),;:":
            node.name = tokens[idx]
            idx += 1

        # Parse branch length
        if idx < len(tokens) and tokens[idx] == ":":
            idx += 1
            if idx < len(tokens):
                node.length = float(tokens[idx]) * scale
                idx += 1

        return node, idx

    def find_deepest_leaf(node, depth=0, path=[]):
        """Find deepest leaf and its path"""
        current_depth = depth + node.length
        current_path = path + [(node.name if node.name else "internal", node.length)]

        if not node.children:  # Leaf
            return current_depth, current_path

        max_depth = 0
        max_path = []
        for child in node.children:
            child_depth, child_path = find_deepest_leaf(child, current_depth, current_path)
            if child_depth > max_depth:
                max_depth = child_depth
                max_path = child_path

        return max_depth, max_path

    def make_tree_ultrametric(node, target_depth, current_depth=0):
        """Adjust terminal branch lengths to make tree ultrametric"""
        current_depth += node.length

        if not node.children:  # Leaf node
            # Adjust terminal branch to reach target depth
            adjustment = target_depth - current_depth
            node.length += adjustment
        else:
            # Recursively process children
            for child in node.children:
                make_tree_ultrametric(child, target_depth, current_depth)

    newick = newick.strip().rstrip(';')
    tokens = tokenize(newick)
    tree, _ = build_tree(tokens)

    max_depth, path = find_deepest_leaf(tree)

    # Make tree ultrametric if requested
    if make_ultrametric:
        make_tree_ultrametric(tree, max_depth)
        # Recalculate after making ultrametric
        max_depth, path = find_deepest_leaf(tree)

    return max_depth, path, tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate tree height and optionally rescale branch lengths',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate tree height
  python3 tree_height.py tree.newick

  # Rescale branches (e.g., convert years to substitutions per site)
  python3 tree_height.py tree.newick --scale 0.0000000004

  # Rescale and save to new file
  python3 tree_height.py tree.newick --scale 0.0000000004 --output rescaled.newick

  # Make tree ultrametric (fixes rounding errors)
  python3 tree_height.py tree.newick --scale 180 --ultrametric --output ultrametric.newick

  # From stdin
  cat tree.newick | python3 tree_height.py --scale 2.0
        """
    )

    parser.add_argument('input', nargs='?', help='Input Newick file (or use stdin)')
    parser.add_argument('-s', '--scale', type=float, default=1.0,
                        help='Scaling factor for branch lengths (default: 1.0)')
    parser.add_argument('-u', '--ultrametric', action='store_true',
                        help='Force tree to be ultrametric by adjusting terminal branches')
    parser.add_argument('-o', '--output', type=str,
                        help='Output file for rescaled tree (Newick format)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress tree height output (useful with --output)')

    args = parser.parse_args()

    # Read tree from file or stdin
    if args.input:
        with open(args.input, 'r') as f:
            tree_str = f.read()
    else:
        tree_str = sys.stdin.read()

    max_depth, path, tree = parse_newick_detailed(tree_str, scale=args.scale, make_ultrametric=args.ultrametric)

    if not args.quiet:
        if args.scale != 1.0:
            print(f"Scaling factor: {args.scale}")
        if args.ultrametric:
            print(f"Ultrametric adjustment: ENABLED (terminal branches adjusted)")
        if args.scale != 1.0 or args.ultrametric:
            print()

        # Choose formatting based on tree height magnitude
        if max_depth < 0.01:
            # Use scientific notation for small values
            fmt = "{:.8e}"
            print(f"Tree height: {fmt.format(max_depth)}")
            print(f"\nDeepest leaf: {path[-1][0]}")
            print(f"\nPath from root to deepest leaf:")
            print(f"{'Node':<20} {'Branch Length':<20} {'Cumulative Distance':<20}")
            print("-" * 60)

            cumulative = 0
            for i, (name, length) in enumerate(path):
                cumulative += length
                if i == 0:
                    print(f"{'Root':<20} {fmt.format(length):<20} {fmt.format(cumulative):<20}")
                else:
                    node_label = name if name else f"internal_{i}"
                    print(f"{node_label:<20} {fmt.format(length):<20} {fmt.format(cumulative):<20}")
        else:
            # Use comma-separated format for large values
            print(f"Tree height: {max_depth:,.2f}")
            print(f"\nDeepest leaf: {path[-1][0]}")
            print(f"\nPath from root to deepest leaf:")
            print(f"{'Node':<20} {'Branch Length':<20} {'Cumulative Distance':<20}")
            print("-" * 60)

            cumulative = 0
            for i, (name, length) in enumerate(path):
                cumulative += length
                if i == 0:
                    print(f"{'Root':<20} {length:<20,.2f} {cumulative:<20,.2f}")
                else:
                    node_label = name if name else f"internal_{i}"
                    print(f"{node_label:<20} {length:<20,.2f} {cumulative:<20,.2f}")

    # Save rescaled tree if output file specified
    if args.output:
        rescaled_newick = tree.to_newick() + ";"
        with open(args.output, 'w') as f:
            f.write(rescaled_newick + "\n")
        if not args.quiet:
            print(f"\nRescaled tree saved to: {args.output}")
