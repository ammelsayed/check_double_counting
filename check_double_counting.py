#!/usr/bin/env python3

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

PDG = {
    "d": 1, "d~": -1, "u": 2, "u~": -2, "s": 3, "s~": -3,
    "c": 4, "c~": -4, "b": 5, "b~": -5, "t": 6, "t~": -6,
    "e-": 11, "e+": -11, "ve": 12, "ve~": -12,
    "mu-": 13, "mu+": -13, "vm": 14, "vm~": -14,
    "ta-": 15, "ta+": -15, "vt": 16, "vt~": -16,
    "g": 21, "a": 22, "z": 23, "w+": 24, "w-": -24, "h": 25,
}

PROCESS_RE = re.compile(r"^Process:\s*(.*?)\s+QCD<=99")
DECAY_RE = re.compile(r"^\s*Decay:\s*(.*?)\s+WEIGHTED<=")
DIAGRAM_RE = re.compile(r"^\s*(\d+)\s+(\(\(.*\)\))\s+\(QCD=")
TOKEN_RE = re.compile(r"\d+\((-?\d+)\)")


def split_top_level(s):
    out, start, depth = [], 0, 0
    for i, c in enumerate(s):
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == "," and depth == 0:
            if s[start:i].strip():
                out.append(s[start:i].strip())
            start = i + 1
    if s[start:].strip():
        out.append(s[start:].strip())
    return out


def parse_vertices(topology):
    vertices = []
    for vertex in split_top_level(topology[1:-1]):
        body = vertex[1:-1]
        connection, interaction_id = body.rsplit(",id:", 1)
        if ">" in connection:
            left, right = connection.split(">", 1)
        else:
            left, right = connection, ""
        left = [int(x) for x in TOKEN_RE.findall(left)]
        right = [int(x) for x in TOKEN_RE.findall(right)]
        vertices.append((left, right, int(interaction_id)))
    return vertices


def parse_file(path):
    processes = []
    current = None
    for line in Path(path).read_text().splitlines():
        match = PROCESS_RE.match(line)
        if match:
            current = {
                "process": match.group(1).strip(),
                "decays": [],
                "diagrams": [],
            }
            processes.append(current)
            continue
        if current is None:
            continue

        match = DECAY_RE.match(line)
        if match:
            parent, daughters = match.group(1).split(">", 1)
            current["decays"].append(
                (parent.strip(), daughters.strip().split())
            )
            continue

        match = DIAGRAM_RE.match(line)
        if match:
            current["diagrams"].append(
                (int(match.group(1)), match.group(2))
            )
    return processes


def process_parts(process):
    left, right = process["process"].split(">", 1)
    initial = tuple(sorted(left.split()))
    final = tuple(right.split())
    return initial, final


def add_cp_modes(modes):
    result = defaultdict(list)
    for parent, daughters in modes.items():
        for decay in daughters:
            if decay not in result[parent]:
                result[parent].append(decay)
            cp_parent = -parent
            cp_decay = tuple(sorted(-x for x in decay))
            if cp_decay not in result[cp_parent]:
                result[cp_parent].append(cp_decay)
    return result


def make_modes(process):
    modes = defaultdict(list)
    for parent, daughters in process["decays"]:
        if parent not in PDG or any(x not in PDG for x in daughters):
            raise ValueError(
                f"Unknown particle name in decay: {parent} -> {' '.join(daughters)}"
            )
        modes[PDG[parent]].append(tuple(PDG[x] for x in daughters))
    return add_cp_modes(modes)


def leaf_options(modes, particle, stack=()):
    if particle in stack:
        return {(particle,)}
    if particle not in modes:
        return {(particle,)}

    result = set()
    for daughters in modes[particle]:
        options = {()}
        for daughter in daughters:
            child_options = leaf_options(modes, daughter, stack + (particle,))
            options = {
                tuple(sorted(a + b))
                for a in options
                for b in child_options
            }
        result.update(options)
    return result


def final_state_options(process):
    modes = make_modes(process)
    options = {()}
    for name in process_parts(process)[1]:
        if name not in PDG:
            raise ValueError(f"Unknown particle name: {name}")
        options = {
            tuple(sorted(a + b))
            for a in options
            for b in leaf_options(modes, PDG[name])
        }
    return options


def decay_signatures(process, combined_modes):
    signatures = set()
    for parent_name, daughter_names in process["decays"]:
        parent = PDG[parent_name]
        options = {()}
        for name in daughter_names:
            options = {
                tuple(sorted(a + b))
                for a in options
                for b in leaf_options(combined_modes, PDG[name])
            }
        for leaves in options:
            signatures.add((parent, leaves))
    return signatures


def realizes(diagram, signature, combined_modes):
    parent, target_leaves = signature

    for left, right, _ in parse_vertices(diagram[1]):
        particles = left + right
        if parent not in particles:
            continue

        particles = particles.copy()
        particles.remove(parent)
        options = {()}
        for particle in particles:
            options = {
                tuple(sorted(a + b))
                for a in options
                for b in leaf_options(combined_modes, particle)
            }

        if target_leaves in options:
            return True

    return False


def same_initial_state(process1, process2):
    return process_parts(process1)[0] == process_parts(process2)[0]


def find_overlaps(file1, file2):
    processes1 = parse_file(file1)
    processes2 = parse_file(file2)
    results = []

    for p1 in processes1:
        for p2 in processes2:
            if not same_initial_state(p1, p2):
                continue

            final1 = final_state_options(p1)
            final2 = final_state_options(p2)
            common_final_states = final1 & final2
            if not common_final_states:
                continue

            modes1 = make_modes(p1)
            modes2 = make_modes(p2)
            combined_modes = defaultdict(list)
            for modes in (modes1, modes2):
                for parent, decays in modes.items():
                    for decay in decays:
                        if decay not in combined_modes[parent]:
                            combined_modes[parent].append(decay)

            required = (
                decay_signatures(p1, combined_modes)
                | decay_signatures(p2, combined_modes)
            )

            matches1 = [
                n for n, d in p1["diagrams"]
                if all(realizes((n, d), sig, combined_modes) for sig in required)
            ]
            matches2 = [
                n for n, d in p2["diagrams"]
                if all(realizes((n, d), sig, combined_modes) for sig in required)
            ]

            if matches1 and matches2:
                results.append({
                    "process1": p1["process"],
                    "process2": p2["process"],
                    "final_state": next(iter(common_final_states)),
                    "file1_diagrams": matches1,
                    "file2_diagrams": matches2,
                })

    return results


def format_state(state):
    names = {v: k for k, v in PDG.items()}
    return " ".join(names.get(x, str(x)) for x in state)


def main():
    parser = argparse.ArgumentParser(
        description="Check two MadGraph diagram text files for physical double counting."
    )
    parser.add_argument("file1")
    parser.add_argument("file2")
    args = parser.parse_args()

    results = find_overlaps(args.file1, args.file2)

    total1 = sum(len(x["file1_diagrams"]) for x in results)
    total2 = sum(len(x["file2_diagrams"]) for x in results)

    print(f"File 1: {args.file1}")
    print(f"File 2: {args.file2}")
    print()
    print(f"Overlapping process pairs: {len(results)}")
    print(f"Duplicated diagrams in file 1: {total1}")
    print(f"Duplicated diagrams in file 2: {total2}")
    print()

    if not results:
        print("No duplicated diagrams found.")
        return

    print("Duplicated diagrams:")
    for result in results:
        print()
        print(f"  {result['process1']}  <->  {result['process2']}")
        print(f"  Final state: {format_state(result['final_state'])}")
        print(f"  File 1 diagrams: {', '.join(map(str, result['file1_diagrams']))}")
        print(f"  File 2 diagrams: {', '.join(map(str, result['file2_diagrams']))}")


if __name__ == "__main__":
    main()
