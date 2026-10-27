"""
One possible haplotype for a locus with several joined breakpoint pairs.

The breakpoints of a locus, ordered by position, cut it into segments: 0 is
the left flank, k the segment between breakpoints k-1 and k, the last one the
right flank. Joined pairs (deletions and tandem duplications) are applied
outermost first. An event nested in a duplicated segment goes to one of its
copies, chosen at random, or to both copies for a deletion whose segment is
less covered than the flanks (it is then missing from the haplotype).
Events that overlap without nesting are not modelled.
"""

from __future__ import annotations

import random

BOTH_COPIES_DEPTH = 0.8


def _first_run(pieces: list[dict], segs: set[int], tags: tuple) -> tuple[int, int] | None:
    """First contiguous run of pieces from *segs* carrying all *tags*."""
    start = None
    for i, p in enumerate(pieces):
        inside = p["seg"] in segs and all(t in p["tags"] for t in tags)
        if inside and start is None:
            start = i
        elif not inside and start is not None:
            return start, i
    return (start, len(pieces)) if start is not None else None


def build_model(n_bkps: int, events: list[tuple], depths: list[float] | None = None,
                seed: int = 0) -> tuple[list[dict], list[str]]:
    """Apply *events* ``[(label, sv_type, i, j)]``, with *i* < *j* the indices
    of the two breakpoints in positional order, to the reference segments.

    Returns the pieces of the haplotype, ``{"seg", "copy", "junction"}`` with
    *copy* the copy of the outermost duplication the piece belongs to and
    *junction* the label of a new junction right before it, and notes on how
    the model was built.
    """
    rng = random.Random(seed)
    pieces = [{"seg": k, "tags": (), "junction": None} for k in range(n_bkps + 1)]
    notes: list[str] = []
    applied: list[set[int]] = []
    dups: dict[int, tuple[str, set[int]]] = {}

    for label, sv_type, i, j in sorted(events, key=lambda e: e[3] - e[2], reverse=True):
        segs = set(range(i + 1, j + 1))
        if sv_type not in ("deletion", "duplication"):
            notes.append(f"{label}: SV type unknown, not modelled")
            continue
        if any(segs & other and not (segs <= other or other <= segs) for other in applied):
            notes.append(f"{label} overlaps another junction without nesting: not modelled")
            continue
        applied.append(segs)

        targets: list[tuple] = [()]
        containers = [(d, s) for d, (_, s) in dups.items() if segs < s]
        if containers:
            dup_id = min(containers, key=lambda c: len(c[1]))[0]
            dup_label = dups[dup_id][0]
            if (sv_type == "deletion" and depths is not None
                    and max(depths[k] for k in segs) < BOTH_COPIES_DEPTH):
                targets = [((dup_id, 1),), ((dup_id, 2),)]
                notes.append(f"{label} deletion on both copies of {dup_label}: "
                             "its segment is less covered than the flanks")
            else:
                copy = rng.choice((1, 2))
                targets = [((dup_id, copy),)]
                notes.append(f"{label} {sv_type} placed on copy {copy} of {dup_label}, "
                             "chosen at random")

        for tags in targets:
            run = _first_run(pieces, segs, tags)
            if run is None:
                continue
            start, end = run
            if sv_type == "deletion":
                del pieces[start:end]
                if start < len(pieces):
                    pieces[start] = {**pieces[start], "junction": label}
            else:
                dup_id = len(dups)
                dups[dup_id] = (label, segs)
                copy1 = [{**p, "tags": p["tags"] + ((dup_id, 1),)} for p in pieces[start:end]]
                copy2 = [{**p, "tags": p["tags"][:-1] + ((dup_id, 2),), "junction": None}
                         for p in copy1]
                copy2[0]["junction"] = label
                pieces[start:end] = copy1 + copy2

    for p in pieces:
        p["copy"] = p["tags"][0][1] if p["tags"] else None
    return pieces, notes
