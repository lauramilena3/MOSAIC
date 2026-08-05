#!/usr/bin/env python3

"""Repair wrapped Pharokka CDS identifiers without shortening them.

The script:

1. Collapses wrapped /ID qualifier values by removing whitespace.
2. Removes existing /protein_id qualifiers.
3. Adds exactly one /protein_id equal to each repaired /ID.
4. Atomically replaces the original pharokka.gbk.

Usage:
    python scripts/genbank_cleanup.py path/to/pharokka_output
    python scripts/genbank_cleanup.py path/to/pharokka.gbk
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
from pathlib import Path


QUALIFIER_START = re.compile(
    r'^(?P<indent>\s*)/(?P<name>[A-Za-z0-9_]+)="(?P<value>.*)$'
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Repair wrapped Pharokka CDS identifiers."
    )
    parser.add_argument(
        "input_path",
        type=Path,
        help="Pharokka output directory or direct path to pharokka.gbk",
    )
    return parser.parse_args()


def resolve_genbank_path(input_path: Path) -> Path:
    input_path = input_path.resolve()

    if input_path.is_file():
        return input_path

    if not input_path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    expected = input_path / "pharokka.gbk"

    if expected.is_file():
        return expected

    candidates = sorted(input_path.glob("*.gbk"))

    if not candidates:
        raise FileNotFoundError(
            f"No pharokka.gbk or other *.gbk file found in {input_path}"
        )

    if len(candidates) > 1:
        names = ", ".join(path.name for path in candidates)
        raise RuntimeError(
            f"Multiple GenBank files found: {names}. "
            "Pass the intended file directly."
        )

    return candidates[0]


def clean_identifier(value: str) -> str:
    """Remove whitespace introduced by GenBank continuation lines."""
    cleaned = re.sub(r"\s+", "", value)

    if not cleaned:
        raise ValueError("Encountered an empty CDS identifier")

    if '"' in cleaned:
        raise ValueError(
            f"Unexpected quote inside CDS identifier: {cleaned!r}"
        )

    return cleaned


def read_quoted_qualifier(
    lines: list[str],
    start_index: int,
) -> tuple[str, str, str, int]:
    """Read a possibly multiline quoted GenBank qualifier.

    Returns:
        indentation, qualifier name, complete value, next line index
    """
    first_line = lines[start_index].rstrip("\r\n")
    match = QUALIFIER_START.match(first_line)

    if match is None:
        raise ValueError(
            f"Invalid qualifier at line {start_index + 1}: {first_line!r}"
        )

    indent = match.group("indent")
    name = match.group("name")
    value_parts = [match.group("value")]
    index = start_index + 1

    while True:
        current = value_parts[-1]

        if current.endswith('"'):
            value_parts[-1] = current[:-1]
            break

        if index >= len(lines):
            raise ValueError(
                f"Unterminated /{name} qualifier starting at "
                f"line {start_index + 1}"
            )

        continuation = lines[index].strip()
        value_parts.append(continuation)
        index += 1

    return indent, name, "".join(value_parts), index


def repair_genbank_text(text: str) -> tuple[str, int]:
    lines = text.splitlines(keepends=True)
    output: list[str] = []

    repaired_ids: set[str] = set()
    repaired_count = 0
    index = 0

    while index < len(lines):
        stripped = lines[index].lstrip()

        if stripped.startswith('/protein_id="'):
            # Remove the existing protein_id, including continuation lines.
            _, _, _, index = read_quoted_qualifier(lines, index)
            continue

        if stripped.startswith('/ID="'):
            indent, _, raw_id, index = read_quoted_qualifier(lines, index)
            cleaned_id = clean_identifier(raw_id)

            if cleaned_id in repaired_ids:
                raise ValueError(
                    f"Duplicate CDS ID after cleanup: {cleaned_id!r}"
                )

            repaired_ids.add(cleaned_id)
            repaired_count += 1

            # Keep the complete identifier on one physical line.
            output.append(f'{indent}/ID="{cleaned_id}"\n')
            output.append(f'{indent}/protein_id="{cleaned_id}"\n')
            continue

        output.append(lines[index])
        index += 1

    return "".join(output), repaired_count


def validate_output(text: str, expected_count: int) -> None:
    id_values = re.findall(
        r'^[ \t]*/ID="([^"\r\n]+)"',
        text,
        flags=re.MULTILINE,
    )
    protein_values = re.findall(
        r'^[ \t]*/protein_id="([^"\r\n]+)"',
        text,
        flags=re.MULTILINE,
    )

    if len(id_values) != expected_count:
        raise RuntimeError(
            f"Expected {expected_count} ID qualifiers, "
            f"but found {len(id_values)}"
        )

    if len(protein_values) != expected_count:
        raise RuntimeError(
            f"Expected {expected_count} protein_id qualifiers, "
            f"but found {len(protein_values)}"
        )

    if id_values != protein_values:
        raise ValueError(
            "ID and protein_id qualifier values do not match"
        )

    for identifier in id_values:
        if any(character.isspace() for character in identifier):
            raise ValueError(
                f"Whitespace remains in identifier: {identifier!r}"
            )

    if len(set(id_values)) != len(id_values):
        raise ValueError("Duplicate CDS identifiers remain after cleanup")


def clean_genbank(gbk_path: Path) -> int:
    original_text = gbk_path.read_text(encoding="utf-8")

    cleaned_text, repaired_count = repair_genbank_text(original_text)

    if repaired_count == 0:
        raise ValueError(f"No /ID qualifiers found in {gbk_path}")

    validate_output(cleaned_text, repaired_count)

    temporary_path: Path | None = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{gbk_path.name}.",
            suffix=".tmp",
            dir=gbk_path.parent,
            delete=False,
        ) as temporary_file:
            temporary_path = Path(temporary_file.name)
            temporary_file.write(cleaned_text)
            temporary_file.flush()
            os.fsync(temporary_file.fileno())

        os.replace(temporary_path, gbk_path)
        temporary_path = None

    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()

    return repaired_count


def main() -> int:
    args = parse_arguments()

    try:
        gbk_path = resolve_genbank_path(args.input_path)
        repaired_count = clean_genbank(gbk_path)
    except Exception as error:
        print(f"GenBank cleanup failed: {error}", file=sys.stderr)
        return 1

    print(f"Cleaned GenBank file: {gbk_path}")
    print(f"CDS identifiers repaired: {repaired_count}")
    print(f"protein_id qualifiers written: {repaired_count}")
    print("Full Pharokka identifiers were preserved.")
    print("Original GenBank file replaced successfully.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
