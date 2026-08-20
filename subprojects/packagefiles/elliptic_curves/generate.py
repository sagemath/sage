"""Generate the databases shipped by the ``elliptic_curves`` package."""

import argparse
import sqlite3
from contextlib import ExitStack
from pathlib import Path


def create_cremona_database(allcurves: Path, output: Path) -> None:
    """Create the mini Cremona SQLite database from ``allcurves``."""
    output.unlink(missing_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    class_data = []
    curve_data = []
    with allcurves.open(encoding="ascii") as source:
        for line in source:
            conductor, isogeny_class, number, equation, rank, torsion = line.split()
            label = conductor + isogeny_class
            curve = label + number
            if number == "1":
                class_data.append((conductor, label, rank))
            curve_data.append((curve, label, equation, torsion))

    with sqlite3.connect(output) as connection:
        connection.executescript(
            """
            CREATE TABLE t_class(
                rank INTEGER,
                class TEXT PRIMARY KEY,
                conductor INTEGER
            );
            CREATE TABLE t_curve(
                curve TEXT PRIMARY KEY,
                class TEXT,
                tors INTEGER,
                eqn TEXT UNIQUE
            );
            CREATE INDEX i_t_class_conductor ON t_class(conductor);
            CREATE INDEX i_t_curve_class ON t_curve(class);
            """
        )
        connection.executemany(
            "INSERT INTO t_class(conductor, class, rank) VALUES (?, ?, ?)",
            class_data,
        )
        connection.executemany(
            "INSERT INTO t_curve(curve, class, eqn, tors) VALUES (?, ?, ?, ?)",
            curve_data,
        )


def create_rank_databases(inputs: list[Path], outputs: list[Path]) -> None:
    """Create the text databases of elliptic curves grouped by rank."""
    allcurves, *rank_sources = inputs
    output_by_rank = {output.name: output for output in outputs}
    source_ranks = {source.name for source in rank_sources}
    if not source_ranks.issubset(output_by_rank):
        missing = ", ".join(sorted(source_ranks - output_by_rank.keys()))
        raise ValueError(f"missing output files for {missing}")

    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)

    with ExitStack() as stack:
        destinations = {
            rank: stack.enter_context(path.open("w", encoding="ascii"))
            for rank, path in output_by_rank.items()
        }
        with allcurves.open(encoding="ascii") as source:
            for line in source:
                rank = "rank" + line.split()[4]
                destinations[rank].write(line)
        for source_path in rank_sources:
            with source_path.open(encoding="ascii") as source:
                destinations[source_path.name].write(source.read())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    cremona = subparsers.add_parser("cremona")
    cremona.add_argument("allcurves", type=Path)
    cremona.add_argument("output", type=Path)

    ranks = subparsers.add_parser("ranks")
    ranks.add_argument("inputs", nargs="+", type=Path)
    ranks.add_argument("--outputs", nargs="+", required=True, type=Path)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "cremona":
        create_cremona_database(args.allcurves, args.output)
    else:
        create_rank_databases(args.inputs, args.outputs)


if __name__ == "__main__":
    main()
