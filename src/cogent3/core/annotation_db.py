from __future__ import annotations

import abc
import json
import sqlite3
from typing import Optional, TypeVar

from cogent3 import open_
from cogent3.parse.genbank import MinimalGenbankParser
from cogent3.parse.gff import gff_parser

T = Optional[str]
R = TypeVar(
    "R"
)  # define properly, a typed dict, we need to ensure appropriate keys in a dict for constructing cogent3.annotation.Feature instances


class AnnotationDbBase(abc.ABC):
    @abc.abstractmethod
    def find_records(
        self, *,
        seq_name: T = None,
        bio_type: T = None,
        identifier: T = None,
        start: int = None,
        end: int = None,
        strand: T = None,
    ) -> R:
        """returns matching records"""

    @abc.abstractmethod
    def distinct(
        self, *,
        seq_name: bool = False,
        bio_type: bool = False,
        identifier: bool = False,
    ) -> set[str]:
        """returns the number of distinct values matching records"""


def _ordered_values(values: dict) -> list:
    # doc strings for all methods and functions
    # c3dev instructions for docstrings --numpy docstring format
    # vscode extension

    keys = (
        "SeqID",
        "Source",
        "Type",
        "Start",
        "End",
        "Score",
        "Strand",
        "Phase",
        "Attributes",
    )
    return [values[k] for k in keys]


def _make_gff3_db():
    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row

    c = conn.cursor()
    #add a primary key
    c.execute(
        """ CREATE TABLE GFF (
                SeqID text,
                Source text,
                Type text,
                Start integer,
                End integer, 
                Score text,
                Strand text,
                Phase text,
                Attributes text
            
            )"""
    )

    return c


class GffAnnotationDb(AnnotationDbBase):
    def __init__(self, path):
        self.db = _make_gff3_db()
        self._populate_from_file(path)

    def _populate_from_file(self, path):
        for record in gff_parser(path):
            data = _ordered_values(dict(record))
            data[-1] = json.dumps(data[-1])
            self.db.execute(
                "INSERT INTO GFF VALUES (?,?,?,?,?,?,?,?,?)",
                data,
            )

    def db_query(
        self, *, seq_name: T = None, bio_type: T = None, identifier: T = None
    ):  # what return type
        query, values = self._make_sql_query(
            seq_name=seq_name, bio_type=bio_type, identifier=identifier
        )
        self.db.execute(query, tuple(values))
        return self.db.fetchall()

    def _make_sql_query(
        self,
        *,
        seq_name: T = None,
        bio_type: T = None,
        identifier: T = None,
        start: int = None,
        end: int = None,
    ):
        # The user must provide one of the following arguments
        if not any([bio_type, identifier]):
            raise ValueError("no arguments provided")

        query = "SELECT * FROM GFF "
        clauses = []
        values = []

        if seq_name:
            clauses.append("SeqID == ?")
            values.append(seq_name)

        if bio_type:
            clauses.append("Type == ?")
            values.append(bio_type)

        if identifier:
            # identifier can be stored as "ID" or "Parent"
            clauses.append("Attributes like ?")
            values.append(f"%{identifier}%")

        if start is not None and end is not None:
            clauses.append("Start >= ? AND End < ?")
            values.extend((start, end))
        elif start is not None:
            clauses.append("Start >= ?")
            values.append(start)
        elif end is not None:
            clauses.append("End < ?")
            values.append(end)

        return query + " AND ".join(clauses), values

    def find_records(
        self, *, seq_name: T = None, bio_type: T = None, identifier: T = None
    ) -> list[R]:
        # this needs to be a typed dict
        rowdict = {}
        for row in self.db_query(
            seq_name=seq_name, bio_type=bio_type, identifier=identifier
        ):
            attr = json.loads(row["Attributes"])
            id_ = attr["ID"]

            if id_ not in rowdict:
                rowdict[id_] = {
                    "name": id_,
                    "type": row["Type"],
                    "spans": [(row["Start"], row["End"])],
                }
            else:
                rowdict[id_]["spans"].append((row["Start"], row["End"]))

        return list(rowdict.values())

    def distinct(
        self, *,
        seq_name: bool = False,
        bio_type: bool = False,
        identifier: bool = False,
    ) -> set:

        if not any ([seq_name,bio_type,identifier]):
            return set()

        components = []

        if seq_name:
            components.append(seq_name)

        if bio_type:
            components.append(bio_type)

        if identifier:
            components.append(identifier)

        query_template = 'SELECT DISTINCT {} FROM gff ORDER BY MIN(rowid)'
        query_template = query_template.format(','.join(components))
        
        

    


def _fetch_from_features(feature):
    location = [(f.first() - 1, f.last()) for f in feature["location"]]
    start = location[0][0]
    end = location[-1][1]
    strand = feature["location"].strand()
    return [
        feature["type"],
        json.dumps(location),
        feature["locus_tag"],
        start,
        end,
        strand,
    ]


def _make_genbank_db():
    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row

    c = conn.cursor()

    c.execute(
        """ CREATE TABLE GENBANK (
                LocusID text,
                Type text,
                Spans text,
                Locus_Tag text,
                Start integer,
                End integer,
                Strand text
            )"""
    )

    return c


class GenbankAnnotationDb(AnnotationDbBase):
    def __init__(self, path):
        self.db = _make_genbank_db()
        self._populate_from_file(path)

    def _populate_from_file(self, path):
        with open_(path) as infile:
            data = list(MinimalGenbankParser(infile.readlines()))

        record = data[0]
        locus_id = record["locus"]
        sql_template = "INSERT INTO GENBANK VALUES {}"
        sql = None
        for feature in record["features"][1:]:
            if "locus_tag" not in list(feature.keys()):
                continue

            data = _fetch_from_features(feature)
            data.insert(0, locus_id)

            if sql is None:
                data_placeholder = f"({','.join(('?',)*len(data))})"
                sql = sql_template.format(data_placeholder)

            self.db.execute(
                sql,
                data,
            )

    def db_query(
        self,
        *,
        seq_name: T = None,
        bio_type: T = None,
        identifier: T = None,
        start: int = None,
        end: int = None,
        strand: T = None,
    ):
        query, values = self._make_sql_query(
            bio_type=bio_type, identifier=identifier, start=start, end=end
        )
        self.db.execute(""" SELECT LocusID FROM GENBANK LIMIT 1""")
        values.insert(0, seq_name)
        self.db.execute(query, tuple(values))
        return self.db.fetchall()

    def _make_sql_query(
        self,
        *,
        bio_type: T = None,
        identifier: T = None,
        start: int = None,
        end: int = None,
        strand: T = None,
    ):
        # The user must provide one of the following arguments
        if not any([bio_type, identifier]):
            raise ValueError("no arguments provided")

        query = "SELECT * FROM GENBANK WHERE LocusID == ? AND "
        clauses = []
        values = []

        if bio_type:
            clauses.append("Type == ?")
            values.append(bio_type)

        if identifier:
            clauses.append("Locus_Tag like ?")
            values.append(f"%{identifier}%")

        if start is not None and end is not None:
            clauses.append("Start >= ? AND End < ?")
            values.extend((start, end))
        elif start is not None:
            clauses.append("Start >= ?")
            values.append(start)
        elif end is not None:
            clauses.append("End < ?")
            values.append(end)

        return query + " AND ".join(clauses), values

    def find_records(
        self,
        *,
        seq_name: T = None,
        bio_type: T = None,
        identifier: T = None,
        start: int = None,
        end: int = None,
        strand: T = None,
    ) -> R:
        rowdict = {}
        for row in self.db_query(
            seq_name=seq_name,
            bio_type=bio_type,
            identifier=identifier,
            start=start,
            end=end,
            strand=strand,
        ):
            id_ = json.loads(row["Locus_Tag"])[0]

            dict_value = id_ + row["Type"]

            rowdict[dict_value] = {
                "name": id_,
                "type": row["Type"],
                "spans": json.loads(row["Spans"]),
            }

        return list(rowdict.values())
