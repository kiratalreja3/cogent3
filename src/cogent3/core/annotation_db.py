"""Contains an abstract base class for the implementation
of storing sequence annotations in an in-memory SQLite database. The classes for 
GenBank and GFF files have been implemented, along with all the required helper functions"""

from __future__ import annotations

import abc
import json
import sqlite3

from typing import Optional, TypedDict

from cogent3 import open_
from cogent3.parse.genbank import MinimalGenbankParser
from cogent3.parse.gff import gff_parser

__author__ = "Kirat Alreja, Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Kirat Alreja, Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Prototype"

T = Optional[str]


class R(TypedDict):
    name: str
    type: str
    spans: list[(int, int)]


class AnnotationDbBase(abc.ABC):
    @abc.abstractmethod
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
        """returns matching records"""

    @abc.abstractmethod
    def describe(
        self,
        *,
        seq_name: bool = False,
        bio_type: bool = False,
        identifier: bool = False,
    ) -> dict:
        """returns the number of distinct values matching records"""


def _ordered_values(values: dict) -> list:

    """converts the parsed GFF dictionary values into a relevant list for the database"

    Returns:
        list: values for the in-memory database
    """

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

    """initialises an in-memory database

    Returns:
    SQLite cursor : the in-memory database

    """
    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row

    c = conn.cursor()
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
        """loads the in-memory database with values from the parsed GFF file

        Args:
            path : user-provided path to the GFF file
        """

        sql_template = "INSERT INTO GFF VALUES {}"
        sql = None

        for record in gff_parser(path):
            data = _ordered_values(dict(record))
            data[-1] = json.dumps(data[-1])

            if sql is None:
                data_placeholder = f"({','.join(('?',)*len(data))})"
                sql = sql_template.format(data_placeholder)

            self.db.execute(sql, data)

    def db_query(
        self,
        *,
        seq_name: T = None,
        bio_type: T = None,
        identifier: T = None,
        start: T = None,
        end: T = None,
        strand: T = None,
    ):
        """executes the SQL query generated by "_make_sql_query" on the in-memory database

        Args:
            seq_name (T, optional): sequence name, eg "sequence001". Defaults to None.
            bio_type (T, optional): bio type of the annotation, eg "CDS". Defaults to None.
            identifier (T, optional): the unique identifier. Defaults to None.
            start (int, optional): start coordinate value of the annotation. Defaults to None.
            end (int, optional): end coordinate value of the annotation. Defaults to None.
            strand (T, optional): strand value, unused for now. Defaults to None.

        Returns:
            SQLite rows: resulting rows from the executed query
        """

        query, values = self._make_sql_query(
            seq_name=seq_name,
            bio_type=bio_type,
            identifier=identifier,
            start=start,
            end=end,
            strand=strand,
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
        strand: T = None,
    ):
        """generates an SQL query with user-provided arguments

        Args:
            seq_name (T, optional): sequence name, eg "sequence001". Defaults to None.
            bio_type (T, optional): bio type of the annotation, eg "CDS". Defaults to None.
            identifier (T, optional): the unique identifier. Defaults to None.
            start (int, optional): start coordinate value of the annotation. Defaults to None.
            end (int, optional): end coordinate value of the annotation. Defaults to None.
            strand (T, optional): strand value, unused for now. Defaults to None.

        Raises:
            ValueError: when either bio_type or identifier is not provided

        Returns:
            (string,list): returns a tuple of the SQL query string and a list of values corresponding to the query
        """
        if not any([bio_type, identifier]):
            raise ValueError("no arguments provided")

        query = "SELECT * FROM GFF WHERE "
        clauses = []
        values = []

        if seq_name:
            clauses.append("SeqID == ?")
            values.append(seq_name)

        if bio_type:
            clauses.append("Type == ?")
            values.append(bio_type)

        if identifier:
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
        self,
        *,
        seq_name: T = None,
        bio_type: T = None,
        identifier: T = None,
        start: T = None,
        end: T = None,
        strand: T = None,
    ) -> list[R]:

        """converts the results of queries into relevant format for constructing _Annotatable objects

        Args:
        seq_name (T, optional): sequence name, eg "sequence001". Defaults to None.
        bio_type (T, optional): bio type of the annotation, eg "CDS". Defaults to None.
        identifier (T, optional): the unique identifier. Defaults to None.
        start (int, optional): start coordinate value of the annotation. Defaults to None.
        end (int, optional): end coordinate value of the annotation. Defaults to None.
        strand (T, optional): strand value, unused for now. Defaults to None.

        Returns:
            list: a list with each element that can be made into a _Annotatable object
        """
        rowdict = {}
        for row in self.db_query(
            seq_name=seq_name,
            bio_type=bio_type,
            identifier=identifier,
            start=start,
            end=end,
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

    def describe(
        self,
        *,
        seq_name: bool = False,
        bio_type: bool = False,
        identifier: bool = False,
    ) -> dict:
        """returns the unique values for user-provided arguments from the database

        Args:
            seq_name (bool, optional): _description_. Defaults to False.
            bio_type (bool, optional): _description_. Defaults to False.
            identifier (bool, optional): _description_. Defaults to False.

        Returns:
            dict: _description_
        """

        if not any([seq_name, bio_type, identifier]):
            return dict()

        value_dict = dict()

        if seq_name:
            seq_name_set = set()
            self.db.execute("SELECT DISTINCT SeqID FROM gff")
            result = self.db.fetchall()
            for row in result:
                seq_name_set.add(row["SeqID"])
            value_dict["SeqID"] = seq_name_set

        if bio_type:
            type_set = set()
            self.db.execute("SELECT DISTINCT Type FROM gff")
            result = self.db.fetchall()
            for row in result:
                type_set.add(row["Type"])
            value_dict["Type"] = type_set

        if identifier:
            identifier_set = set()
            self.db.execute("SELECT DISTINCT Attributes FROM gff")
            result = self.db.fetchall()
            for row in result:
                identifier_set.add(json.loads(row["Attributes"])["ID"])
            value_dict["identifier"] = identifier_set

        return value_dict


def _fetch_from_features(feature):

    """converts the parsed GenBank records into a revelant list for the database

    Returns:
        list: relevant values in the correct format for the database
    """

    location = [(f.first() - 1, f.last()) for f in feature["location"]]
    start = location[0][0]
    end = location[-1][1]
    strand = feature["location"].strand()
    return [
        feature["type"],
        json.dumps(location),
        feature["locus_tag"][0],
        start,
        end,
        strand,
    ]


def _make_genbank_db():

    """initialises an in-memory database

    Returns:
    SQLite cursor : the in-memory database

    """

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
                Strand integer
            )"""
    )

    return c


class GenbankAnnotationDb(AnnotationDbBase):
    def __init__(self, path):
        self.db = _make_genbank_db()
        self._populate_from_file(path)

    def _populate_from_file(self, path):

        """loads the in-memory database with values from the parsed GenBank file

        Args:
            path : user-provided path to the GenBank file
        """
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

            self.db.execute(sql, data)

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
        """executes the query constructed by "_make_sql_query" on the in-memory database

        Args:
            seq_name (T, optional): sequence name, eg "sequence001". Defaults to None.
            bio_type (T, optional): bio type of the annotation, eg "CDS". Defaults to None.
            identifier (T, optional): the unique identifier. Defaults to None.
            start (int, optional): start coordinate value of the annotation. Defaults to None.
            end (int, optional): end coordinate value of the annotation. Defaults to None.
            strand (T, optional): strand value, unused for now. Defaults to None.

        Returns:
            SQLite rows: resulting rows from the executed query
        """
        query, values = self._make_sql_query(
            seq_name=seq_name,
            bio_type=bio_type,
            identifier=identifier,
            start=start,
            end=end,
            strand=strand,
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
        strand: T = None,
    ):

        """generates an SQL query with user-provided arguments

        Args:
            seq_name (T, optional): sequence name, eg "sequence001". Defaults to None.
            bio_type (T, optional): bio type of the annotation, eg "CDS". Defaults to None.
            identifier (T, optional): the unique identifier. Defaults to None.
            start (int, optional): start coordinate value of the annotation. Defaults to None.
            end (int, optional): end coordinate value of the annotation. Defaults to None.
            strand (T, optional): strand value, unused for now. Defaults to None.

        Raises:
            ValueError: when either bio_type or identifier is not provided

        Returns:
            (string,list): returns a tuple of the SQL query string and a list of values corresponding to the query
        """

        if not any([bio_type, identifier]):
            raise ValueError("no arguments provided")

        query = "SELECT * FROM GENBANK WHERE "
        clauses = []
        values = []

        if seq_name:
            clauses.append("LocusID == ?")
            values.append(seq_name)

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
    ) -> list[R]:

        """converts the results of queries into relevant format for constructing _Annotatable objects

        Args:
        seq_name (T, optional): sequence name, eg "sequence001". Defaults to None.
        bio_type (T, optional): bio type of the annotation, eg "CDS". Defaults to None.
        identifier (T, optional): the unique identifier. Defaults to None.
        start (int, optional): start coordinate value of the annotation. Defaults to None.
        end (int, optional): end coordinate value of the annotation. Defaults to None.
        strand (T, optional): strand value, unused for now. Defaults to None.

        Returns:
            list: a list with each element that can be made into a _Annotatable object
        """

        rowdict = {}
        for row in self.db_query(
            seq_name=seq_name,
            bio_type=bio_type,
            identifier=identifier,
            start=start,
            end=end,
            strand=strand,
        ):
            id_ = row["Locus_Tag"]

            dict_value = id_ + row["Type"]

            rowdict[dict_value] = {
                "name": id_,
                "type": row["Type"],
                "spans": [tuple(l) for l in json.loads(row["Spans"])],
            }

        return list(rowdict.values())

    def describe(
        self,
        *,
        seq_name: bool = False,
        bio_type: bool = False,
        identifier: bool = False,
    ) -> dict:

        """returns the unique values for user-provided arguments from the database

        Args:
            seq_name (bool, optional): _description_. Defaults to False.
            bio_type (bool, optional): _description_. Defaults to False.
            identifier (bool, optional): _description_. Defaults to False.

        Returns:
            dict: _description_
        """

        if not any([seq_name, bio_type, identifier]):
            return dict()

        value_dict = dict()

        if seq_name:
            seq_name_set = set()
            self.db.execute("SELECT DISTINCT LocusID FROM genbank")
            result = self.db.fetchall()
            for row in result:
                seq_name_set.add(row["LocusID"])
            value_dict["LocusID"] = seq_name_set

        if bio_type:
            type_set = set()
            self.db.execute("SELECT DISTINCT Type FROM genbank")
            result = self.db.fetchall()
            for row in result:
                type_set.add(row["Type"])
            value_dict["Type"] = type_set

        if identifier:
            identifier_set = set()
            self.db.execute("SELECT DISTINCT Locus_Tag FROM genbank")
            result = self.db.fetchall()
            for row in result:
                identifier_set.add(row["Locus_Tag"])
            value_dict["identifier"] = identifier_set

        return value_dict
