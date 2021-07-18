#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
#
# This file is a part of the ProteinDF software package.
#
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import sys
import sqlite3
import logging
logger = logging.getLogger(__name__)


class DbManager(object):
    """
    >>> db = DbManager()
    >>> db.create_table('table1', ['id', 'coulumn1', 'coulumn2'], 'id')
    >>> db.get_table_names()
    [u'table1']
    >>> db.has_table('table1')
    True

    >>> db.get_field_names('table1')
    ['id', 'coulumn1', 'coulumn2']
    >>> db.insert('table1', {'id':1, 'coulumn1':'Aichi', 'coulumn2':'Nagoya'})
    >>> db.insert('table1', {'id':2, 'coulumn1':'Miyagi', 'coulumn2':'Sendai'})
    >>> db.insert('table1', {'id':3, 'coulumn1':'Tokyo', 'coulumn2':'Tokyo'})

    >>> db.select('table1', where='coulumn1 = "Tokyo"')
    [{'coulumn1': u'Tokyo', 'coulumn2': u'Tokyo', 'id': 3}]

    >>> db.update('table1', contents={'coulumn2':'Shinjuku'}, \
        where='coulumn1 = "Tokyo"')

    >>> db.select('table1', where='coulumn1 = "Tokyo"')
    [{'coulumn1': u'Tokyo', 'coulumn2': u'Shinjuku', 'id': 3}]

    >>> db.select('table1')
    [{'coulumn1': u'Aichi', 'coulumn2': u'Nagoya', 'id': 1},\
 {'coulumn1': u'Miyagi', 'coulumn2': u'Sendai', 'id': 2},\
 {'coulumn1': u'Tokyo', 'coulumn2': u'Shinjuku', 'id': 3}]

    >>> db.delete('table1', where='id = 1')
    >>> db.select('table1')
    [{'coulumn1': u'Miyagi', 'coulumn2': u'Sendai', 'id': 2},\
 {'coulumn1': u'Tokyo', 'coulumn2': u'Shinjuku', 'id': 3}]
    """

    def __init__(self, db=':memory:', sql_debugout=False):
        self._connection = sqlite3.connect(db)
        self._cursor = self._connection.cursor()
        self._sql_debugout = sql_debugout

    def __del__(self):
        self._connection.close()

    def __getitem__(self, key):
        answer = None
        if (self.has_table(key)):
            answer = DbTable(db_manager=self,
                             table_name=key)
        return answer

    # table ==================================================================
    def create_table(self, table_name, field_names, primary_key=None):
        """
        TABLEを作成する

        table_name: テーブル名
        field_names: フィールド(カラム)名
        型を指定しない場合は、field_namesをリストにする。
        型を指定する場合は、field_namesは{名前:型}の辞書型にする。
        """
        if not self.has_table(table_name):
            fields = []
            if isinstance(field_names, dict):
                for k, v in field_names:
                    fields.append('{name} {type}'.format(k, v))
                field_names = fields
            fields_str = ', '.join(field_names)

            primary_key_str = ''
            if primary_key:
                if isinstance(primary_key, str):
                    primary_key_str = ', PRIMARY KEY({0})'.format(primary_key)
                if isinstance(primary_key, list):
                    primary_key_str += ', PRIMARY KEY ('
                    primary_key_str += ', '.join(primary_key) + ')'

            sql = 'CREATE TABLE {table_str} ({fields_str} {primary_key_str});'
            sql = sql.format(table_str=table_name,
                             fields_str=fields_str,
                             primary_key_str=primary_key_str)
            self.execute(sql)
        else:
            sys.stderr.write('already exist table: %s\n' % (table_name))

    def get_table_names(self):
        """
        TABLEの名前をリストで返す
        """
        table_names = []
        sql = "SELECT name FROM sqlite_master WHERE type='table'"
        self.execute(sql)
        results = self._cursor.fetchall()
        for row in results:
            table_name = row[0]
            table_names.append(table_name)
        return table_names

    def has_table(self, table_name):
        """
        指定されたTABLEを保持しているかどうかを返す
        """
        table_names = self.get_table_names()
        return (table_name in table_names)

    # field ==================================================================
    def get_field_names(self, table_name, fields="*"):
        """
        指定されたTABLE内のフィールドをリストで返す
        """
        field_names = None
        if (self.has_table(table_name) == True):
            field_names = []
            sql = "SELECT {0} FROM {1} LIMIT 1;".format(fields, table_name)
            self.execute(sql)
            for col, field_description in enumerate(self._cursor.description):
                field_name = field_description[0]
                field_names.append(field_name)
        return field_names

    def get_primary_keys(self, table_name):
        """
        primary key制約のフィールド名のリストを返す
        """
        answer = []
        if (self.has_table(table_name) == True):
            sql = "PRAGMA table_info({0})".format(table_name)
            self.execute(sql)

            data = self._cursor.fetchall()
            fld_info = [{}] * len(data)
            fields = []
            if self._cursor.description:
                for col, field_description in enumerate(self._cursor.description):
                    field_name = field_description[0]
                    fields.append(field_name)
                for row_index, row_data in enumerate(data):
                    entry = {}
                    for col_index, item in enumerate(row_data):
                        field_name = fields[col_index]
                        entry[field_name] = item
                    fld_info[row_index] = entry

            for info in fld_info:
                if info.get('pk') != 0:
                    answer.append(info.get('name'))

        return answer

    def insert(self, table, contents):
        """
        データレコードを追加

        contentsは、filedをキーとした辞書型
        """
        fields = []
        values = []
        for field, value in contents.items():
            fields.append(field)
            values.append(value)
        fields_str = ", ".join(fields)
        values_str = ", ".join('?' for v in values)
        sql = "INSERT INTO {table}({fields}) VALUES({values});"
        sql = sql.format(table=table,
                         fields=fields_str,
                         values=values_str)
        self.execute(sql, values)
        self._connection.commit()

    def update(self, table, contents, where):
        """
        データレコードを更新

        contents、whereは、filedをキーとした辞書型
        """
        parameters = []
        set_sections = []
        for field, value in contents.items():
            set_sections.append('%s=?' % (field))
            parameters.append(value)
        set_str = ', '.join(set_sections)

        where_str = ''
        if isinstance(where, str):
            where_str = 'WHERE ' + where
        elif isinstance(where, dict):
            where_sections = []
            for key, value in where.items():
                where_sections.append('%s=?' % (key))
                parameters.append(value)
            where_str = 'WHERE ' + ', '.join(where_sections)
        sql = 'UPDATE {table} SET {set_str} {where_str};'
        sql = sql.format(table=table,
                         set_str=set_str,
                         where_str=where_str)
        self.execute(sql, parameters)
        self._connection.commit()

    def delete(self, table, where):
        """
        データレコードを削除
        """
        parameters = []

        if isinstance(where, str):
            where_str = 'WHERE ' + where
        elif isinstance(where, dict):
            where_sections = []
            for key, value in where.items():
                where_sections.append('%s=?' % (key))
                parameters.append(value)
            where_str = 'WHERE ' + ', '.join(where_sections)

        sql = 'DELETE FROM {table} {where_str};'
        sql = sql.format(table=table,
                         where_str=where_str)
        self.execute(sql, parameters)
        self._connection.commit()

    def select(self, table, fields=None, where=None):
        """
        データを取得する
        where句はANDのみサポート
        """
        # make SQL
        parameters = []

        field_str = '*'
        if fields:
            field_str = ', '.join(fields)

        where_str = ''
        if isinstance(where, dict):
            where_sections = []
            for key, value in where.items():
                where_sections.append('{}=?'.format(key))
                parameters.append(value)
            where_str = 'WHERE ' + ' and '.join(where_sections)
        elif where != None:
            where_str = 'WHERE ' + str(where)

        sql = 'SELECT {field_str} FROM {table} {where_str};'
        sql = sql.format(table=table,
                         field_str=field_str,
                         where_str=where_str)

        # execute
        logger.debug("DB select: {}".format(sql))
        self.execute(sql, parameters)

        data = self._cursor.fetchall()
        logger.debug(data)
        answer = [{}] * len(data)
        fields = []
        if self._cursor.description:
            for col, field_description in enumerate(self._cursor.description):
                field_name = field_description[0]
                fields.append(field_name)
            for row_index, row_data in enumerate(data):
                entry = {}
                for col_index, item in enumerate(row_data):
                    field_name = fields[col_index]
                    entry[field_name] = item
                answer[row_index] = entry
        return answer

    # SQL ====================================================================
    def execute(self, sql, parameters=None):
        """
        SQLを実行する
        """
        logger.debug("sql> {0}".format(sql))
        if (parameters != None):
            return self._cursor.execute(sql, parameters)
        else:
            return self._cursor.execute(sql)

    def get_results(self, sql):
        """
        SQLを実行し、結果をリストで返す
        """
        self.execute(sql)
        data = self._cursor.fetchall()
        field_names = []
        answer = []
        if self._cursor.description:
            for col, field_description in enumerate(self._cursor.description):
                field_name = field_description[0]
                field_names.append(field_name)

            for row in data:
                row_items = {}
                for index, item in enumerate(row):
                    row_items[field_names[index]] = item

            answer.append(row_items)
        return answer

    # etc ====================================================================
    def set_user_version(self, version):
        """
        ユーザーバージョンを設定する
        """
        version = int(version)
        self.execute('PRAGMA user_version = %d;' % (version))

    def get_user_version(self):
        """
        ユーザーバージョンを返す
        """
        answer = 0
        results = self.get_results('PRAGMA user_version;')
        if results != None:
            answer = int(results[0].get('user_version', 0))
        return answer

    # output =================================================================
    def __str__(self):
        answer = ''
        tables = self.get_table_names()
        for table in tables:
            answer += self.pp_table(table)
        return answer

    def pp_table(self, table_name):
        """
        pretty print for table
        """
        answer = ''
        sql = "SELECT * FROM {0};".format(table_name)
        self.execute(sql)
        answer += 'TABLE: %s\n' % (table_name)
        answer += self.pp()
        #answer += '\n'
        return answer

    def pp(self, data=None, check_row_lengths=True):
        """
        pretty print for cursor data
        """
        if not data:
            data = self._cursor.fetchall()
        names = []
        lengths = []
        rules = []
        answer = ""
        if self._cursor.description:
            for col, field_description in enumerate(self._cursor.description):
                # print(field_description)
                field_name = field_description[0]
                names.append(field_name)
                field_length = field_description[2] or 12
                field_length = max(field_length, len(field_name))
                if check_row_lengths:
                    data_length = max([len(str(row[col])) for row in data])
                    field_length = max(field_length, data_length)
                lengths.append(field_length)
                rules.append('-' * field_length)
            format = " ".join(["%%-%ss" % l for l in lengths])
            result = [format % tuple(names), format % tuple(rules)]
            for row in data:
                result.append(format % tuple(row))
            answer = "\n".join(result)
        return answer


if __name__ == '__main__':
    import doctest
    doctest.testmod()
