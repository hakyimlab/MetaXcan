import gzip
import io
import pandas

from .. import Utilities

def try_parse(string, fail=None):
    try:
        return float(string)
    except Exception:
        return fail;

def load_data(path, key_name, value_name, white_list=None, numeric=True):
    def _ogz(p):
        return  io.TextIOWrapper(gzip.open(p, "r"), newline="")
    _o = _ogz if ".gz" in path else open
    data = {}
    c_key=None
    c_value=None
    with _o(path) as file:
        for i, line in enumerate(file):
            if i==0:
                header = line.strip().split()
                c_key = header.index(key_name)
                c_value = header.index(value_name)
                continue

            comps = line.strip().split()
            key = comps[c_key]
            if white_list:
                if not key in white_list:
                    continue

            value = comps[c_value]
            if numeric:
                value = try_parse(value)

            if value:
                data[key] = value
    return data

def load_data_column(path, column_name):
    def _ogz(p):
        return  io.TextIOWrapper(gzip.open(p, "r"), newline="")
    _o = _ogz if ".gz" in path else open
    data = []
    c_key=None
    with _o(path) as file:
        for i, line in enumerate(file):
            if i==0:
                header = line.strip().split()
                c_key = header.index(column_name)
                continue

            comps = line.strip().split()
            value = comps[c_key]
            data.append(value)
    return data

def to_data_frame(data, keys, key_column, value_column, to_numeric=None):
    ids = [x for x in keys if x in data]
    data = [(x, data[x]) for x in ids]
    if len(data) == 0:
        return pandas.DataFrame({key_column: [], value_column: []})
    data = Utilities.to_dataframe(data, [key_column, value_column], to_numeric)
    return data