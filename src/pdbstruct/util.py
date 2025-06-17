import ast
import os


this_dir = os.path.dirname(__file__)

class AttrDict(dict):
    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            raise AttributeError(item)

    def __setattr__(self, key, value):
        self[key] = value


def read_parameters(fname):
    with open(fname, "r") as f:
        return AttrDict(ast.literal_eval(f.read()))
