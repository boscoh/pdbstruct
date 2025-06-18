import ast
import logging
import os
import sys

import click
import tqdm

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


config = read_parameters(os.path.join(this_dir, "hollow.defaults.txt"))
config.is_background = True


def init_console_logging():
    logging.basicConfig(
        stream=sys.stdout, level=logging.INFO, format="%(levelname)s: %(message)s"
    )
    config.is_background = False


def click_validate_positive(ctx, param, value):
    if value is not None and value < 0:
        raise click.BadParameter("Value must be positive.")
    return value


def tqdm_range(*args):
    return tqdm.trange(*args, disable=config.is_background)


def tqdm_iter(*args):
    return tqdm.tqdm(*args, disable=config.is_background)
