import time
from datetime import datetime
from pytz import timezone
import pytz
from functools import wraps
from typing import Union, List, Tuple


def timer(f):
    "Decorator used for timing functions."

    @wraps(f)
    def wrap(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        print(f"`{f.__name__}` took: {te-ts:2.4f} sec")
        return result

    return wrap


def timing(f):
    "Decorator used for timing functions and also displaying arguments"

    @wraps(f)
    def wrap(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        print(f"func:{f.__name__} args:[{args, kw}] took: {te-ts:2.4f} sec")
        return result

    return wrap


def timestamped(print_pdt: bool = False, date_format="%Y-%m%d"):
    date = datetime.now(tz=pytz.utc)
    date = date.astimezone(timezone("US/Pacific"))

    if print_pdt:
        print("append HMS")
        date_fmt = date_format + "_%H%M%S"
        STAMP = date.strftime(date_fmt) + "PDT"
    else:
        date_fmt = date_format
        STAMP = date.strftime(date_fmt)
    return STAMP


def stringify_list(items_list: list):
    """
    converts list of items to string separated by a vertical line `|`
    """
    return "|".join(items_list)


def hyphen_to_underscore(names: Union[str, List[str]]):
    if isinstance(names, str):
        return names.replace("-", "_")
    else:
        return [name.replace("-", "_") for name in names]


def underscore_to_hyphen(names: Union[str, List[str]]):
    if isinstance(names, str):
        return names.replace("_", "-")
    else:
        return [name.replace("_", "-") for name in names]


def extract_neighbors_leiden(string: str) -> Tuple[int, float]:
    "Extract n_neighbors and leiden resolution out of string name."
    s = string.split("_")
    return (int(s[1]), float(s[3]))
