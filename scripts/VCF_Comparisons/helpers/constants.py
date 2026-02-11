from enum import Enum, auto
import os



class SETTINGS(Enum):
    STRAGLR = (auto(), True, 0, 0)
    VAMOS = (auto(), False, 0, 0)
    DEFAULT = (auto(), False, 0, 0)
    OFFSET_START = (auto(), False, -1, 0)

    def __init__(self, value, pos_only, start_offset, end_offset):
        self._value_ = value
        self.pos_only = pos_only
        self.start_offset = start_offset
        self.end_offset = end_offset

class COMP_ORDER(Enum):
    VERTICAL = auto()
    CROSS = auto()

class COMP_METHOD(Enum):
    LEVENSHTEIN = auto()
    LENGTH = auto()
    STRAGLR_LENGTH = auto()

class ORDER_METHOD(Enum):
    ASCII = auto()
    NUMERIC = auto()
