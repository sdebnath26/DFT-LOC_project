from __future__ import annotations

from pprint import pformat
from typing import Any


def print_dict(data: dict[str, Any]) -> None:
    for key, value in data.items():
        print(f"{key}: {pformat(value)}")
