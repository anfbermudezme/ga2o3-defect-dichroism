from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterator, Optional

import matplotlib as mpl


@contextmanager
def mpl_rc(rc: Optional[Dict[str, object]] = None) -> Iterator[None]:
    """A tiny helper around matplotlib.rc_context.

    This keeps plotting functions from permanently mutating global rcParams.
    """
    if rc is None:
        yield
        return
    with mpl.rc_context(rc):
        yield


def ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
