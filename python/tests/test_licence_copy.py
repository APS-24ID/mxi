"""python/LICENSE, which the wheels carry (pyproject's license-files), is the
repository's LICENSE: setuptools takes no path outside python/, so it is a
copy, and this keeps it one."""

import os

HERE = os.path.dirname(os.path.abspath(__file__))


def test_the_wheels_licence_is_the_repositorys():
    with open(os.path.join(HERE, "..", "LICENSE")) as a, open(
        os.path.join(HERE, "..", "..", "LICENSE")
    ) as b:
        assert a.read() == b.read()
