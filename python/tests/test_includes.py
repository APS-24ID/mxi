"""Every source includes the standard header for what it uses.

A standard header may bring in another, so that a source using std::sort
without <algorithm> compiles with one compiler's library and not the next:
apps/mxi_index.cc did, with GCC 13 here and Ubuntu's GCC 11, and failed with
RHEL 9's -- "std::sort and std::unique not found" (Graeme's report, 1.1.1).
"""

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]

# Names, and the header that declares them.
NEEDS = {
    "<algorithm>": r"std::(sort|stable_sort|unique|min_element|max_element|"
    r"minmax_element|lower_bound|upper_bound|equal_range|binary_search|find_if|"
    r"find|count_if|count|any_of|all_of|none_of|fill|fill_n|copy_if|reverse|"
    r"rotate|partition|nth_element|partial_sort|transform|remove_if|remove|"
    r"replace|for_each|max|min|clamp|set_difference|set_intersection|merge|"
    r"is_sorted|equal|mismatch|generate|shuffle|includes|adjacent_find|search)\b",
    "<numeric>": r"std::(accumulate|iota|inner_product|partial_sum|"
    r"adjacent_difference|reduce|transform_reduce|gcd|lcm|exclusive_scan|"
    r"inclusive_scan)\b",
}

SOURCES = sorted(
    str(p.relative_to(ROOT))
    for pattern in (
        "src/*.cc",
        "src/*.hh",
        "src/*.cu",
        "apps/*.cc",
        "tests/*.cc",
        "tests/*.hh",
        "tools/*.cc",
    )
    for p in ROOT.glob(pattern)
)


@pytest.mark.parametrize("source", SOURCES)
def test_a_source_includes_what_it_uses(source):
    text = (ROOT / source).read_text()
    code = re.sub(r"//[^\n]*|/\*.*?\*/", "", text, flags=re.S)
    code = re.sub(r'"(?:[^"\\\n]|\\.)*"', '""', code)  # not inside strings
    missing = {
        header: sorted(set(m.group(0) for m in re.finditer(pattern, code)))
        for header, pattern in NEEDS.items()
        if re.search(pattern, code) and f"#include {header}" not in text
    }
    assert not missing, f"{source} uses, without their headers: {missing}"


# What every file using it needs whatever the platform: never included only
# under a condition. A script of Claude's put <numbers> under test_combine.cc's
# "#if defined(__linux__)" -- built on Linux, GCC's library bringing it in
# anyway; not on Graeme's Mac (7 November).
ALWAYS = {"<cstdint>", "<cinttypes>", "<numbers>", '"output.hh"', '"platform.hh"'}


@pytest.mark.parametrize("path", SOURCES)
def test_what_every_platform_needs_is_not_included_under_a_condition(path):
    lines = (ROOT / path).read_text().split("\n")
    conditions, under = [], []
    for n, line in enumerate(lines):
        s = line.strip()
        guard = re.match(r"#\s*ifndef\s+(\w+)", s)
        if (
            guard
            and n + 1 < len(lines)
            and re.match(
                r"#\s*define\s+" + guard.group(1) + r"\b", lines[n + 1].strip()
            )
        ):
            conditions.append(None)  # an include guard: the whole file
            continue
        if re.match(r"#\s*if", s):
            conditions.append(s)
            continue
        if re.match(r"#\s*endif", s):
            if conditions:
                conditions.pop()
            continue
        m = re.match(r'#\s*include\s*(<[^>]+>|"[^"]+")', s)
        held = [c for c in conditions if c is not None]
        if m and held and m.group(1) in ALWAYS:
            under.append(f"line {n + 1}: {m.group(1)} under {held[-1]}")
    assert not under, under
