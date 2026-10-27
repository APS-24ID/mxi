#!/bin/sh
# The complete source for a release: the repository at the tag and each
# submodule at its pinned commit, in one archive -- GitHub's own source
# archives leave the submodules out, and without them nothing builds. Also the
# source a conda-forge recipe can point at.
#
#   scripts/source_archive.sh v1.0.6     # writes mxi-1.0.6-source.tar.gz
set -eu
tag=${1:?usage: source_archive.sh TAG}
version=${tag#v}
name="mxi-${version}"
root=$(git rev-parse --show-toplevel)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cd "$root"
git archive --format=tar --prefix="${name}/" "$tag" | tar -x -C "$work"
# Each submodule at the commit the tag pins, from the submodule's own checkout.
git ls-tree -r "$tag" | awk '$2 == "commit" {print $3, $4}' |
while read -r commit path; do
    if [ ! -d "$path/.git" ] && [ ! -f "$path/.git" ]; then
        echo "submodule $path is not checked out: git submodule update --init --recursive" >&2
        exit 1
    fi
    (cd "$path" && git archive --format=tar --prefix="${name}/${path}/" "$commit") |
        tar -x -C "$work"
done
tar -C "$work" -czf "$root/${name}-source.tar.gz" "$name"
echo "${name}-source.tar.gz"
