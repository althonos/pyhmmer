import argparse
import re


_HEADER_PATTERN = re.compile(r"^@@ -(\d+),?(\d+)? \+(\d+),?(\d+)? @@.*$")


def _apply_patch(s, patch, revert=False):
    # see https://stackoverflow.com/a/40967337
    s = s.splitlines(keepends=True)
    p = patch.splitlines(keepends=True)
    
    i = 0
    while not p[i].startswith(("---", "+++")):
        i += 1
    p = p[i:]
    
    t = []
    i = 0
    sl = 0
    midx, sign = (1, "+") if not revert else (3, "-")
    while i < len(p) and p[i].startswith(("---", "+++")):
        i += 1  # skip header lines

    while i < len(p):
        match = _HEADER_PATTERN.match(p[i])
        if not match:
            raise ValueError("Invalid line in patch: {!r}".format(p[i]))
        i += 1
        l = int(match.group(midx)) - 1 + (match.group(midx + 1) == "0")
        t.extend(s[sl:l])
        sl = l
        while i < len(p) and p[i][0] != "@":
            if i + 1 < len(p) and p[i + 1][0] == "\\":
                line = p[i][:-1]
                i += 2
            else:
                line = p[i]
                i += 1
            if len(line) > 0:
                if line[0] == sign or line[0] == " ":
                    t += line[1:]
                sl += line[0] != sign

    t.extend(s[sl:])
    return "".join(t)


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-p", "--patch", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

with open(args.input, "r") as f:
    in_ = f.read()
with open(args.patch, "r") as f:
    patch = f.read()

patched = _apply_patch(in_, patch)

with open(args.output, "w") as dst:
    dst.write(patched)
