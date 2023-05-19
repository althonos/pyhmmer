import argparse
import itertools
import json
import os
import re
import math

import numpy
import matplotlib.pyplot as plt
from palettable.cartocolors.qualitative import Bold_3

plt.rcParams["svg.fonttype"] = "none"


parser = argparse.ArgumentParser()
parser.add_argument("--folder", required=True)
args = parser.parse_args()




with open(os.path.join(args.folder, "hmmsearch.json")) as f:
    data = json.load(f)

CPU_RX = re.compile(r"(?:--cpu|--jobs) (\d+)")
for result in data["results"]:
    if result["command"].startswith(("hmmsearch", "hmmscan")):
        result["tool"] = result["command"].split(" ")[0]
    elif "hmmscan" in result["command"]:
        result["tool"] = "pyhmmer.hmmscan"
    elif "hmmsearch" in result["command"]:
        result["tool"] = "pyhmmer.hmmsearch"
    else:
        raise ValueError(f"Unknown command: {result['command']!r}")
    if "Pfam.v33-1.hmm " in result["command"]:
        result["format"] = "text"
    elif "Pfam.v33-1.pressed.hmm.h3m " in result["command"]:
        result["format"] = "binary"
    elif "Pfam.v33-1.pressed.hmm " in result["command"]:
        result["format"] = "pressed"
    else:
        raise ValueError(f"Unknown format: {result['command']!r}")
    result["cpu"] = int(CPU_RX.search(result["command"]).group(1))

plt.figure(1, figsize=(12, 12))
plt.subplot(2, 2, 1)

data["results"].sort(key=lambda r: (r["tool"], r["format"], r["cpu"]))
for color, (tool, group) in zip(
    Bold_3.hex_colors[1:], 
    filter(
       lambda x: "hmmsearch" in x[0], 
       itertools.groupby(data["results"], key=lambda r: r["tool"])
    )
):
    group = list(group)

    group_text = [r for r in group if r["format"] == "text"]
    if any(group_text):
        X = numpy.array([r["cpu"] for r in group_text])
        Y = numpy.array([r["mean"] for r in group_text])
        ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group_text]
        plt.plot(X, Y, color=color, linestyle="--", marker="o")
        plt.fill_between(X, Y - ci, Y + ci, color=color, alpha=0.1)

    group_pressed = [r for r in group if r["format"] == "pressed"]
    if any(group_pressed):
        X = numpy.array([r["cpu"] for r in group_pressed])
        Y = numpy.array([r["mean"] for r in group_pressed])
        ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group_pressed]
        plt.plot(X, Y, label=tool, color=color, marker="o")
        plt.fill_between(X, Y - ci, Y + ci, color=color, alpha=0.1)


plt.legend()
plt.xlabel("CPUs")
plt.ylabel("Time (s)")

plt.figure(1, figsize=(24, 6))
plt.subplot(2, 2, 2)

data["results"].sort(key=lambda r: (r["tool"], r["format"], r["cpu"]))
for color, (tool, group) in zip(
    Bold_3.hex_colors[1:], 
    filter(
       lambda x: "hmmscan" in x[0], 
       itertools.groupby(data["results"], key=lambda r: r["tool"])
    )
):
    group = list(group)

    group_text = [r for r in group if r["format"] == "text"]
    if any(group_text):
        X = numpy.array([r["cpu"] for r in group_text])
        Y = numpy.array([r["mean"] for r in group_text])
        ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group_text]
        plt.plot(X, Y, color=color, linestyle="--", marker="o")
        plt.fill_between(X, Y - ci, Y + ci, color=color, alpha=0.1)

    group_pressed = [r for r in group if r["format"] == "pressed"]
    if any(group_pressed):
        X = numpy.array([r["cpu"] for r in group_pressed])
        Y = numpy.array([r["mean"] for r in group_pressed])
        ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group_pressed]
        plt.plot(X, Y, label=tool, color=color, marker="o")
        plt.fill_between(X, Y - ci, Y + ci, color=color, alpha=0.1)


plt.legend()
plt.xlabel("CPUs")
plt.ylabel("Time (s)")


with open(os.path.join(args.folder, "phmmer.json")) as f:
    data = json.load(f)

CPU_RX = re.compile(r"(?:--cpu|--jobs) (\d+)")
for result in data["results"]:
    if result["command"].startswith("phmmer"):
        result["tool"] = result["command"].split(" ")[0]
    else:
        result["tool"] = "pyhmmer.phmmer"
    result["cpu"] = int(CPU_RX.search(result["command"]).group(1))

# plt.figure(2, figsize=(6,6))
plt.subplot(2, 2, 3)
data["results"].sort(key=lambda r: (r["tool"], r["cpu"]))
for color, (tool, group) in zip(
    Bold_3.hex_colors[1:], itertools.groupby(data["results"], key=lambda r: r["tool"])
):
    group = list(group)
    X = numpy.array([r["cpu"] for r in group])
    Y = numpy.array([r["mean"] for r in group])
    ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group]
    plt.plot(X, Y, label=tool, color=color, marker="o")
    plt.fill_between(X, Y - ci, Y + ci, color=color, alpha=0.1)

plt.legend()
plt.xlabel("CPUs")
plt.ylabel("Time (s)")



with open(os.path.join(args.folder, "nhmmer.json")) as f:
    data = json.load(f)

CPU_RX = re.compile(r"(?:--cpu|--jobs) (\d+)")
for result in data["results"]:
    if result["command"].startswith("nhmmer"):
        result["tool"] = result["command"].split(" ")[0]
    else:
        result["tool"] = "pyhmmer.nhmmer"
    result["cpu"] = int(CPU_RX.search(result["command"]).group(1))

# plt.figure(2, figsize=(6,6))
plt.subplot(2, 2, 4)
data["results"].sort(key=lambda r: (r["tool"], r["cpu"]))
for color, (tool, group) in zip(
    Bold_3.hex_colors[1:], itertools.groupby(data["results"], key=lambda r: r["tool"])
):
    group = list(group)
    X = numpy.array([r["cpu"] for r in group])
    Y = numpy.array([r["mean"] for r in group])
    ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group]
    plt.plot(X, Y, label=tool, color=color, marker="o")
    plt.fill_between(X, Y - ci, Y + ci, color=color, alpha=0.1)

plt.legend()
plt.xlabel("CPUs")
plt.ylabel("Time (s)")

plt.tight_layout()
plt.savefig(os.path.join(args.folder, "plot.svg"), transparent=True)
plt.show()
