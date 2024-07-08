import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

with open(args.input, "r") as f:
    lines = f.readlines()
with open(args.output, "w") as f:
    f.writelines(l.replace("static ", "") for l in lines)