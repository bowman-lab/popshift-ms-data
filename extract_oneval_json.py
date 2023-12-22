import json
from sys import argv


with open(argv[1]) as f:
    obj = json.load(f)

print(obj[argv[2]])
