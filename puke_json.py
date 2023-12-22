import json
from sys import argv


with open(argv[1]) as f:
    obj = json.load(f)

for k in obj.keys():
    print(k, obj[k])

