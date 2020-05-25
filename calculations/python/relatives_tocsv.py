from calculations.python.paths import *
import ast

def load_relatives_groups(filename: str) -> {int: [[str]]}:
    result = {}
    with open(filename, 'r') as input:
        for line in input:
            line = ast.literal_eval(line)
            result[line["distance"]] = line["groups"]
    return result

def load_relatives_comparisons(filename: str):
    with open(filename, 'r') as input:
        for line in input:
            yield ast.literal_eval(line)

def relative_comparison_to_csv_line(groups: {int: [[str]]}, comparison) -> str:
    headers_compared = groups[comparison["distance_compared"]][comparison["group_id_compared"]]
    headers_reference = groups[comparison["distance_reference"]][comparison["group_id_reference"]]
    mutation = comparison["mutation"]
    return f'{repr(str(headers_compared))}, {repr(str(headers_reference))}, {mutation["position"]}, {mutation["arg"]}, {mutation["type"]}'

with open(data_path("relatives.csv"), 'w') as output:
    groups = load_relatives_groups(data_path("relatives_groups.txt"))
    print(','.join(["headers_compared", "headers_reference", "mutation_position", "mutation_arg", "mutation_type"]), file=output)
    for comparison in load_relatives_comparisons(data_path("relatives_comparisons.txt")):
        print(relative_comparison_to_csv_line(groups, comparison), file=output)