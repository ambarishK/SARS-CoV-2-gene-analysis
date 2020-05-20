import sys
from calculations.python.find_tests import CovidTest, find_covid_test_in_reference
from calculations.python.paths import *

with open(data_path(REFERENCE_GENOME), 'r') as file:
    reference = ''.join(file.read().split()[1:])

parts = find_covid_test_in_reference(CovidTest("", "", sys.argv[1], sys.argv[2], sys.argv[3]))
input = {"F": sys.argv[1], "P": sys.argv[2], "R": sys.argv[3]}

print('<!DOCTYPE html><html><head><meta charset="UTF-8"><title>Virus test verification tool</title></head><body>')

SHOW_NEIGHBOURHOOD = 10

for part in ["F", "P", "R"]:
    print(f"Part {part} found at index {parts[part]} in the reference genome:<br>")
    print(f"Your input: <b>{input[part]}</b><br>")
    print(f"Found in the reference: ...{reference[parts[part] - SHOW_NEIGHBOURHOOD:parts[part]]}<b>{reference[parts[part]:parts[part] + len(input[part])]}</b>{reference[parts[part] + len(input[part]):parts[part] + len(input[part]) + SHOW_NEIGHBOURHOOD]}...<br>")
    print("<br>")

print("Please note that for the R part complementary sequence has been used.<br>")

print('<form action="results.php" method="get">')
for part in ["F", "P", "R"]:
    print(f'\t<label for="{part}">Index of {part}:</label><br>')
    print(f'\t<input type="text" id="{part}" name="{part}"><br>')
print('\t<input type="submit" value="Submit">')


print('</body></html>')