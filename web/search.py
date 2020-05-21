import sys
from calculations.python.find_tests import CovidTest, find_covid_test_in_reference, complementary
from calculations.python.paths import *

with open(data_path(REFERENCE_GENOME), 'r') as file:
    reference = ''.join(file.read().split()[1:])

for i in range(1, 4):
    sys.argv[i] = sys.argv[i].upper()
parts = find_covid_test_in_reference(CovidTest("", "", sys.argv[1], sys.argv[2], sys.argv[3]), reference)
sequence = {"F": sys.argv[1], "P": sys.argv[2], "R": sys.argv[3]}

print('<!DOCTYPE html><html><head><meta charset="UTF-8"><title>Virus test verification tool</title></head><body>')

SHOW_NEIGHBOURHOOD = 10

for part in ["F", "P", "R"]:
    print(f"Part {part} found at index {parts[part]} in the reference genome:<br>")
    print(f"Your input: <b>{sequence[part]}</b><br>")
    print(f"Found in the reference: ...{reference[parts[part] - SHOW_NEIGHBOURHOOD:parts[part]]}<b>{reference[parts[part]:parts[part] + len(sequence[part])]}</b>{reference[parts[part] + len(sequence[part]):parts[part] + len(sequence[part]) + SHOW_NEIGHBOURHOOD]}...<br>")
    print("<br>")

print(f"Please note that for the R part complementary sequence <b>{complementary(sequence['R'])}</b> has been used.<br><br>")

print("Do you want to continue with those values?<br><br>")
print('<form action="results.php" method="get">')
for part in ["F", "P", "R"]:
    print(f'\t<label for="{part}_index">Index of {part}:</label><br>')
    print(f'\t<input type="text" id="{part}_index" name="{part}_index" value="{parts[part]}"><br>')
    print(f'\t<input type="hidden" id="{part}_input" name="{part}_input" value="{sequence[part]}">')
print(f'\t<label for="hybridization_temperature">Hybridization Temperature (°C):</label><br>')
print(f'\t<input type="text" id="hybridization_temperature" name="hybridization_temperature" value="37"><br>')
print('\t<input type="submit" value="Continue">')


print('</body></html>')