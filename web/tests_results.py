import sys, subprocess
import statistics
from web.chart_data import GOOGLE_CHARTS_COUNTRIES_NAMES
from calculations.python.find_tests import load_genomes as load_genomes_raw, CovidTest, find_covid_test_in_genome, CovidTestPartResult
from calculations.python.paths import *
from multiprocessing import Pool, cpu_count
import numpy as np
from web.mutations import icolorbar2

np.seterr(divide='ignore', invalid='ignore')

sequence = {"F": sys.argv[1], "P": sys.argv[2], "R": sys.argv[3]}
index = {"F": int(sys.argv[4]), "P": int(sys.argv[5]), "R": int(sys.argv[6])}
hybridization_temperature = float(sys.argv[7].replace(",", "."))

SEARCH_RADIUS = 20
LOG_HG_RATIO_LEN = 6

def load_genomes():
    required_len = max([v for _, v in index.items()]) + SEARCH_RADIUS - 1 + max(LOG_HG_RATIO_LEN, max((len(v) for _, v in sequence.items())))
    return ((h, c) for h, c in load_genomes_raw() if len(c) >= required_len)

with open(data_path(REFERENCE_GENOME), 'r') as file:
    reference = ''.join(file.read().split('\n')[1:])

test = CovidTest("Your", "Test", sequence["F"], sequence["P"], sequence["R"])

results_dictionary = {}
def get_free_energy(test: CovidTest, genome: str, test_result: {str: CovidTestPartResult}, temperature: float) -> {str: float}:
    global results_dictionary
    result = {}
    for test_part in ["F", "P", "R"]:
        l = getattr(test, test_part)
        r = genome[test_result[test_part].begin: test_result[test_part].end]
        if (l, r) in results_dictionary:
            value = results_dictionary[(l, r)]
        else:
            value = float(subprocess.check_output(f'UNAFOLDDAT=data/oligoarrayaux {data_path("hybrid-min")} -q -t {temperature} -T {temperature} {l} {r}',shell=True).split()[0])
            results_dictionary[(l, r)] = value
        result[test_part] = value
    return result

def perform_test_and_get_free_energy(genome: (str, str)) -> ((str, {str: CovidTestPartResult}), {str: float}):
    global test
    tr = find_covid_test_in_genome(test, genome[1], reference, index, 0, {test_part: (max(index[test_part] - SEARCH_RADIUS, 0), min(len(genome[1]), index[test_part] + SEARCH_RADIUS)) for test_part in ["F", "P", "R"]})[0]
    return ((genome[0], tr), get_free_energy(test, genome[1], tr, hybridization_temperature))

try:
    pool = Pool(cpu_count())
    test_results_ = pool.map(perform_test_and_get_free_energy, load_genomes())
finally:
    pool.close()
    pool.join()

test_results = [tr for tr, _ in test_results_]
free_energies = [dg for _, dg in test_results_]

del test_results_

LOG_LH_RATIO_BETA_0 = -5.62
LOG_LH_RATIO_BETA_1 = -1.55
LOG_LH_RATIO_BETA_2 = 0.33
LOG_LH_RATIO_BETA_3 = 0.18

def get_log_lh_ratio(genome: str, test: CovidTest, test_result: {str: CovidTestPartResult}, free_energy: {str: float}) -> {str: float}:
    def calc(test_part: str, gen_part: str, free_energy: float):
        test_part = test_part[-LOG_HG_RATIO_LEN:]
        gen_part = gen_part[-LOG_HG_RATIO_LEN:]
        ix = 0
        for i in reversed(range(LOG_HG_RATIO_LEN)):
            if test_part[i] != gen_part[i]:
                ix = i + 1
                break
        return LOG_LH_RATIO_BETA_0 + LOG_LH_RATIO_BETA_1 * free_energy + LOG_LH_RATIO_BETA_2 * ix + LOG_LH_RATIO_BETA_3 * free_energy * ix
    return {"F": calc(test.F, genome[test_result["F"].begin: test_result["F"].end], free_energy["F"]), "R": calc(test.R[::-1], genome[test_result["R"].begin: test_result["R"].end][::-1], free_energy["F"])}

log_lh_ratios = [get_log_lh_ratio(g[1], test, tr[1], dG) for g, tr, dG in zip(load_genomes(), test_results, free_energies)]

'''
def get_csv_header(test: CovidTest):
    ret = ["header"]
    for part in ["F", "P", "R"]:
        for e in ["begin", "end", "mutations", "relaxed-mutations", "free-energy", "log-lh-ratio"]:
            ret.append(f'{test.gene}_{test.country}_{part}_{e}')
    ret.remove(f'{test.gene}_{test.country}_P_log-lh-ratio')
    return ret

def get_csv_line(header: str, test_result: {str: CovidTestPartResult}, free_energy: {str: float}, log_lh_ratio: {str: float}) -> [str]:
    ret = [header]
    for part, tr in test_result.items():
        ret.append(str(tr.begin))
        ret.append(str(tr.end))
        ret.append(repr(repr(tr.mutations)))
        ret.append(repr(repr(tr.relaxed_mutations)))
        ret.append(str(free_energy[part]))
        if part != "P":
            ret.append(str(log_lh_ratio[part]))
    return ret

print(','.join(get_csv_header(test)))

for tr, free_energy, log_lh_ratio in zip(test_results, free_energies, log_lh_ratios):
    header, test_result = tr
    print(','.join(get_csv_line(header, test_result, free_energy, log_lh_ratio)))
'''


def prob_part(log_lh: float):
    try:
        return 1 / (2.7182 ** (-log_lh) + 1)
    except OverflowError:
        return 0.0 #(2.7182 ** (-log_lh) + 1) was too big, 1/<big number> is about zero

def get_countries_probabilities(test_results: [(str, {str: CovidTestPartResult})], log_lh_ratios: [{str: float}]) -> {str: float}:
    ret = {}
    for tr, log_lh in zip(test_results, log_lh_ratios):
        header, _ = tr
        country = header.split('/')[1]
        if country in GOOGLE_CHARTS_COUNTRIES_NAMES:
            country = GOOGLE_CHARTS_COUNTRIES_NAMES[country]
        else:
            country = GOOGLE_CHARTS_COUNTRIES_NAMES[header.split('/')[2]]
        if country not in ret:
            ret[country] = []
        ret[country].append(prob_part(log_lh["F"] * log_lh["R"]))
    for country in ret:
        ret[country] = statistics.mean(ret[country])
    return ret

HISTOGRAM_BUCKETS = 11

print('<!DOCTYPE html><html><head><meta charset="UTF-8"><title>Virus test verification tool</title>')
print('<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>')
print('<script type="text/javascript">')
print('')
print('let map_data = [["Country", "Percent"]', end='')
for country, probability in get_countries_probabilities(test_results, log_lh_ratios).items():
    print(f', [{repr(country)}, {probability}]')
print('];')
print('let histograms_data = {')
probabilities = {"F": [prob_part(r["F"]) for r in log_lh_ratios], "R": [prob_part(r["R"]) for r in log_lh_ratios]}
probabilities["Total"] = [a * b for a, b in zip(probabilities["F"], probabilities["R"])]
for graph in probabilities:
    probabilities[graph] = (np.histogram(probabilities[graph], np.linspace(0.0, 1.0, HISTOGRAM_BUCKETS))[0]*100/len(probabilities[graph])).tolist()
histogram_buckets = np.linspace(0.0, 1.0, HISTOGRAM_BUCKETS)
for graph, data in probabilities.items():
    print(f'{repr(graph)}: [["Bucket", "Probability"],')
    for i, probability in enumerate(data):
        print(f'["{round(histogram_buckets[i], 2)} - {round(histogram_buckets[i + 1], 2)}", {probability}],')
    print('],')
print('};')
print('''
google.charts.load('current', {
    'packages': ['geochart', 'corechart'],
});
google.charts.setOnLoadCallback(drawMap);
google.charts.setOnLoadCallback(drawHistograms);
function drawMap() {
    let data = google.visualization.arrayToDataTable(map_data);
    let options = {
        colorAxis: {minValue: 0,  maxValue: 1, colors: ['#ff0000', '#00ff00']},
        minValue: 0,
        maxValue: 1,
        backgroundColor: '#81d4fa',
        datalessRegionColor: '#ffffff',
    };
    let chart = new google.visualization.GeoChart(document.getElementById('map'));
    chart.draw(data, options);
}
function drawHistograms() {
    let names = {F: "forward primer", R: "reverse primer", Total: "test"}
    Object.keys(histograms_data).forEach(function(graph) {
        let data = google.visualization.arrayToDataTable(histograms_data[graph]);
        let options = {
            width: 600,
            height: 500,
            chartArea: {width: '80%', height: '60%'},
            vAxis: {minValue: 0, maxValue: 100, title: "% of genomes", format: "#'%'"},
            legend: {position: "none"},
            hAxis: {slantedText: true, slantedTextAngle: 70, title: "Estimated likelihood of amplification for " + names[graph]},
            title: "Results of " + names[graph],
        };
        let chart = new google.visualization.ColumnChart(document.getElementById('hist_' + graph));
        chart.draw(data, options);
    });
}
''')
print('</script>')
print('</head><body>')
print('<table><tr>')
print('<td><div id="hist_F"></div></td>')
print('<td><div id="hist_Total"></div></td>')
print('<td><div id="hist_R"></div></td>')
print('</tr></table>')
print('<div id="map" style="width: 900px; height: 500px; text-align: center;"></div>')

def calc_hist(test_results: [(str, {str: CovidTestPartResult})]) -> {str: [int]}:
    histogram_dict = {"F": {}, "P": {}, "R": {}}
    for _, test_result in test_results:
        for test_part in ["F", "P", "R"]:
            for mutation in test_result[test_part].mutations:
                if mutation.position in histogram_dict[test_part]:
                    histogram_dict[test_part][mutation.position] += 1
                else:
                    histogram_dict[test_part][mutation.position] = 1
    return {tp: [data[i] if i in data else 0 for i in range(max(data))] for tp, data in histogram_dict.items()}

icolorbar2(calc_hist(test_results))
print('''
<script type="text/javascript">
for(let element of document.getElementsByClassName("updatemenu-header")) {
	element.__onclick()
	document.getElementsByClassName("updatemenu-dropdown-button")[0].__onclick(0)
}
</script>
''')
print('</body></html>')
