import sys, subprocess
from typing import Iterable

from calculations.python.find_tests import load_genomes, CovidTest, find_covid_test_in_genome, CovidTestPartResult, complementary, tests_to_print
from calculations.python.paths import *
from multiprocessing import Pool, cpu_count

sequence = {"F": sys.argv[1], "P": sys.argv[2], "R": sys.argv[3]}
index = {"F": int(sys.argv[4]), "P": int(sys.argv[5]), "R": int(sys.argv[6])}
hybridization_temperature = float(sys.argv[7].replace(",", "."))

SEARCH_RADIUS = 20

with open(data_path(REFERENCE_GENOME), 'r') as file:
    reference = ''.join(file.read().split()[1:])

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
    tr = find_covid_test_in_genome(test, genome[1], reference, index, 0,
                                       {test: (index[test] - SEARCH_RADIUS, index[test] + SEARCH_RADIUS) for test in
                                        ["F", "P", "R"]})[0]
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
        test_part = test_part[-6:]
        gen_part = gen_part[-6:]
        ix = 0
        for i in reversed(range(6)):
            if test_part[i] != gen_part[i]:
                ix = i + 1
                break
        return LOG_LH_RATIO_BETA_0 + LOG_LH_RATIO_BETA_1 * free_energy + LOG_LH_RATIO_BETA_2 * ix + LOG_LH_RATIO_BETA_3 * free_energy * ix

    return {"F": calc(test.F, genome[test_result["F"].begin: test_result["F"].end], free_energy["F"]), "R": calc(test.R[::-1], genome[test_result["R"].begin: test_result["R"].end][::-1], free_energy["F"])}

log_lh_ratios = [get_log_lh_ratio(g[1], test, tr[1], dG) for g, tr, dG in zip(load_genomes(), test_results, free_energies)]

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
