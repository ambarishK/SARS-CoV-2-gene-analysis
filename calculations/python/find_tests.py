import Levenshtein
from calculations.python.paths import *
from multiprocessing import Pool, cpu_count, Value

class CovidTest:
    def __init__(self, gene: str, country: str, F: str, P: str, R: str) -> None:
        self.gene = gene
        self.country = country
        self.F = F
        self.P = P
        self.R = R
    @staticmethod
    def parse(filename: str) -> ['CovidTest']:
        with open(filename, 'r') as file:
            contents = file.read().split('>')
            if len(contents) == 0 or len(contents[0]) > 0:
                raise ValueError("Invalid header")
            contents.pop(0)
        entries = []
        if len(contents) % 3 != 0:
            raise ValueError("Invalid number of entries " + str(len(contents)))
        for entry in contents:
            entry = entry.split()
            if len(entry) == 0:
                raise ValueError("Invalid empty header")
            header = entry[0].split('_')
            if len(header) != 3:
                raise ValueError("Invalid header or incorrect size")
            header = tuple(header)
            entries.append((header, ''.join(entry[1:])))
        result = []
        expected_parts = ["F", "P", "R"]
        for i in range(0, len(entries), 3):
            if entries[i][0][0] != entries[i+1][0][0] or entries[i][0][0] != entries[i+2][0][0]:
                raise ValueError("Invalid header; different genes")
            if entries[i][0][1] != entries[i+1][0][1] or entries[i][0][1] != entries[i+2][0][1]:
                raise ValueError("Invalid header; different countries")
            parts = {}
            for j in range(i, i+3):
                for part in expected_parts:
                    if entries[j][0][2] == part:
                        parts[part] = entries[j][1]
            if len(parts) != len(expected_parts):
                raise ValueError("Too few parts")
            result.append(CovidTest(entries[i][0][0], entries[i][0][1], *[parts[p] for p in expected_parts]))
        return result
    def __str__(self):
        return str(self.__dict__)
    def __repr__(self):
        return str(self.__dict__)

def find_covid_test_in_genome(test: CovidTest, genome: str, reference: str, reference_begin: {str: int}, neighbourhood_radius: int=10, max_diff: float=0.6) -> {str: int}:
    result = {}
    prev = -1
    for part in reference_begin:
        test_part = getattr(test, part)
        reference_pre = reference[max(0, reference_begin[part] - neighbourhood_radius): reference_begin[part]]
        reference_post = reference[reference_begin[part] + len(test_part): reference_begin[part] + len(test_part) + neighbourhood_radius]
        test_data = reference_pre + test_part + reference_post
        begin = min(range(prev, len(genome)), key=lambda x: Levenshtein.distance(genome[max(0, x - neighbourhood_radius): x+len(test_part)+neighbourhood_radius], test_data))
        prev = begin
        result[part] = begin
        distance = Levenshtein.distance(genome[max(0, begin - neighbourhood_radius): begin+len(test_part)+neighbourhood_radius], test_data)
        if distance > max_diff * len(test_data):
            raise ValueError(f"Difference {distance} too big compared to length {len(test_data)} {test} part {part}")
    return result

def find_covid_test_in_reference(test: CovidTest, reference: str) -> {str: int}:
    result = {}
    prev = -1
    for part in ["F", "P", "R"]:
        test_data = getattr(test, part)
        begin = min(range(prev, len(reference)), key=lambda x: Levenshtein.distance(reference[x:x+len(test_data)], test_data))
        if begin > prev:
            prev = begin
            result[part] = begin
        else:
            raise ValueError(f"Incorrect test {test} part {part} prev {prev} begin {begin}")
    return result

def find_tests(genome: str, tests: [CovidTest], reference: str, tests_in_reference: [{str: int}]) -> [{str: int}]:
    result = []
    for test, ref in zip(tests, tests_in_reference):
        result.append(find_covid_test_in_genome(test, genome, reference, ref))
    return result

def tests_to_print(header: str, tests_results: [{str: int}], tests: [CovidTest]):
    return {"header": header, "tests": [{"name": test.gene + "_" + test.country + "_" + part, "begin": result[part], "end": result[part] + len(getattr(test, part))} for result, test in zip(tests_results, tests) for part in result]}

def load_genomes():
    genomes = []
    with open(data_path(FILTERED_GENOMES), 'r') as file:
        next(file)
        header = None
        contents = ''
        for line in file:
            if line.startswith(">"):
                if header is not None:
                    genomes.append((header, contents))
                contents = ''
                header = line[:-1]
                continue
            else:
                contents += line.rstrip()
        genomes.append((header, contents))
    return genomes

tests = CovidTest.parse(data_path("primers_public-3.fas"))

with open(data_path(REFERENCE_GENOME), 'r') as file:
    reference = ''.join(file.read().split()[1:])
tests_in_reference = [find_covid_test_in_reference(test, reference) for test in tests]

genomes = load_genomes()
print(f"{len(genomes)} genomes loaded.", flush=True)

counter = None

def init(arg):
    global counter
    counter = arg

def perform_tests(genome):
    header, contents = genome
    global counter
    r = (header, find_tests(contents, tests, reference, tests_in_reference))
    with counter.get_lock():
        counter.value += 1
        c = counter.value
    if int(c * 100 / len(genomes)) != int((c + 1) * 100 / len(genomes)):
        print(str(int((c + 1) * 100 / len(genomes))) + "%", flush=True)
    return r

try:
    pool = Pool(cpu_count(), initializer=init, initargs=(Value('i', 0),))
    result = pool.map(perform_tests, genomes)
finally:
    pool.close()
    pool.join()

with open(data_path(TESTS_RESULTS), 'w') as file:
    print(len(genomes), file=file)
    for header, test_results in result:
        print(tests_to_print(header, test_results, tests), file=file)