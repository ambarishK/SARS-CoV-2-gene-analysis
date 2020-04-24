
from numpy.random import random
import matplotlib.pyplot as plt
import time
import random
start = time.time()

class Patient():
     # default state is susceptible
     def __init__(self, state = 'susceptible'): self.state = state
     def infect(self): self.state = 'infected'
     def isolate(self): self.state = 'isolated'
     def exposed(self): self.state = 'exposed'
     def recover(self):
        self.state = 'recovered'
        if False: # set to true to explore alternative model
            if random() < 0.8:
                self.state = 'susceptible'
class PatientList():
     # create lists for each type of agent
     def __init__(self):
         self.susceptible_agents = []
         self.infected_agents = []
         self.isolated_agents = []
         self.recovered_agents = []
         self.exposed_agents = []
     def append(self, agent):
         if agent.state == 'susceptible': self.susceptible_agents.append(agent)
         elif agent.state == 'infected': self.infected_agents.append(agent)
         elif agent.state == 'recovered': self.recovered_agents.append(agent)
         elif agent.state == 'isolated': self.isolated_agents.append(agent)
         elif agent.state == 'exposed': self.exposed_agents.append(agent)
         else: print("error: must be one of the three valid states")
     def exposed(self):
         #shuffle(self.susceptible_agents) # shuffle list to random order
         perm = random.randint(0,len(self.susceptible_agents)-1)
         temp = self.susceptible_agents[-1]
         self.susceptible_agents[len(self.susceptible_agents)-1]= self.susceptible_agents[perm]
         self.susceptible_agents[perm]=temp
         patient = self.susceptible_agents.pop() # remove patient from list
         patient.exposed()
         self.append(patient) # handle appropriate list

     def infected(self):
         # shuffle(self.susceptible_agents) # shuffle list to random order
         perm = random.randint(0, len(self.exposed_agents) - 1)
         temp = self.exposed_agents[-1]
         self.exposed_agents[len(self.exposed_agents) - 1] = self.exposed_agents[perm]
         self.exposed_agents[perm] = temp
         patient = self.exposed_agents.pop()  # remove patient from list
         patient.infect()
         self.append(patient)  # handle appropriate list
     def isolate(self):
         # shuffle(self.susceptible_agents) # shuffle list to random order
         perm = random.randint(0, len(self.infected_agents) - 1)
         temp = self.infected_agents[-1]
         self.infected_agents[len(self.infected_agents) - 1] = self.infected_agents[perm]
         patient = self.infected_agents.pop()
         patient.isolate()
         self.append(patient)

     def recover(self):
         #shuffle(self.infected_agents)
         perm = random.randint(0, len(self.isolated_agents)-1)
         temp = self.isolated_agents[-1]
         self.isolated_agents[len(self.isolated_agents)-1] = self.isolated_agents[perm]
         patient = self.isolated_agents.pop()
         patient.recover()
         self.append(patient)
     def get_num_susceptible(self): return len(self.susceptible_agents)
     def get_num_infected(self): return len(self.infected_agents)
     def get_num_recovered(self): return len(self.recovered_agents)
     def get_num_isolated(self): return len(self.isolated_agents)
     def get_num_exposed(self):return len(self.exposed_agents)
     def get_num_total(self): return len(self.susceptible_agents)+len(self.infected_agents)+len(self.recovered_agents)
beta = 2.4# susceptibility rate
Isolation_rate = 0.72
sickness_rate = 0.18
gamma = 0.14 # recovery rate
susceptible_count = 3300
infected_count = 5
recovered_count = 0
isolated_count = 0
exposed_count = 40
S = [] # lists to record output
Inf = []
Is =[]
R = []
exp = []

t = []
timex = 0.0
patients = PatientList()
# create the individual patients
for indiv in range(susceptible_count):

     agent = Patient() # by default all new patients are susceptible
     patients.append(agent) # add to list

for indiv in range(infected_count):
     agent = Patient(state='infected')
     patients.append(agent)

for indiv in range(recovered_count):
     agent = Patient(state='recovered')
     patients.append(agent)

for indiv in range(isolated_count):
     agent = Patient(state='isolated')
     patients.append(agent)

for indiv in range(exposed_count):
     agent = Patient(state='exposed')
     patients.append(agent)

counter =0
while patients.get_num_infected() > 0 or patients.get_num_isolated()>0 :

    counter += 1
    for susc in range(patients.get_num_susceptible()):
        if random.random() < beta * (patients.get_num_infected() / \
            float(patients.get_num_total())):
            patients.exposed() # infect patient
    print('exposed')
    print(patients.get_num_exposed())
    for exposed in range(patients.get_num_exposed()):
        if random.random() < sickness_rate:
            patients.infected() # recover patient
    print('infec')
    print(patients.get_num_infected())
    for infected in range(patients.get_num_infected()):
        if random.random() < Isolation_rate:
            patients.isolate() # recover patient
    print(patients.get_num_isolated())
    for isolated_count in range(patients.get_num_isolated()):
        if random.random() < gamma:
            patients.recover() # recover patient
    t.append(timex)# record values for plotting
    S.append(patients.get_num_susceptible())
    Inf.append(patients.get_num_infected())
    Is.append(patients.get_num_isolated())
    R.append(patients.get_num_recovered())
    exp.append(patients.get_num_exposed())
    timex += 1 # update time

fig1 = plt.figure() # plot output
print(t)
plt.xlim(0, max(t))
plt.ylim(0, susceptible_count+infected_count+recovered_count)
plt.xlabel('time')
plt.ylabel('# patients')

#plt.plot(t, S, label="S")
plt.plot(t, Inf, label="Infected")
plt.plot(t, R, label="R")
plt.plot(t,Is, label = 'Isolated')
plt.plot(t,exp, label = 'exposed')
plt.legend()
plt.show()

end = time.time()

print(end - start)