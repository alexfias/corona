import pandas as pd
import matplotlib.pyplot as plt
import yaml
import numpy as np

def load_params():
    
    with open('params.yaml') as file:
        
        params = yaml.load(file)
        for key in params:
            params[key]=float(params[key])
        print(params)
        
    return params

def derivatives(S, E, I, params):
    
    #compute derivatives 
    nu = params['beta']/params['R0']
    
    dS = params['b']*(1.-S)-params['beta']*S*I
    dE = params['beta']*S*I-params['mu']*E
    dI = params['mu']*E-nu*I-params['b']*I
    
    return dS, dE, dI

def solve(nodes=nodes, pop=pop, t=0, params=0,snapshots=1):
    S=pop[node + ' S']
    E=pop[node + ' E']
    I=pop[node + ' I']
    R=pop[node + ' R']
    dt = 1

    for d_t in range(len(snapshots)-1):
        
        dS, dE, dI = derivatives(S[d_t], E[d_t], I[d_t], params)
        S = np.append(S,[S[d_t]+dt*dS])
        E = np.append(E,[E[d_t]+dt*dE])
        I = np.append(I,[I[d_t]+dt*dI])
        R = np.append(R,[1-I[d_t]-S[d_t]-E[d_t]])

    
    return S,E,I,R
    
def initialise_pop(nodes=np.array(['0']),snapshots=np.array(['0']),total_pop={'0':100000.}):
    #initialise populations of different nodes
    
    N=100000.
    index = ['S','E','I','R']
    pop_index={'S':1.,'E':0.,'I':1.,'R':0}
    nodes_t = np.array([])
    
    for i in index:
        nodes_t=np.append(nodes_t,[node+ ' ' +i for node in nodes])
    pop = pd.DataFrame(columns=nodes_t,index=snapshots)
    for node in nodes:
        pop[node +' '+ 'S'].loc[snapshots[0]]=pop_index['S']
        pop[node +' '+ 'E'].loc[snapshots[0]]=pop_index['E']
        pop[node +' '+ 'I'].loc[snapshots[0]]=pop_index['I']/total_pop[node]
        pop[node +' '+ 'E'].loc[snapshots[0]]=pop_index['R']
        
    return pop    

def plot_results(S,E,I,R,snapshots):
    
    fig, ax = plt.subplots()
    ax.plot(snapshots,E,label='Exposed')
    ax.plot(snapshots,S,label='Susceptible')
    ax.plot(snapshots,I,label='Infected')
    ax.plot(snapshots,R,label='Recovered')

    ax.set_xlim([snapshots[0],snapshots[len(snapshots)-1]])
    plt.xticks(rotation=90)
    plt.legend()

    plt.show()



#main
params = load_params()

snapshots = pd.date_range('2020-01-01', periods=500, freq='d')

nodes=['0']
pop=initialise_pop(nodes=nodes)

for node in nodes:
    S,E,I,R=solve(pop=pop,params=params,snapshots=snapshots)

plot_results(S,E,I,R,snapshots)
