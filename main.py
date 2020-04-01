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

def solve(nodes, pop, t=0, params=0,snapshots=1):
    
    dt = 1
    d_t_minus = snapshots[0]
    for d_t in snapshots[1:]:
        for node in nodes:  
            dS, dE, dI = derivatives(pop[node+' S'][d_t_minus], pop[node+' E'][d_t_minus], pop[node+' I'][d_t_minus], params)
            pop[node+' S'][d_t] = pop[node+' S'][d_t_minus]+dt*dS
            pop[node+' E'][d_t] = pop[node+' E'][d_t_minus]+dt*dE
            pop[node+' I'][d_t] = pop[node+' I'][d_t_minus]+dt*dI
            pop[node+' R'][d_t] = 1-pop[node+' I'][d_t]-pop[node+' S'][d_t]-pop[node+' E'][d_t]
        migration(pop,params)
        d_t_minus = d_t
    return pop
    
def initialise_pop(nodes=np.array(['0']),snapshots=np.array(['0']),total_pop={'0':100000.}):
    #initialise populations of different nodes
    
    index = ['S','E','I','R']
    pop_index={'S':1.,'E':0.,'I':1.,'R':0}
    nodes_t = np.array([])
    
    for i in index:
        nodes_t=np.append(nodes_t,[node+ ' ' +i for node in nodes])
    pop = pd.DataFrame(columns=nodes_t,index=snapshots)
    for node in nodes:
        pop[node +' '+ 'S'].loc[snapshots[0]]=pop_index['S']-pop_index['I']/total_pop[node]
        pop[node +' '+ 'E'].loc[snapshots[0]]=pop_index['E']
        pop[node +' '+ 'I'].loc[snapshots[0]]=pop_index['I']/total_pop[node]
        pop[node +' '+ 'R'].loc[snapshots[0]]=pop_index['R']
        
    return pop    

def plot_results(pop,snapshots):
    for node in nodes:
        fig, ax = plt.subplots()
        ax.plot(snapshots,pop[node +' '+ 'E'],label='Exposed')
        ax.plot(snapshots,pop[node +' '+ 'S'],label='Susceptible')
        ax.plot(snapshots,pop[node +' '+ 'I'],label='Infected')
        ax.plot(snapshots,pop[node +' '+ 'R'],label='Recovered')

        ax.set_xlim([snapshots[0],snapshots[len(snapshots)-1]])
        ax.set_title(node)
        plt.xticks(rotation=90)
        plt.legend()

        plt.show()
    
def migration(pop,params):
    #internal and external migration
    
    return pop

#main
params = load_params()

snapshots = pd.date_range('2020-01-01', periods=500, freq='d')

nodes=['Berlin','Bayern','Baden Württemberg']
pop=initialise_pop(nodes=nodes,snapshots=snapshots,total_pop={'Berlin':100000.,'Bayern':100000.,'Baden Württemberg':100000.})


pop=solve(nodes=nodes,pop=pop,params=params,snapshots=snapshots)

plot_results(pop,snapshots)
