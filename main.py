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

def solve(S=[0], E=[0], I=[0], R=[0], t=0, params=0,snapshots=1):
    
    dt = 1

    for d_t in range(len(snapshots)-1):
        
        dS, dE, dI = derivatives(S[d_t], E[d_t], I[d_t], params)
        S = np.append(S,[S[d_t]+dt*dS])
        E = np.append(E,[E[d_t]+dt*dE])
        I = np.append(I,[I[d_t]+dt*dI])
        R = np.append(R,[1-I[d_t]-S[d_t]-E[d_t]])

    
    return S,E,I,R
    
    
#main
params = load_params()

snapshots = pd.date_range('2020-01-01', periods=500, freq='d')

#initial values with a fully susceptible population and one infected person
S = np.array([1.])
E = np.array([0.])
I = np.array([1./params['N']])
R = np.array([0.])

S,E,I,R=solve(S=S,E=E,I=I,R=R,params=params,snapshots=snapshots)

print(S)
print(I)
print(E)
print(R)

fig, ax = plt.subplots()
ax.plot(snapshots,E,label='Exposed')
ax.plot(snapshots,S,label='Susceptible')
ax.plot(snapshots,I,label='Infected')
ax.plot(snapshots,R,label='Recovered')

ax.set_xlim([snapshots[0],snapshots[len(snapshots)-1]])
plt.xticks(rotation=90)
plt.legend()

plt.show()
