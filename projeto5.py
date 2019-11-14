  
#%matplotlib inline
import math 
import random       #para gerar um numero aleatorio na medida
from matplotlib import pyplot as plt    #para plotar o grafico da medida
import numpy as np

#seta o vetor psi em um estado de base dado na entrada
def seta_base(Psi, estado):
    count=0
    Psi[estado]=1
    for c in Psi:    #passa o loop pelo vetor Psi 
        if count != estado:
            Psi[count]=0
        count=count+1   #enumera os componentes de psi
        
#dado um vetor psi de amplitudes de probabilidade de um sitema de bits quanticos com numero arbitrario de bits. Simula uma medida colapsando psi para um dos seus componentes, seguindo a distribuicao de probabilidade adequada.
def medir(Psi):
    r=random.random()   #r recebe um numero aleatorio float de 0 ate 1
    q=0
    count=0
    for c in Psi:       #esse loop vai somando as probabilidades da medida colapsar para algum componente de psi ate r estar contido no intervalo de 0 ate q
        q=q+abs(c)**2   #q recebe o valor absoluto da amplitude de probabilidade desse componente de psi
        if q >r:
            break
        count=count+1   #enumera os componentes de psi
    return count



#mede a funcao de onda Psi n vezes e imprime o grafico com a frequencia de obtencao de cada resultado de medida
#qubits e o numero de qubits
def medir_n(Psi, n, qubits):
    medidas=[0 for c in Psi]
    i=0
    while i<n:
        i=i+1
        resultado=medir(Psi)
        medidas[resultado]=medidas[resultado]+1

    plt.xlabel('Resultado')
    plt.ylabel('Numero de medidas')
    plt.title('Medidas')
    
    index=np.arange(2**qubits) #index de numeros inteiros com numero de resultados possiveis
    rotulo=[bin(c)[2:].zfill(qubits) for c in index] #cria um vetor com rotulos em binario para o eixo x
    plt.xticks(index, rotulo)   #para o grafico usar os rotulos
    plt.bar(np.arange(2**qubits), medidas) #plota o grafico



from sympy import I, Matrix, symbols            #para os produtos de tensores
from sympy.physics.quantum import TensorProduct #para operar gates em registradores com varios bits
import math
import cmath
from sympy.matrices import *

def tensorNqubits(portao, alvo, qubits): #muda o portao de 1 bit para um de N qubits no alvo a escolha
    i=0
    alvo=alvo-1
    while i < alvo: #faz produtos tensoriais para os qubits anteriores ao alvo
        portao=TensorProduct(eye(2), portao) 
        i = i + 1
    i=0
    while i < qubits - alvo - 1: #faz produtos tensoriais para os qubits posteriores ao alvo
        portao=TensorProduct(portao, eye(2)) 
        i = i + 1
    return portao 

#aplica o operador de hadamard no qubit alvo de Psi
def hadamard(Psi, alvo, qubits):
    
    h1q=ones(2) #faz uma matriz 2x2 portao de hadamard para 1 qubit
    h1q[1,1]=-1
    
    had=tensorNqubits(h1q, alvo, qubits)
    
    psi=np.matrix(Psi)  #para poder multiplicar o vetor de estado pelo operador
    
    a=1/math.sqrt(2)
    psi=psi*had*a         #opera a multiplicacao de matrizes mais o fator a que normaliza a matriz do portao de hadamard
    
    return psi.tolist()[0] #transforma a matriz de volta em lista



#operador de difusao de Groover
def bloco_groover (psi, qubits, resposta):
    
    #oraculo##################################################
    psi[resposta]=-psi[resposta]  #muda o sinal da amplitude da resposta correta
    
    #operador de difusao de groover###########################
    #loop para fazer hadamard em todos os qubits
    count=1
    while count <= qubits:    
        psi=hadamard(psi, count, qubits) #hadamard em todos os qubits
        count = count + 1
        
    #Operador J (matriz identidade com primeiro elemento da diagonal -1)
    psi[0]=-psi[0] 
    
    #hadamard em todos os qubits de novo
    count=1
    while count <= qubits:
        psi=hadamard(psi, count, qubits)
        count = count + 1
    return psi



######################################################################

qubits=input()
resposta=1
Psi=[0]*2**qubits
seta_base(Psi,0)
count=1
while count <= qubits:    #passa o loop pelo vetor psi
        Psi=hadamard(Psi, count, qubits) #hadamard em todos os qubits
        count = count + 1

#blocos do operador de difusao de groover
count=1
while count <= round(math.pi/4*math.sqrt(2**qubits)):    #passa o loop pelo vetor psi
        Psi=bloco_groover(Psi, qubits, resposta)
        count = count + 1

print("done")
#medir_n(Psi, 10, qubits)
