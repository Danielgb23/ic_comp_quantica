%matplotlib inline
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

#retorna o estado Psi com o portao de mundanca de fase aplicado no qubit alvo
def muda_fase(Psi, alvo, fase, qubits):
    
    p1q=eye(2) #faz uma matriz 2x2 portao de hadamard para 1 qubit
    p1q[1,1]=cmath.exp( complex(0, fase) )
    
    phase=tensorNqubits(p1q , alvo, qubits)
    
    psi=np.matrix(Psi)  #para poder multiplicar o vetor de estado pelo operador
    
    psi=psi*phase         #opera a multiplicacao de matrizes mais o fator a que normaliza a matriz do portao de hadamard
    
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

#cnot de 3 quibts com o qubit 2 como controlador
def cnot3(psi, qubit_negado):
    #matriz e produto tensorial se o bit negado for o 1
    if qubit_negado == 1:
        cnotm2= Matrix([[1, 0, 0, 0], #matriz de 2 qubits onde primeiro controla o posterior
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 0]])
        cnot=TensorProduct(cnotm2, eye(2))#produto tensorial para operar 3 qubits
        
    #matriz e produto tensorial se o bit negado for o 3
    else:
        cnotm2= Matrix([[1, 0, 0, 0], #matriz de 2 qubits onde um qubit controla o anterior
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]])
        cnot=TensorProduct(eye(2), cnotm2)#produto tensorial para operar 3 qubits
    
    psi=np.matrix(psi)  #para poder multiplicar o vetor de estado pelo operador
    
    psi=psi*cnot       #opera a multiplicacao de matrizes
    
    return psi.tolist()[0] #transforma a matriz de volta em lista
    
import math 
import cmath
import scipy.sparse as scp


def sp_tensorNqubits(portao, alvo, qubits): #muda o portao de 1 bit para um de N qubits no alvo a escolha
    i=0
    alvo=alvo-1
    while i < alvo: #faz produtos tensoriais para os qubits anteriores ao alvo
        portao=scp.kron(scp.eye(2), portao) 
        i = i + 1
    i=0
    while i < qubits - alvo - 1: #faz produtos tensoriais para os qubits posteriores ao alvo
        portao=scp.kron(portao, scp.eye(2)) 
        i = i + 1
    return portao 

#aplica o operador de hadamard no qubit alvo de Psi
def sp_hadamard(Psi, alvo, qubits):
    
    h1q=scp.csr_matrix( #faz uma matriz 2x2 portao de hadamard para 1 qubit
   [[1, 1],
     [1, -1]])

    had=sp_tensorNqubits(h1q, alvo, qubits)
    
    psi=scp.csr_matrix(Psi)  #para poder multiplicar o vetor de estado pelo operador
    
    a=1/math.sqrt(2)
    psi=psi*had*a         #opera a multiplicacao de matrizes mais o fator a que normaliza a matriz do portao de hadamard
    return psi #transforma a matriz de volta em lista

#retorna o estado Psi com o portao de mundanca de fase aplicado no qubit alvo



#operador de difusao de Groover
def sp_bloco_groover (psi, qubits, resposta):
    
    #oraculo##################################################
    psi[0, resposta]=-psi[0, resposta]  #muda o sinal da amplitude da resposta correta
    
    #operador de difusao de groover###########################
    #loop para fazer hadamard em todos os qubits
    count=1
    while count <= qubits:    
        psi=sp_hadamard(psi, count, qubits) #hadamard em todos os qubits
        count = count + 1
        
    #Operador J (matriz identidade com primeiro elemento da diagonal -1)
    psi[0,0]=-psi[0, 0] 
    
    #hadamard em todos os qubits de novo
    count=1
    while count <= qubits:
        psi=sp_hadamard(psi, count, qubits)
        count = count + 1
    return psi



