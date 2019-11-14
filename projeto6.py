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



######################################################################

qubits=input()
resposta=1

Psi=[0]*2**qubits
Psi[0]=1

count=1
while count <= qubits:    #passa o loop pelo vetor psi
        Psi=sp_hadamard(Psi, count, qubits) #hadamard em todos os qubits
        count = count + 1

#blocos do operador de difusao de groover
count=1
while count <= round(math.pi/4*math.sqrt(2**qubits)):    #passa o loop pelo vetor psi
        Psi=sp_bloco_groover(Psi, qubits, resposta)
        count = count + 1
print(Psi)
print("done")
