# First research internship at LFDQ (quantum devices physics laboratory) UNICAMP

[The material of this internship was made in portuguese.]

[3 minutes video in portuguese](https://youtu.be/jOJm0fjOHDE?si=EniDWa5tUfLl2Tvt)

On this internship I made several python simulations of a quantum computer from these 
[project suggestions by D. Candela](https://github.com/Danielgb23/ic_comp_quantica/blob/master/Candela%20-%202015%20-%20Undergraduate%20computational%20physics%20projects%20on%20qu.pdf.pdf).

My simulations are [in this notebook](https://github.com/Danielgb23/ic_comp_quantica/blob/master/Projeto_de_computacao_quantica.ipynb) where
I went from the basic of simulating one or several qubits with a vector of complex coefficients, as well as the gates with matrices, to more complicated simulations
like the Grover search algorithm and even the complete Shor algorithm with the quantum and classical parts.
I also used sparse matrices to test the computing limits of my machine where I did increasingly bigger runs of the Grover algorithm.

![Shor's algorithm period finding result](https://github.com/user-attachments/assets/6eb6ef9d-a802-4c78-9705-c1a0425dc271)

_Shor's algorithm period finding step (quantum step) result for the factorization o 21 with a=10 (with a resulting period of 6)._


[In the second notebook](https://github.com/Danielgb23/ic_comp_quantica/blob/master/caderno_qiskit_comp_real.ipynb) I did the Grover search algorithm and the Shor factorization algorithm
in the IBM's real quantum device on the cloud using qiskit.

![image](https://github.com/user-attachments/assets/0938b5f8-00c2-4714-8770-0025803a1cfb)

_f(x) gate of the Shor's algorithm made with qiskit gates._
