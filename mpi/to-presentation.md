comentários:
- de maneira simples, pegamos o que foi desenvolvido para o t01 de threads e passamos para mpi
- implementação em mestre-escravo onde o mestre também realiza computação
- mudanças realizadas:
    - vars que precisam ser enviadas entre os processos:
        - centroids_tmp (global)
            - só serve de "ligação", cada processo tem o centroids, o centroids_tmp é a "união" de todos os eles
        - map_tmp (global)
            - uma var especial utilizada somente pelo mestre para receber o map de cada escravo
        - dirty_tmp (global)
            - só serve de "ligação", cada processo tem o dirty, o dirty_tmp é a "união" de todos os eles
        - too_far_tmp (global)
            - só serve de "ligação", cada processo tem o too_far, o too_far_tmp é a "união" de todos os eles
        - population_tmp
            - só serve de "ligação", cada processo tem o population, o population_tmp é a "união" de todos os eles
    - size, rank, baseCalc e finalCalc também são novas vars globais
    - population e compute_centroids possuem a computação dividida para cada processo com base nos pontos (cálculo do baseCalc e finalCalc é feito no main)
    - cálculo de base e de final não pode ter precedência explícita com () para conseguir realizar arredondamentos corretos e assim pegar toda a faixa de valores
    - quando necessário se aloca e libera dinamicamente memória para as novas vars que são ponteiros
    - sempre que vai ser dado um exit, antes se finaliza o ambiente mpi
    - após cada chamada do populate:
        - todos os processos são barrados
        - todos os processos ficam com o mesmo too_far e com o mesmo dirty
    - para cada partição no maior for do compute_centroids:
        - todos os processos são barrados
        - todos os processos ficam com o mesmo population e o mesmo centroids antes do cálculo do centro da partição
        - todos os processos calculam o centro da partição que possuirá o mesmo valor final para todos os cálculos
    - map de todo mundo é atualizado após a execução do do {...} while da seguinte forma:
        1. todos escravos enviam o map para o mestre
        2. mestre faz um map com a informações de todos os maps recebidos
        3. mestre envia para todos os escravos o map final através de um broadcast
- possíveis melhorias:
    - testar possibilidade de diminuir o número de vars novas e verificar a real necessidade delas serem globais ou não
    - para deixar o código mais "bonito" continuar utilizando o too_far como variável de controle no do {...} while. De maneira similar, atualizar o dirty de cada processo com o valor final de dirty_tmp
    - verificar melhoria no desempenho se somente um processo realizasse o cálculo do centro da partição e enviasse o resultado para todos os outros processos (estilo map)
    - pensando de maneira mais genérica, o ideal seria fazer somente uma barreira no do {...} while que ficaria após o compute_centroids. Assim,teria que se pensar em um modo de todos terem o mesmo dirty no começo do compute_centroids
    - também pensando de maneira mais genérica, seria melhor se não houvesse uma barreira para cada partição no compute_centroids, o ideal seria somente uma barreira no final da função, mas isso impactaria no cálculo do centro da partição, que talvez tivesse que ser desmembrado de lá e colocado em outro lugar
